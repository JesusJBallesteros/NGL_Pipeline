"""Physically-constrained cleaning of DeepLabCut pose.

Body-part LABELS are not hardcoded: the reader discovers every part in the csv and
the caller chooses any subset. The algorithms need anatomical MEANING, so labels are
mapped to ROLES (beak/head/back/left_wing/right_wing/tail). By default a role maps to
the identically-named label; use cfg['roles'] to map different names, e.g.
    roles = {'beak': 'bill', 'head': 'nape'}
Any selected label without a role is still cleaned and returned - it just gets generic
treatment (temporal interpolation, no rigid-bone constraint).

A joint in a frame is treated as IMPLAUSIBLE (re-estimated from surrounding frames) when:
  1. likelihood < p_cut                          (DLC uncertainty)
  2. local deviation > dev_fac x body-length     (jump / noise spike)
  3. a rigid bone length is grossly off          (size is invariant)
  4. head is posterior to the wing line          (anatomical ordering)
Rules 3-4 only apply to the roles actually present.

Public API:
    available_parts(csv_path)          -> list of every label in the file
    read_dlc(csv_path, parts=None)     -> XY, P, n, available
    clean(csv_path, cfg)               -> dict(X, Y, A, info)
"""
import csv
import os
import numpy as np

# anatomical template, expressed in ROLES (not labels)
ROLES = ['beak', 'head', 'back', 'left_wing', 'right_wing', 'tail']
RIGID_ROLES = [('beak', 'head'), ('head', 'back'), ('back', 'tail'),
               ('back', 'left_wing'), ('back', 'right_wing'), ('left_wing', 'right_wing')]
PARENT_ROLES = {'beak': 'head', 'tail': 'back', 'left_wing': 'back', 'right_wing': 'back'}
DEGREE_ROLES = {'beak': 1, 'tail': 1, 'left_wing': 2, 'right_wing': 2, 'head': 2, 'back': 4}

DEFAULTS = dict(
    parts=None,            # labels to use; None -> every part in the csv
    roles=None,            # {role: label} overrides; default identity by name
    body_px=None,          # explicit body scale if no rigid bone can supply one
    fps_in=59.94,
    downsample_step=2,
    target_fps=None,
    start_time=None,
    end_time=None,
    p_cut=0.5,
    w_med=9,
    dev_fac=0.6,
    smooth=5,
    bone_tol_frac=0.4,
    bone_tol_mad=5.0,
    order_margin=0.10,
)


# ------------------------------------------------------------------ I/O
def _isnum(s):
    try:
        float(s); return True
    except (ValueError, TypeError):
        return False


def _header_rows(rows, path=''):
    """Locate the DLC 'coords' header row (cells are x/y/likelihood); the
    'bodyparts' row is the one directly above it. Name-agnostic."""
    for i, r in enumerate(rows[:10]):
        vals = [c.strip().lower() for c in r[1:] if c.strip()]
        if vals and set(vals) <= {'x', 'y', 'likelihood'}:
            if i == 0:
                break
            return i - 1, i
    raise ValueError('could not find the DLC "coords"/"bodyparts" header rows in %s' % path)


def _colmap(rows, bi, ci):
    """{label: {'x': col, 'y': col, 'likelihood': col}} plus label order."""
    cols, order = {}, []
    for j, (name, coord) in enumerate(zip(rows[bi], rows[ci])):
        c = coord.strip().lower()
        n = name.strip()
        if c in ('x', 'y', 'likelihood') and n:
            if n not in cols:
                cols[n] = {}; order.append(n)
            cols[n][c] = j
    return cols, [n for n in order if len(cols[n]) == 3]


def available_parts(csv_path):
    rows = list(csv.reader(open(csv_path, newline='')))
    bi, ci = _header_rows(rows, csv_path)
    return _colmap(rows, bi, ci)[1]


def read_dlc(csv_path, parts=None):
    rows = list(csv.reader(open(csv_path, newline='')))
    bi, ci = _header_rows(rows, csv_path)
    cols, avail = _colmap(rows, bi, ci)

    use = [str(p).strip() for p in parts] if parts else list(avail)
    missing = [p for p in use if p not in cols or len(cols[p]) < 3]
    if missing:
        raise ValueError('body part(s) %s not found in %s.\nAvailable: %s'
                         % (missing, os.path.basename(csv_path), avail))

    ncol = len(rows[ci])
    data = np.array([[float(x) for x in r] for r in rows[ci + 1:]
                     if len(r) == ncol and _isnum(r[-1])], dtype=float)
    XY = {p: data[:, (cols[p]['x'], cols[p]['y'])].copy() for p in use}
    P = {p: data[:, cols[p]['likelihood']].copy() for p in use}
    return XY, P, len(data), avail


# ------------------------------------------------------------------ helpers
def _smooth1d(a, w):
    if w <= 1:
        return a
    p = w // 2
    return np.convolve(np.pad(a, (p, p), 'edge'), np.ones(w)/w, 'valid')[:len(a)]


def _roll_median(xy, valid, w):
    n, pad = len(xy), w // 2
    out = xy.copy()
    for i in range(n):
        lo, hi = max(0, i-pad), min(n, i+pad+1)
        sel = valid[lo:hi]
        seg = xy[lo:hi][sel] if sel.any() else xy[lo:hi]
        out[i] = np.median(seg, axis=0)
    return out


def _fill(xy, valid, smooth):
    idx = np.arange(len(xy))
    out = xy.astype(float).copy()
    if valid.sum() >= 2:
        out[:, 0] = np.interp(idx, idx[valid], xy[valid, 0])
        out[:, 1] = np.interp(idx, idx[valid], xy[valid, 1])
    elif valid.sum() == 1:
        out[:] = xy[valid][0]
    if smooth > 1:
        out[:, 0] = _smooth1d(out[:, 0], smooth)
        out[:, 1] = _smooth1d(out[:, 1], smooth)
    return out


def _dev(xy, valid, w):
    return np.hypot(*(xy - _roll_median(xy, valid, w)).T)


def _polar_reconstruct(J, validJ, anchor, Lc, smooth):
    """Interpolate DIRECTION from surrounding frames, keep observed length where
    valid, snap to canonical length Lc where flagged -> bone length invariant."""
    n = len(J); idx = np.arange(n)
    rel = J - anchor
    ang = np.arctan2(rel[:, 1], rel[:, 0]); rad = np.hypot(rel[:, 0], rel[:, 1])
    v = validJ
    if v.sum() >= 2:
        ang_f = np.interp(idx, idx[v], np.unwrap(ang[v]))
        rad_f = np.interp(idx, idx[v], rad[v])
    elif v.sum() == 1:
        ang_f = np.full(n, ang[v][0]); rad_f = np.full(n, rad[v][0])
    else:
        ang_f = np.zeros(n); rad_f = np.full(n, Lc)
    rad_f = rad_f.copy(); rad_f[~v] = Lc
    ang_s = _smooth1d(ang_f, smooth); rad_s = _smooth1d(rad_f, smooth)
    return anchor[:, 0] + rad_s*np.cos(ang_s), anchor[:, 1] + rad_s*np.sin(ang_s)


def _resolve_roles(parts, roles_cfg):
    """{role: label}; a role is present only if its label is among `parts`."""
    out = {}
    given = {str(k).strip().lower(): str(v).strip() for k, v in (roles_cfg or {}).items() if v}
    for r in ROLES:
        lab = given.get(r, r if r in parts else None)
        if lab is None:
            continue
        if lab not in parts:
            raise ValueError('role %r maps to label %r, which is not among the selected '
                             'parts %s' % (r, lab, parts))
        out[r] = lab
    return out


def _canon_len(canon, a, b):
    return canon.get((a, b), canon.get((b, a)))


def _usable_hint(P_raw, labels, p_cut, fps):
    """Where in the WHOLE recording does tracking actually pass the gate?"""
    cnt = np.sum([P_raw[k] >= p_cut for k in labels], axis=0)
    good = cnt >= max(1, (len(labels) + 1) // 2)      # at least half the parts
    if not good.any():
        return ('no part of the whole recording passes this gate either - check the '
                'DLC model or lower pCut')
    idx = np.where(good)[0]
    return ('tracking looks usable from t=%.0fs to t=%.0fs (%.0f%% of the recording passes)'
            % (idx[0]/fps, idx[-1]/fps, 100.0*good.mean()))


def _slice_and_downsample(XY, P, cfg, labels):
    """Returns the windowed/downsampled data plus the ABSOLUTE raw frame index of
    every kept sample (needed for timestamps: the downsampler picks the better
    frame of each pair, so the spacing is not perfectly uniform)."""
    n = len(XY[labels[0]])
    fps = cfg['fps_in']
    s = 0 if cfg['start_time'] is None else max(0, int(round(cfg['start_time']*fps)))
    e = n if cfg['end_time'] is None else min(n, int(round(cfg['end_time']*fps)))
    idx = np.arange(s, e)
    XY = {k: v[s:e] for k, v in XY.items()}
    P = {k: v[s:e] for k, v in P.items()}
    step = cfg['downsample_step']
    if cfg.get('target_fps'):
        step = max(1, int(round(fps / cfg['target_fps'])))
    if step > 1:
        lik = np.sum([P[k] for k in labels], axis=0)
        keep = [max(range(a, min(a+step, len(lik))), key=lambda j: lik[j])
                for a in range(0, len(lik), step)]
        XY = {k: v[keep] for k, v in XY.items()}
        P = {k: v[keep] for k, v in P.items()}
        idx = idx[keep]
    return XY, P, step, idx


# ------------------------------------------------------------------ main
def clean(csv_path, cfg=None):
    cfg = {**DEFAULTS, **(cfg or {})}

    def log(m):
        if cfg.get('verbose', True):
            print('[clean] ' + m, flush=True)

    XY, P, n_raw, avail = read_dlc(csv_path, cfg.get('parts'))
    labels = list(XY.keys())
    log('read %d frames from %s' % (n_raw, os.path.basename(csv_path)))
    log('parts in file: %s' % ', '.join(avail))
    log('parts used   : %s' % ', '.join(labels))

    role2lab = _resolve_roles(labels, cfg.get('roles'))
    lab2role = {v: k for k, v in role2lab.items()}
    have = lambda r: r in role2lab
    L = lambda r: role2lab[r]
    unroled = [p for p in labels if p not in lab2role]
    log('roles: ' + (', '.join('%s->%s' % (r, role2lab[r]) for r in ROLES if have(r)) or 'none')
        + (' | no role (generic): %s' % ', '.join(unroled) if unroled else ''))

    P_raw = P                     # full-length likelihoods, before windowing
    XY, P, step, frame_index = _slice_and_downsample(XY, P, cfg, labels)
    n = len(XY[labels[0]])
    log('window + downsample (step %d) -> %d frames @ %.2f fps' % (step, n, cfg['fps_in']/step))

    P_CUT, W_MED, DEV_FAC, SMOOTH = cfg['p_cut'], cfg['w_med'], cfg['dev_fac'], cfg['smooth']
    rigid = [(L(a), L(b)) for a, b in RIGID_ROLES if have(a) and have(b)]
    degree = {p: DEGREE_ROLES.get(lab2role.get(p), 1) for p in labels}

    valid = {k: P[k] >= P_CUT for k in labels}                          # rule 1
    reason = {k: np.array(['lik' if not v else '' for v in valid[k]], dtype=object)
              for k in labels}

    # Refuse to invent a skeleton out of untracked noise: without a single frame above
    # the gate there is nothing to interpolate from, and we would silently draw raw
    # detections scattered all over the arena.
    nval = {k: int(valid[k].sum()) for k in labels}
    log('usable frames (likelihood >= %.2f): %s' % (P_CUT, ', '.join(
        '%s=%d/%d (%.0f%%)' % (k, nval[k], n, 100.0*nval[k]/max(n, 1)) for k in labels)))
    win = ('t=%s..%s s' % (cfg['start_time'] if cfg['start_time'] is not None else 'start',
                           cfg['end_time'] if cfg['end_time'] is not None else 'end'))
    if max(nval.values()) == 0:
        best = max(float(P[k].max()) for k in labels)
        raise ValueError(
            'No frame in the selected window (%s) passes the likelihood gate '
            '(p_cut=%.2f): the best likelihood found was %.3f, so the animal is '
            'effectively untracked there.\n  hint: %s\n'
            'Fix: pick a different startTime/endTime, or lower pCut.'
            % (win, P_CUT, best, _usable_hint(P_raw, labels, P_CUT, cfg['fps_in'])))
    dead = [k for k in labels if nval[k] == 0]
    if dead:
        log('WARNING: no usable frames for %s in %s - these parts will be left blank '
            'rather than drawn from raw noise' % (', '.join(dead), win))

    # body scale: prefer head-back, else the first available rigid bone, else cfg
    L_body = None
    if have('head') and have('back'):
        d = np.hypot(*(XY[L('head')] - XY[L('back')]).T)
        m = valid[L('head')] & valid[L('back')]
        if m.any():
            L_body = float(np.median(d[m]))
    if L_body is None and rigid:
        a, b = rigid[0]
        d = np.hypot(*(XY[a] - XY[b]).T)
        m = valid[a] & valid[b]
        if m.any():
            L_body = float(np.median(d[m]))
    if L_body is None or not np.isfinite(L_body) or L_body <= 0:
        L_body = float(cfg.get('body_px') or 30.0)
        log('WARNING: no rigid bone available for the body scale; using %.1f px '
            '(set cfg["body_px"] to control this)' % L_body)
    T_DEV = DEV_FAC * L_body
    log('detecting implausible joints (body=%.1fpx, jump>%.1fpx)' % (L_body, T_DEV))

    pos_sigma = {}
    for k in labels:                                                   # rule 2: jumps
        d = _dev(XY[k], valid[k], W_MED)
        jump = valid[k] & (d > T_DEV)
        for i in np.where(jump)[0]:
            reason[k][i] = 'jump'
        valid[k] &= ~jump
        # positional scatter of the good detections about their own smooth path,
        # i.e. a robust 1-sigma tracking noise in px (MAD -> sigma)
        good = d[valid[k]]
        pos_sigma[k] = float(1.4826 * np.median(good)) if good.size else float('nan')

    est = {k: _fill(XY[k], valid[k], SMOOTH) for k in labels}
    devc = {k: _dev(est[k], np.ones(n, bool), W_MED) for k in labels}

    def attribute(cands, i):
        return max(cands, key=lambda k: (round(devc[k][i], 1), 1 - P[k][i], -degree[k]))

    canon = {}
    for a, b in rigid:                                                 # rule 3: bone length
        both = valid[a] & valid[b]
        if both.sum() < 3:
            continue
        d = np.hypot(*(XY[a][both] - XY[b][both]).T)
        med = np.median(d); mad = np.median(np.abs(d - med)); canon[(a, b)] = float(med)
        tol = max(cfg['bone_tol_frac']*med, cfg['bone_tol_mad']*mad)
        Lb = np.hypot(*(est[a] - est[b]).T)
        for i in np.where(both & (np.abs(Lb - med) > tol))[0]:
            c = attribute((a, b), i)
            if valid[c][i]:
                valid[c][i] = False; reason[c][i] = 'bone'

    ordering_ok = all(have(r) for r in ('head', 'back', 'left_wing', 'right_wing'))
    if ordering_ok:                                                    # rule 4: ordering
        u = est[L('head')] - est[L('back')]
        un = np.linalg.norm(u, axis=1, keepdims=True)
        u = np.divide(u, un, out=np.zeros_like(u), where=un > 1e-6)
        Wmid = 0.5*(est[L('left_wing')] + est[L('right_wing')])
        proj_h = np.sum((est[L('head')]-est[L('back')])*u, axis=1)
        proj_w = np.sum((Wmid-est[L('back')])*u, axis=1)
        for i in np.where(proj_w > proj_h + cfg['order_margin']*L_body)[0]:
            c = attribute((L('head'), L('left_wing'), L('right_wing')), i)
            if valid[c][i]:
                valid[c][i] = False; reason[c][i] = 'order'
    else:
        log('ordering rule skipped (needs head, back and both wings)')

    log('flagged frames/joint: ' + ', '.join('%s=%d' % (k, int((~valid[k]).sum())) for k in labels))

    # ---- reconstruction: polar where a rigid parent exists, else temporal ----
    parent = {L(r): L(pr) for r, pr in PARENT_ROLES.items() if have(r) and have(pr)}
    X = {}; Y = {}
    for k in labels:
        if k not in parent:
            f = _fill(XY[k], valid[k], SMOOTH)
            X[k], Y[k] = f[:, 0], f[:, 1]
    for k, pk in parent.items():
        anchor = np.column_stack([X[pk], Y[pk]])
        Lc = _canon_len(canon, k, pk)
        if Lc is None:
            f = _fill(XY[k], valid[k], SMOOTH)
            X[k], Y[k] = f[:, 0], f[:, 1]
        else:
            X[k], Y[k] = _polar_reconstruct(XY[k], valid[k], anchor, Lc, SMOOTH)

    if ordering_ok:            # anatomical clamp: wing line never ahead of the head
        hd, bk, lw, rw = L('head'), L('back'), L('left_wing'), L('right_wing')
        ux = X[hd]-X[bk]; uy = Y[hd]-Y[bk]
        Ln = np.hypot(ux, uy); ok = Ln > 1e-6
        ux = np.where(ok, ux/np.where(ok, Ln, 1), 0.0); uy = np.where(ok, uy/np.where(ok, Ln, 1), 0.0)
        wmx = 0.5*(X[lw]+X[rw]); wmy = 0.5*(Y[lw]+Y[rw])
        ph = (X[hd]-X[bk])*ux + (Y[hd]-Y[bk])*uy
        pw = (wmx-X[bk])*ux + (wmy-Y[bk])*uy
        margin = cfg['order_margin']*L_body
        delta = np.where(ok & (pw > ph - margin), pw - (ph - margin), 0.0)
        for k in (lw, rw):
            X[k] = X[k] - delta*ux; Y[k] = Y[k] - delta*uy
    log('reconstruction + anatomical clamp complete')

    for k in dead:            # NaN keeps matplotlib from drawing markers or bones
        X[k] = np.full(n, np.nan); Y[k] = np.full(n, np.nan)

    ALPHA = {}; FLAG = {}
    for k in labels:
        FLAG[k] = ~valid[k]
        a = np.clip(P[k], 0.12, 0.98)
        a[FLAG[k]] = np.minimum(a[FLAG[k]], 0.22)
        ALPHA[k] = a

    info = dict(n=n, L_body=L_body, T_DEV=T_DEV, canon=canon, step=step,
                out_fps=cfg['fps_in']/step, FLAG=FLAG, reason=reason,
                order=labels, roles=role2lab, available=avail,
                P=P, frame_index=frame_index, fps_in=cfg['fps_in'],
                pos_sigma=pos_sigma, dead=dead)
    return dict(X=X, Y=Y, A=ALPHA, info=info)


if __name__ == '__main__':
    import sys
    if len(sys.argv) > 2 and sys.argv[2] == '--list':
        print('parts:', ', '.join(available_parts(sys.argv[1])))
    else:
        res = clean(sys.argv[1])
        i = res['info']
        print('frames=%d  out_fps=%.2f  body=%.1fpx' % (i['n'], i['out_fps'], i['L_body']))
