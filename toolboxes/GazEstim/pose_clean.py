"""Physically-constrained cleaning of DeepLabCut pose for the jackdaw skeleton.

Parameterized, importable version of the development `clean_pose.py`. A joint in a
frame is treated as IMPLAUSIBLE (re-estimated from surrounding frames) when any of:
  1. likelihood < p_cut                               (DLC uncertainty)
  2. local deviation > dev_fac x body-length          (jump / noise spike)
  3. a rigid bone length is grossly off               (size is invariant)
  4. head is posterior to the wing line               (anatomical ordering)
Reconstruction is by polar anchoring (wings & tail -> back, beak -> head):
interpolate DIRECTION from surrounding frames, snap LENGTH to canonical.

Public API:
    XY, P, n_raw = read_dlc(csv_path)
    result = clean(csv_path, cfg)      # cfg: dict of overrides (see DEFAULTS)
    result -> dict(X, Y, A, info)      # info has n, canon, L_body, FLAG, reason ...
"""
import csv
import os
import numpy as np

REQUIRED = ['beak', 'tail', 'left_wing', 'right_wing', 'head', 'back']
DEGREE   = {'beak': 1, 'tail': 1, 'left_wing': 2, 'right_wing': 2, 'head': 2, 'back': 4}
RIGID    = [('beak', 'head'), ('head', 'back'), ('back', 'tail'),
            ('back', 'left_wing'), ('back', 'right_wing'), ('left_wing', 'right_wing')]

DEFAULTS = dict(
    fps_in=59.94,          # raw video fps
    downsample_step=2,     # keep 1 frame per window (higher total likelihood) -> fps_in/step
    target_fps=None,       # if set, overrides downsample_step = round(fps_in/target_fps)
    start_time=None,       # seconds; None -> from first frame
    end_time=None,         # seconds; None -> to last frame
    p_cut=0.5,             # likelihood gate
    w_med=9,               # rolling-median window for the smooth reference
    dev_fac=0.6,           # jump threshold = dev_fac x body length
    smooth=5,              # final temporal smoothing window
    bone_tol_frac=0.4,     # bone violation if |len-canon| > max(frac*canon, mad*MAD)
    bone_tol_mad=5.0,
    order_margin=0.10,     # head-behind-wing clamp margin (x body length)
)


# ------------------------------------------------------------------ I/O
def read_dlc(csv_path):
    """Read a DLC .csv robustly: auto-detect header rows, an optional leading
    frame-index column, and map each required part to its (x,y,likelihood) columns
    by name from the `bodyparts` header row (handles any column order)."""
    rows = list(csv.reader(open(csv_path, newline='')))
    bodyparts = None
    for r in rows[:8]:
        if any(c.strip() in REQUIRED for c in r):
            bodyparts = [c.strip() for c in r]
            break
    if bodyparts is None:
        raise ValueError('Could not find a bodyparts header row containing %s' % REQUIRED)

    def is_data(r):
        try:
            float(r[-1]); float(r[-2]); return len(r) == len(bodyparts)
        except (ValueError, IndexError):
            return False
    data = np.array([[float(x) for x in r] for r in rows if is_data(r)], dtype=float)

    cols = {}
    for part in REQUIRED:
        idx = [j for j, name in enumerate(bodyparts) if name == part]
        if len(idx) < 3:
            raise ValueError('part "%s" not found 3x (x,y,likelihood) in header' % part)
        cols[part] = (idx[0], idx[1], idx[2])
    XY = {p: data[:, (cols[p][0], cols[p][1])].copy() for p in REQUIRED}
    P  = {p: data[:, cols[p][2]].copy() for p in REQUIRED}
    return XY, P, len(data)


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


def _slice_and_downsample(XY, P, cfg):
    n = len(XY['beak'])
    fps = cfg['fps_in']
    s = 0 if cfg['start_time'] is None else max(0, int(round(cfg['start_time']*fps)))
    e = n if cfg['end_time'] is None else min(n, int(round(cfg['end_time']*fps)))
    XY = {k: v[s:e] for k, v in XY.items()}
    P = {k: v[s:e] for k, v in P.items()}
    step = cfg['downsample_step']
    if cfg.get('target_fps'):
        step = max(1, int(round(fps / cfg['target_fps'])))
    if step > 1:
        lik = np.sum([P[k] for k in REQUIRED], axis=0)
        keep = [max(range(a, min(a+step, len(lik))), key=lambda j: lik[j])
                for a in range(0, len(lik), step)]
        XY = {k: v[keep] for k, v in XY.items()}
        P = {k: v[keep] for k, v in P.items()}
    return XY, P, step


# ------------------------------------------------------------------ main
def clean(csv_path, cfg=None):
    cfg = {**DEFAULTS, **(cfg or {})}

    def log(m):
        if cfg.get('verbose', True):
            print('[clean] ' + m, flush=True)

    XY, P, n_raw = read_dlc(csv_path)
    log('read %d frames from %s' % (n_raw, os.path.basename(csv_path)))
    XY, P, step = _slice_and_downsample(XY, P, cfg)
    n = len(XY['beak'])
    log('window + downsample (step %d) -> %d frames @ %.2f fps' % (step, n, cfg['fps_in']/step))
    P_CUT, W_MED, DEV_FAC, SMOOTH = cfg['p_cut'], cfg['w_med'], cfg['dev_fac'], cfg['smooth']

    valid = {k: P[k] >= P_CUT for k in REQUIRED}                       # rule 1
    reason = {k: np.array(['lik' if not v else '' for v in valid[k]], dtype=object)
              for k in REQUIRED}

    hb = np.hypot(*(XY['head'] - XY['back']).T)
    L_body = float(np.median(hb[valid['head'] & valid['back']]))
    T_DEV = DEV_FAC * L_body
    log('detecting implausible joints (body=%.1fpx, jump>%.1fpx)' % (L_body, T_DEV))

    devc = {}
    for k in REQUIRED:                                                 # rule 2: jumps
        d = _dev(XY[k], valid[k], W_MED)
        jump = valid[k] & (d > T_DEV)
        for i in np.where(jump)[0]:
            reason[k][i] = 'jump'
        valid[k] &= ~jump

    est = {k: _fill(XY[k], valid[k], SMOOTH) for k in REQUIRED}
    devc = {k: _dev(est[k], np.ones(n, bool), W_MED) for k in REQUIRED}

    def attribute(cands, i):
        return max(cands, key=lambda k: (round(devc[k][i], 1), 1 - P[k][i], -DEGREE[k]))

    canon = {}
    for a, b in RIGID:                                                 # rule 3: bone length
        both = valid[a] & valid[b]
        d = np.hypot(*(XY[a][both] - XY[b][both]).T)
        med = np.median(d); mad = np.median(np.abs(d - med)); canon[(a, b)] = float(med)
        tol = max(cfg['bone_tol_frac']*med, cfg['bone_tol_mad']*mad)
        L = np.hypot(*(est[a] - est[b]).T)
        for i in np.where(both & (np.abs(L - med) > tol))[0]:
            c = attribute((a, b), i)
            if valid[c][i]:
                valid[c][i] = False; reason[c][i] = 'bone'

    u = est['head'] - est['back']                                     # rule 4: ordering
    un = np.linalg.norm(u, axis=1, keepdims=True)
    u = np.divide(u, un, out=np.zeros_like(u), where=un > 1e-6)
    Wmid = 0.5*(est['left_wing'] + est['right_wing'])
    proj_h = np.sum((est['head']-est['back'])*u, axis=1)
    proj_w = np.sum((Wmid-est['back'])*u, axis=1)
    for i in np.where(proj_w > proj_h + cfg['order_margin']*L_body)[0]:
        c = attribute(('head', 'left_wing', 'right_wing'), i)
        if valid[c][i]:
            valid[c][i] = False; reason[c][i] = 'order'

    log('flagged frames/joint: ' + ', '.join('%s=%d' % (k, int((~valid[k]).sum())) for k in REQUIRED))

    # reconstruction: core anchors world-interp, peripherals polar-anchored
    X = {}; Y = {}
    backW = _fill(XY['back'], valid['back'], SMOOTH)
    headW = _fill(XY['head'], valid['head'], SMOOTH)
    X['back'], Y['back'] = backW[:, 0], backW[:, 1]
    X['head'], Y['head'] = headW[:, 0], headW[:, 1]
    parent = {'beak': headW, 'tail': backW, 'left_wing': backW, 'right_wing': backW}
    Lc = {'beak': canon[('beak', 'head')], 'tail': canon[('back', 'tail')],
          'left_wing': canon[('back', 'left_wing')], 'right_wing': canon[('back', 'right_wing')]}
    for k in ['beak', 'tail', 'left_wing', 'right_wing']:
        X[k], Y[k] = _polar_reconstruct(XY[k], valid[k], parent[k], Lc[k], SMOOTH)

    # anatomical clamp: wing line never ahead of the head
    ux = X['head']-X['back']; uy = Y['head']-Y['back']
    L = np.hypot(ux, uy); ok = L > 1e-6
    ux = np.where(ok, ux/np.where(ok, L, 1), 0.0); uy = np.where(ok, uy/np.where(ok, L, 1), 0.0)
    wmx = 0.5*(X['left_wing']+X['right_wing']); wmy = 0.5*(Y['left_wing']+Y['right_wing'])
    ph = (X['head']-X['back'])*ux + (Y['head']-Y['back'])*uy
    pw = (wmx-X['back'])*ux + (wmy-Y['back'])*uy
    margin = cfg['order_margin']*L_body
    delta = np.where(ok & (pw > ph - margin), pw - (ph - margin), 0.0)
    for k in ('left_wing', 'right_wing'):
        X[k] = X[k] - delta*ux; Y[k] = Y[k] - delta*uy
    log('reconstruction + anatomical clamp complete')

    ALPHA = {}; FLAG = {}
    for k in REQUIRED:
        FLAG[k] = ~valid[k]
        a = np.clip(P[k], 0.12, 0.98)
        a[FLAG[k]] = np.minimum(a[FLAG[k]], 0.22)
        ALPHA[k] = a

    out_fps = cfg['fps_in'] / step
    info = dict(n=n, L_body=L_body, T_DEV=T_DEV, canon=canon, step=step,
                out_fps=out_fps, FLAG=FLAG, reason=reason, order=REQUIRED)
    return dict(X=X, Y=Y, A=ALPHA, info=info)


if __name__ == '__main__':
    import sys
    res = clean(sys.argv[1])
    info = res['info']
    print('frames=%d  out_fps=%.2f  body=%.1fpx' % (info['n'], info['out_fps'], info['L_body']))
    for (a, b), L in info['canon'].items():
        print('  %-22s %.1f' % (a+'-'+b, L))
