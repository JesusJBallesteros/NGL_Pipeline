"""Derived kinematic features from the cleaned pose, saved as a .mat + summary figure.

Computed per frame
------------------
  angle_deg        head direction from the beak-head vector.
                   0 deg = video vertical (up), increasing CLOCKWISE, directional
                   head -> beak. Image y grows downwards, so screen-up is -y.
  angvel_deg_s     angular velocity of that vector. The angle is CIRCULAR, so it is
                   unwrapped before differentiating, and the derivative uses the real
                   (non-uniform) timestamps via central differences.
  centroid_x/y     centre of gravity = unweighted mean of all selected labels.
  speed_px_s       speed of that centroid.

Confidence
----------
Three complementary things, because they mean different things:
  *_likelihood  - DLC likelihood aggregated over the labels a quantity depends on
                  (min for pairs: the weakest link governs; mean for the centroid,
                  where every label contributes equally).
  *_estimated   - the value depends on a joint the cleaner had to re-estimate, i.e.
                  it is interpolated rather than observed.
  *_sd / *_ci95 - APPROXIMATE uncertainty propagated from the per-label positional
                  scatter (`pos_sigma`, a robust px noise estimate):
                    angle : atan2(sigma_pos, lever)      -> short lever = noisier
                    rates : sigma of the two samples / their time separation
                  This is geometric error propagation, NOT a calibrated statistical
                  CI, and it understates error on `estimated` samples.
"""
import os
import numpy as np

ANGVEL_LIMIT_DEG_S = 500.0      # y-axis range of the angular-velocity vs direction panel


# ----------------------------------------------------------------- helpers
def _col(a):
    return np.asarray(a, dtype=float).reshape(-1, 1)


def _neighbour_span(t):
    """t[i+1]-t[i-1] for central differences, with one-sided ends."""
    n = len(t)
    d = np.empty(n)
    if n < 2:
        return np.full(n, np.nan)
    d[1:-1] = t[2:] - t[:-2]
    d[0] = t[1] - t[0]
    d[-1] = t[-1] - t[-2]
    return d


def _neighbour_combine(a, fn):
    """fn applied over (i-1, i, i+1), edges clamped."""
    n = len(a)
    if n == 1:
        return np.asarray(a, float)
    prev = np.r_[a[0], a[:-1]]
    nxt = np.r_[a[1:], a[-1]]
    return fn(np.vstack([prev, a, nxt]), axis=0)


def _pair_sigma(sd):
    """sqrt(sd[i-1]^2 + sd[i+1]^2) for a central difference."""
    n = len(sd)
    if n < 2:
        return np.full(n, np.nan)
    prev = np.r_[sd[0], sd[:-1]]
    nxt = np.r_[sd[1:], sd[-1]]
    return np.hypot(prev, nxt)


# ----------------------------------------------------------------- compute
def compute(cleaned):
    X, Y, info = cleaned['X'], cleaned['Y'], cleaned['info']
    labels = list(info['order'])
    roles = info.get('roles', {})
    P = info['P']
    FLAG = info['FLAG']
    sig = info.get('pos_sigma', {}) or {}
    n = int(info['n'])

    t_ms = np.asarray(info['frame_index'], float) / float(info['fps_in']) * 1000.0
    t_s = t_ms / 1000.0
    span = _neighbour_span(t_s)

    out = {'timestamp_ms': _col(t_ms),
           'frame_index': _col(np.asarray(info['frame_index'], float))}

    # ---------------- head direction + angular velocity ----------------
    has_angle = ('beak' in roles) and ('head' in roles)
    if has_angle:
        bk, hd = roles['beak'], roles['head']
        dx = np.asarray(X[bk], float) - np.asarray(X[hd], float)
        dy = np.asarray(Y[bk], float) - np.asarray(Y[hd], float)
        lever = np.hypot(dx, dy)
        angle = np.degrees(np.arctan2(dx, -dy)) % 360.0

        p_bk, p_hd = np.asarray(P[bk], float), np.asarray(P[hd], float)
        a_lik = np.minimum(p_bk, p_hd)
        a_est = np.asarray(FLAG[bk], bool) | np.asarray(FLAG[hd], bool)
        s_pos = np.hypot(np.nan_to_num(float(sig.get(bk, np.nan))),
                         np.nan_to_num(float(sig.get(hd, np.nan))))
        a_sd = np.degrees(np.arctan2(s_pos, np.maximum(lever, 1e-6)))

        # circular -> unwrap before differentiating, and use the real timestamps
        if n >= 2:
            unw = np.unwrap(np.radians(angle))
            angvel = np.degrees(np.gradient(unw, t_s))
        else:
            angvel = np.full(n, np.nan)
        w_sd = _pair_sigma(a_sd) / np.maximum(span, 1e-9)
        w_lik = _neighbour_combine(a_lik, np.min)
        w_est = _neighbour_combine(a_est.astype(float), np.max) > 0

        out.update({
            'angle_deg': _col(angle), 'angle_sd_deg': _col(a_sd),
            'angle_ci95_deg': _col(1.96*a_sd), 'angle_likelihood': _col(a_lik),
            'angle_estimated': _col(a_est.astype(np.uint8)), 'lever_px': _col(lever),
            'beak_likelihood': _col(p_bk), 'head_likelihood': _col(p_hd),
            'angvel_deg_s': _col(angvel), 'angvel_sd_deg_s': _col(w_sd),
            'angvel_ci95_deg_s': _col(1.96*w_sd), 'angvel_likelihood': _col(w_lik),
            'angvel_estimated': _col(w_est.astype(np.uint8)),
        })

    # ---------------- centre of gravity + speed ----------------
    Xs = np.vstack([np.asarray(X[k], float) for k in labels])
    Ys = np.vstack([np.asarray(Y[k], float) for k in labels])
    Ps = np.vstack([np.asarray(P[k], float) for k in labels])
    Es = np.vstack([np.asarray(FLAG[k], bool) for k in labels])
    with np.errstate(invalid='ignore'):
        cx = np.nanmean(Xs, axis=0)
        cy = np.nanmean(Ys, axis=0)
    c_lik = np.nanmean(Ps, axis=0)                 # every label contributes equally
    c_est_frac = Es.mean(axis=0)

    sig_vals = np.array([float(sig.get(k, np.nan)) for k in labels], float)
    sig_vals = sig_vals[np.isfinite(sig_vals)]
    c_sigma = float(np.sqrt(np.sum(sig_vals**2)) / max(len(labels), 1)) if sig_vals.size else np.nan

    if n >= 2:
        vx = np.gradient(cx, t_s)
        vy = np.gradient(cy, t_s)
        speed = np.hypot(vx, vy)
    else:
        vx = vy = speed = np.full(n, np.nan)
    v_sd = (np.sqrt(2.0) * c_sigma) / np.maximum(span, 1e-9)
    v_lik = _neighbour_combine(c_lik, np.min)

    out.update({
        'centroid_x': _col(cx), 'centroid_y': _col(cy),
        'centroid_likelihood': _col(c_lik),
        'centroid_estimated_frac': _col(c_est_frac),
        'speed_px_s': _col(speed), 'speed_sd_px_s': _col(v_sd),
        'speed_ci95_px_s': _col(1.96*v_sd), 'speed_likelihood': _col(v_lik),
        'velocity_x_px_s': _col(vx), 'velocity_y_px_s': _col(vy),
    })

    # ---------------- metadata ----------------
    out.update({
        'parts': np.array(labels, dtype=object),
        'has_angle': int(has_angle),
        'beak_label': roles.get('beak', ''), 'head_label': roles.get('head', ''),
        'fps_in': float(info['fps_in']), 'fps_out': float(info['out_fps']),
        'centroid_sigma_px': c_sigma, 'n': n,
        'convention': ('angle: 0 deg = video vertical (up), clockwise, directional '
                       'head->beak; angular velocity from the unwrapped angle using '
                       'the real (non-uniform) timestamps'),
        'note': ('*_sd is geometric error propagation, not a calibrated CI, and '
                 'understates error where *_estimated is 1'),
    })
    return out


def save_mat(data, path):
    try:
        from scipy.io import savemat
    except ImportError:
        raise RuntimeError('saving the features .mat needs scipy: '
                           'run "python -m pip install scipy"') from None
    d = os.path.dirname(os.path.abspath(path))
    if d:
        os.makedirs(d, exist_ok=True)
    savemat(path, {'features': data}, do_compression=True)
    return path


def default_path(video_output):
    base = os.path.splitext(video_output or 'output')[0]
    return base + '_estimated_features.mat'


# ----------------------------------------------------------------- figure
def plot(data, path, video_w=None, video_h=None, min_likelihood=0.5):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    g = lambda k: np.asarray(data[k]).ravel() if k in data else None
    t = g('timestamp_ms') / 1000.0
    speed, cx, cy = g('speed_px_s'), g('centroid_x'), g('centroid_y')
    has_angle = bool(data.get('has_angle', 0))

    fig, axes = plt.subplots(2, 2, figsize=(12, 9))
    fig.suptitle('Estimated features', fontsize=14, y=0.98)
    ax_rose = fig.add_subplot(2, 2, 1, projection='polar')
    axes[0, 0].remove()
    ax_wa, ax_hm, ax_v = axes[0, 1], axes[1, 0], axes[1, 1]

    # ---- 1. distribution of head angles (polar rose) ----
    if has_angle:
        ang, alik, aest = g('angle_deg'), g('angle_likelihood'), g('angle_estimated')
        keep = (alik >= min_likelihood) & (aest == 0) & np.isfinite(ang)
        a = ang[keep]
        ax_rose.set_theta_zero_location('N')      # 0 deg at top ...
        ax_rose.set_theta_direction(-1)           # ... increasing clockwise
        if a.size:
            edges = np.deg2rad(np.arange(0, 361, 10))
            cnt, _ = np.histogram(np.deg2rad(a), bins=edges)
            ax_rose.bar(edges[:-1], cnt, width=np.diff(edges), align='edge',
                        color='#4a7fc1', edgecolor='white', linewidth=0.4)
            m = np.degrees(np.arctan2(np.mean(np.sin(np.deg2rad(a))),
                                      np.mean(np.cos(np.deg2rad(a))))) % 360
            ax_rose.plot([np.deg2rad(m)]*2, [0, cnt.max()], color='#c0392b', lw=2)
            ax_rose.set_title('Head direction (n=%d, %.0f%% used)\nmean %.0f deg'
                              % (a.size, 100*keep.mean(), m), fontsize=10, pad=18)
        ax_rose.set_yticklabels([])
    else:
        ax_rose.text(0, 0, 'no beak/head roles', ha='center')

    # ---- 2. angular velocity vs angle ----
    if has_angle:
        w = g('angvel_deg_s')
        keep2 = keep & np.isfinite(w)
        lim = ANGVEL_LIMIT_DEG_S
        # drop out-of-range samples rather than clipping them: clipping would pile
        # every outlier into the edge rows and fake a bright band at +/-lim
        inr = keep2 & (np.abs(w) <= lim)
        if inr.sum() > 10:
            ax_wa.hist2d(ang[inr], w[inr], bins=[36, 40],
                         range=[[0, 360], [-lim, lim]], cmap='magma')
            ax_wa.axhline(0, color='w', lw=0.8, alpha=0.6)
        beyond = 100.0 * (keep2 & (np.abs(w) > lim)).sum() / max(keep2.sum(), 1)
        ax_wa.set_xlabel('head direction (deg, 0=up, CW)')
        ax_wa.set_ylabel('angular velocity (deg/s)')
        ax_wa.set_title('Angular velocity vs direction\n(+/-%d deg/s shown, %.1f%% off-scale)'
                        % (lim, beyond), fontsize=10)
        ax_wa.set_xlim(0, 360); ax_wa.set_xticks(range(0, 361, 90))
        ax_wa.set_ylim(-lim, lim)
    else:
        ax_wa.axis('off')

    # ---- 3. location heatmap (centre of gravity) ----
    ok = np.isfinite(cx) & np.isfinite(cy)
    rng = [[0, video_w], [0, video_h]] if video_w and video_h else None
    if ok.sum() > 10:
        from matplotlib.colors import LogNorm
        # log scale + empty bins masked: perches otherwise saturate the map and hide
        # everywhere the animal merely passed through
        h = ax_hm.hist2d(cx[ok], cy[ok], bins=60, range=rng, cmap='viridis',
                         cmin=1, norm=LogNorm())
        fig.colorbar(h[3], ax=ax_hm, label='frames (log)', fraction=0.046)
    ax_hm.set_aspect('equal')
    if video_w and video_h:
        ax_hm.set_xlim(0, video_w); ax_hm.set_ylim(video_h, 0)   # video orientation
    else:
        ax_hm.invert_yaxis()
    ax_hm.set_xlabel('x (px)'); ax_hm.set_ylabel('y (px)')
    ax_hm.set_title('Occupancy of the centre of gravity', fontsize=10)

    # ---- 4. speed over time ----
    dec = max(1, len(t) // 20000)          # keep the line drawable on long recordings
    ax_v.plot(t[::dec], speed[::dec], lw=0.6, color='#2c6fbb')
    med = np.nanmedian(speed)
    ax_v.axhline(med, color='#c0392b', lw=1.2, ls='--', label='median %.0f px/s' % med)
    ax_v.set_xlabel('time (s from recording start)')
    ax_v.set_ylabel('speed (px/s)')
    ax_v.set_title('Subject speed (centre of gravity)', fontsize=10)
    ax_v.legend(fontsize=8)
    top = np.nanpercentile(speed, 99.5) if np.isfinite(speed).any() else 1
    ax_v.set_ylim(0, max(top, 1e-6))

    fig.tight_layout(rect=[0, 0, 1, 0.96])
    d = os.path.dirname(os.path.abspath(path))
    if d:
        os.makedirs(d, exist_ok=True)
    fig.savefig(path, dpi=130)
    plt.close(fig)
    return path
