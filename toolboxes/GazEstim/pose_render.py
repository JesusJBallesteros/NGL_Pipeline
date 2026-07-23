"""Render the cleaned skeleton + gaze cones over the real arena background.

Body-part labels are not hardcoded. Skeleton links, marker styles and the gaze
geometry are defined per ROLE (beak/head/back/left_wing/right_wing/tail) and mapped
onto whatever labels the data uses (info['roles'], produced by pose_clean). Selected
labels with no role are still drawn, with a generic colour and no bone links. Gaze
cones need at least a beak and a head role; other roles have fallbacks.

    render(cleaned, cfg)   # cleaned = output of pose_clean.clean(); cfg = dict
"""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.lines import Line2D
from matplotlib.patches import Polygon as MplPolygon, Patch
from matplotlib.colors import to_rgb
import matplotlib.animation as animation
from PIL import Image

CONN_ROLES = [('beak', 'head'), ('head', 'back'), ('head', 'left_wing'),
              ('head', 'right_wing'), ('back', 'left_wing'), ('back', 'right_wing'),
              ('back', 'tail')]
ROLE_STYLE = {'beak': ('#d24bd2', '*', 240), 'head': ('#e6e61f', 's', 95),
              'back': ('#e07820', 'o', 75), 'left_wing': ('#5bb6e6', 'o', 75),
              'right_wing': ('#5bb6e6', 'o', 75), 'tail': ('#3346c4', 'o', 75)}
EXTRA_COLORS = ['#8ad17a', '#d9a13b', '#c76ea0', '#7ec8c0', '#9b8ad1', '#b0b0b0']

GAZE_DEFAULTS = dict(enabled=True, mono_fov=170.0, bino_half=15.0,
                     cone_mult=2.5, eye_fwd_frac=1/3, eye_lat_frac=1/5)


def _build_conns(order, roles):
    return [(roles[a], roles[b]) for a, b in CONN_ROLES
            if a in roles and b in roles and roles[a] in order and roles[b] in order]


def _build_style(order, roles):
    lab2role = {v: k for k, v in roles.items()}
    style, extra = {}, 0
    for lab in order:
        r = lab2role.get(lab)
        if r in ROLE_STYLE:
            style[lab] = ROLE_STYLE[r]
        else:
            style[lab] = (EXTRA_COLORS[extra % len(EXTRA_COLORS)], 'o', 75)
            extra += 1
    return style


def _wedge(ax_x, ax_y, a0, a1, R, nseg=44):
    a = np.linspace(a0, a1, nseg)
    xs = np.concatenate(([ax_x], ax_x + R*np.cos(a)))
    ys = np.concatenate(([ax_y], ax_y + R*np.sin(a)))
    return np.column_stack([xs, ys])


def render(cleaned, cfg):
    X, Y, A, info = cleaned['X'], cleaned['Y'], cleaned['A'], cleaned['info']
    order = info['order']
    roles = info.get('roles', {})
    canon = info['canon']
    L_body = info['L_body']
    n = info['n']
    out_fps = cfg.get('out_fps', info['out_fps'])
    VW, VH = float(cfg['video_w']), float(cfg['video_h'])
    dpi = int(cfg.get('dpi', 120))
    if not (VW > 0 and VH > 0):
        raise ValueError('video_w/video_h must be positive, got %sx%s'
                         % (cfg['video_w'], cfg['video_h']))
    if not out_fps or out_fps <= 0:
        raise ValueError('output fps must be positive, got %r' % out_fps)
    nmax = min(n, cfg['max_frames']) if cfg.get('max_frames') else n

    vb = cfg.get('verbose', True)

    def log(m):
        if vb:
            print('[render] ' + m, flush=True)

    conns = _build_conns(order, roles)
    style = _build_style(order, roles)
    gz = {**GAZE_DEFAULTS, **(cfg.get('gaze') or {})}

    # ---- gaze geometry (needs beak + head; the rest has fallbacks) ----
    gaze_on = bool(gz['enabled']) and 'beak' in roles and 'head' in roles
    if gz['enabled'] and not gaze_on:
        log('gaze cones skipped: a beak and a head role are required')
    if gaze_on:
        BK, HD = roles['beak'], roles['head']

        def cl(a, b):
            return canon.get((a, b), canon.get((b, a)))

        OFF_FWD = (cl(BK, HD) or L_body) * gz['eye_fwd_frac']
        bw = []
        if 'back' in roles:
            bw = [v for v in (cl(roles['back'], roles[w])
                              for w in ('left_wing', 'right_wing') if w in roles) if v]
        OFF_LAT = (float(np.mean(bw)) if bw else L_body) * gz['eye_lat_frac']
        if 'tail' in roles:
            L_BIRD = float(np.median(np.hypot(X[BK]-X[roles['tail']], Y[BK]-Y[roles['tail']])))
        else:
            L_BIRD = 4.0 * L_body
            log('no tail role: bird length approximated as 4x body length (%.0f px)' % L_BIRD)
        CONE_R = gz['cone_mult'] * L_BIRD
        BINO_HALF = np.deg2rad(gz['bino_half'])
        MONO_FAR = np.deg2rad(gz['mono_fov']) - BINO_HALF
        hbx, hby = X[BK]-X[HD], Y[BK]-Y[HD]
        nrm = np.hypot(hbx, hby); nrm[nrm < 1e-6] = 1e-6
        fwx, fwy = hbx/nrm, hby/nrm
        px, py = -fwy, fwx
        TH = np.arctan2(fwy, fwx)
        EBX, EBY = X[HD] + fwx*OFF_FWD, Y[HD] + fwy*OFF_FWD
        E1X, E1Y = EBX + px*OFF_LAT, EBY + py*OFF_LAT
        E2X, E2Y = EBX - px*OFF_LAT, EBY - py*OFF_LAT
        log('gaze: eyes fwd=%.1f lat=%.1f px, cone R=%.0f px (mono %.0f deg, bino +/-%.0f deg)'
            % (OFF_FWD, OFF_LAT, CONE_R, gz['mono_fov'], gz['bino_half']))

    # ---- figure + background ----
    plate = np.asarray(Image.open(cfg['background']).convert('RGB'))
    log('background %s (%dx%d) -> stretched to %dx%d' %
        (os.path.basename(cfg['background']), plate.shape[1], plate.shape[0], VW, VH))
    # Canvas keeps the video's aspect, but H.264/yuv420p needs BOTH sides even, so
    # round down by at most 1 px. The +0.25 guards against matplotlib truncating
    # figsize*dpi downwards (e.g. 6.0833 in x 120 dpi -> 729.99 -> 729, an odd height).
    out_w = int(cfg.get('out_width') or 912)
    out_w -= out_w % 2
    out_h = int(round(out_w * VH / VW))
    out_h -= out_h % 2
    log('canvas %dx%d px (even dimensions required by H.264)' % (out_w, out_h))
    fig, ax = plt.subplots(figsize=((out_w + 0.25)/dpi, (out_h + 0.25)/dpi),
                           dpi=dpi, facecolor='black')
    ax.imshow(plate, extent=[0, VW, VH, 0], zorder=0, interpolation='bilinear')
    ax.set_xlim(0, VW); ax.set_ylim(VH, 0)
    ax.set_aspect('equal'); ax.axis('off')
    fig.subplots_adjust(left=0, right=1, top=1, bottom=0)

    cones = []
    if gaze_on:
        def _cone(color, fa, ea, z):
            r, g, b = to_rgb(color)
            poly = MplPolygon(np.zeros((3, 2)), closed=True, facecolor=(r, g, b, fa),
                              edgecolor=(r, g, b, ea), lw=0.9, zorder=z)
            ax.add_patch(poly); return poly
        cone_m1 = _cone('#ff2e2e', 0.12, 0.40, 3)
        cone_m2 = _cone('#ff2e2e', 0.12, 0.40, 3)
        cone_b = _cone('#2f7bff', 0.26, 0.60, 4)
        cones = [cone_m1, cone_m2, cone_b]

    bone_lines = {}
    for a_, b_ in conns:
        ln, = ax.plot([], [], color='#fafafa', lw=1.8, alpha=0.92, solid_capstyle='round',
                      path_effects=[pe.withStroke(linewidth=3.4, foreground='black')], zorder=5)
        bone_lines[(a_, b_)] = ln

    scat = {}
    for name in order:
        col, mk, sz = style[name]
        scat[name] = ax.scatter([], [], s=sz, marker=mk, facecolor=col,
                                edgecolor='#101010', linewidths=0.9, zorder=10)

    frame_txt = ax.text(0.012, 0.985, '', transform=ax.transAxes, ha='left', va='top',
                        color='white', fontsize=11, family='monospace', zorder=20,
                        bbox=dict(boxstyle='round,pad=0.3', facecolor='black',
                                  edgecolor='none', alpha=0.55))

    handles = [Line2D([0], [0], marker=style[k][1], linestyle='none',
                      markerfacecolor=style[k][0], markeredgecolor='white',
                      markersize=14 if style[k][1] == '*' else 11, label=k)
               for k in order]
    if gaze_on:
        handles += [Patch(facecolor=(1, 0.18, 0.18, 0.45), edgecolor='none',
                          label='monocular %d°' % round(gz['mono_fov'])),
                    Patch(facecolor=(0.18, 0.48, 1, 0.6), edgecolor='none',
                          label='binocular ±%d°' % round(gz['bino_half']))]
    leg = ax.legend(handles=handles, loc='upper right', ncol=2, fontsize=8.5,
                    facecolor='black', edgecolor='#777', labelcolor='white',
                    framealpha=1.0, handletextpad=0.45, columnspacing=1.0, borderpad=0.6)

    # Draw order for blitting, matching the z-order of a normal redraw:
    # cones(3,4) < bones(5) < legend(5, added last) < markers(10) < frame text(20)
    animated = cones + list(bone_lines.values()) + [leg] + list(scat.values()) + [frame_txt]

    def update(i):
        if gaze_on:
            t = TH[i]
            cone_m1.set_xy(_wedge(E1X[i], E1Y[i], t - BINO_HALF, t + MONO_FAR, CONE_R))
            cone_m2.set_xy(_wedge(E2X[i], E2Y[i], t - MONO_FAR, t + BINO_HALF, CONE_R))
            cone_b.set_xy(_wedge(EBX[i], EBY[i], t - BINO_HALF, t + BINO_HALF, CONE_R))
        for (a_, b_), ln in bone_lines.items():
            ln.set_data([X[a_][i], X[b_][i]], [Y[a_][i], Y[b_][i]])
        for name in order:
            scat[name].set_offsets([[X[name][i], Y[name][i]]])
            scat[name].set_alpha(float(A[name][i]))
        frame_txt.set_text('frame %4d / %d' % (i, nmax))
        return cones + list(bone_lines.values()) + list(scat.values()) + [frame_txt]

    # ---- preview still or full mp4 ----
    if cfg.get('preview_frame') is not None:
        update(int(cfg['preview_frame']))
        fig.savefig(cfg['output'], facecolor='black')
        plt.close(fig)
        log('saved preview frame %d -> %s' % (int(cfg['preview_frame']),
                                              os.path.basename(cfg['output'])))
        return cfg['output']

    log('rendering %d frames @ %.2f fps -> %s' % (nmax, out_fps, os.path.basename(cfg['output'])))
    if cfg.get('legacy_render'):
        # slow reference path: full figure redraw per frame via matplotlib's writer
        ani = animation.FuncAnimation(fig, update, frames=nmax, blit=False, interval=1000/out_fps)
        w = animation.FFMpegWriter(fps=out_fps, codec='libx264',
                                   extra_args=['-pix_fmt', 'yuv420p',
                                               '-crf', str(cfg.get('crf', 24)),
                                               '-preset', cfg.get('preset', 'veryfast')],
                                   metadata={'title': 'Jackdaw gaze'})
        step_p = max(1, nmax // 20)

        def _progress(i, nn):
            if vb and (i % step_p == 0 or i == nmax - 1):
                print('[render]   frame %d/%d (%d%%)' % (i+1, nmax, round(100.0*(i+1)/nmax)), flush=True)
        ani.save(cfg['output'], writer=w, dpi=dpi, progress_callback=_progress)
    else:
        _write_video_blit(fig, ax, update, animated, nmax, out_fps, cfg, vb)
    plt.close(fig)
    log('saved %s' % cfg['output'])
    return cfg['output']


def _write_video_blit(fig, ax, update, animated, nmax, out_fps, cfg, vb):
    """Fast path: the background (arena, axes) is rasterised once and cached; each
    frame only redraws the handful of moving artists, and the raw RGBA buffer is
    piped straight to ffmpeg. Pixel-identical to a full redraw, ~9x faster."""
    import subprocess
    import threading

    canvas = fig.canvas
    for a in animated:
        a.set_animated(True)          # excluded from the cached background
    canvas.draw()
    bg = canvas.copy_from_bbox(fig.bbox)
    H, W = np.asarray(canvas.buffer_rgba()).shape[:2]
    # Last-resort guard: whatever matplotlib actually rasterised, hand ffmpeg even
    # dimensions (libx264 + yuv420p rejects odd sides). Costs at most 1 px.
    W2, H2 = W - (W % 2), H - (H % 2)
    crop = (W2, H2) != (W, H)
    if crop and vb:
        print('[render] trimming canvas %dx%d -> %dx%d (H.264 needs even sides)'
              % (W, H, W2, H2), flush=True)

    exe = matplotlib.rcParams.get('animation.ffmpeg_path') or 'ffmpeg'
    cmd = [exe, '-v', 'error', '-y',
           '-f', 'rawvideo', '-pix_fmt', 'rgba', '-s', '%dx%d' % (W2, H2),
           '-r', '%.6f' % out_fps, '-i', '-', '-an',
           '-c:v', cfg.get('codec', 'libx264'),
           '-crf', str(cfg.get('crf', 24)),
           '-preset', cfg.get('preset', 'veryfast'),
           '-pix_fmt', 'yuv420p', cfg['output']]
    proc = subprocess.Popen(cmd, stdin=subprocess.PIPE, stderr=subprocess.PIPE)

    # Drain stderr on a thread: a full stderr pipe would block ffmpeg, and we need
    # its message to explain a failure (a dead ffmpeg shows up here only as EINVAL).
    err_lines = []

    def _drain():
        for line in proc.stderr:
            err_lines.append(line.decode('utf8', 'replace'))
    th = threading.Thread(target=_drain, daemon=True)
    th.start()

    def _died(i, exc):
        proc.wait()
        th.join(timeout=2)
        said = ''.join(err_lines).strip() or ('%s: %s' % (type(exc).__name__, exc))
        return RuntimeError(
            'ffmpeg stopped consuming frames after %d frame(s) (exit code %s).\n'
            '  ffmpeg said: %s\n'
            '  output was : %s\n'
            '  usual causes: the output file is open in a video player, the folder is '
            'not writable, or the disk is full.' % (i, proc.returncode, said, cfg['output']))

    step_p = max(1, nmax // 20)
    try:
        for i in range(nmax):
            canvas.restore_region(bg)
            update(i)
            for a in animated:
                ax.draw_artist(a)
            # .tobytes() -> flat, C-contiguous; safest thing to hand to a pipe
            arr = np.asarray(canvas.buffer_rgba())
            if crop:
                arr = arr[:H2, :W2]
            try:
                proc.stdin.write(arr.tobytes())
            except (BrokenPipeError, OSError) as exc:
                raise _died(i, exc) from None
            if vb and (i % step_p == 0 or i == nmax - 1):
                print('[render]   frame %d/%d (%d%%)' % (i+1, nmax, round(100.0*(i+1)/nmax)), flush=True)
    finally:
        try:
            if proc.stdin and not proc.stdin.closed:
                proc.stdin.close()
        except (BrokenPipeError, OSError):
            pass
        rc = proc.wait()
        th.join(timeout=2)
    if rc != 0:
        raise RuntimeError('ffmpeg failed (exit %d): %s'
                           % (rc, ''.join(err_lines).strip()[:600]))
