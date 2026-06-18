"""Render the cleaned jackdaw skeleton + gaze cones over the real arena background.

Parameterized, importable version of the development `animate.py`. Visual design is
preserved: skeleton with dark-halo bones, confidence-faded markers, two red monocular
cones + one blue binocular cone, photo background stretched to the true video pixel
size, legend in a black box top-right, frame counter top-left.

    render(cleaned, cfg)   # cleaned = output of pose_clean.clean(); cfg = dict (see below)
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

CONNS = [('beak', 'head'), ('head', 'back'), ('head', 'left_wing'),
         ('head', 'right_wing'), ('back', 'left_wing'), ('back', 'right_wing'),
         ('back', 'tail')]
STYLE = {'beak': ('#d24bd2', '*', 240), 'head': ('#e6e61f', 's', 95),
         'back': ('#e07820', 'o', 75), 'left_wing': ('#5bb6e6', 'o', 75),
         'right_wing': ('#5bb6e6', 'o', 75), 'tail': ('#3346c4', 'o', 75)}

GAZE_DEFAULTS = dict(enabled=True, mono_fov=170.0, bino_half=15.0,
                     cone_mult=2.5, eye_fwd_frac=1/3, eye_lat_frac=1/5)


def _wedge(ax_x, ax_y, a0, a1, R, nseg=44):
    a = np.linspace(a0, a1, nseg)
    xs = np.concatenate(([ax_x], ax_x + R*np.cos(a)))
    ys = np.concatenate(([ax_y], ax_y + R*np.sin(a)))
    return np.column_stack([xs, ys])


def render(cleaned, cfg):
    X, Y, A, info = cleaned['X'], cleaned['Y'], cleaned['A'], cleaned['info']
    order = info['order']
    n = info['n']
    out_fps = cfg.get('out_fps', info['out_fps'])
    VW, VH = cfg['video_w'], cfg['video_h']
    dpi = cfg.get('dpi', 120)
    nmax = min(n, cfg['max_frames']) if cfg.get('max_frames') else n

    vb = cfg.get('verbose', True)

    def log(m):
        if vb:
            print('[render] ' + m, flush=True)

    gz = {**GAZE_DEFAULTS, **(cfg.get('gaze') or {})}
    canon = info['canon']

    # ---- gaze geometry ----
    if gz['enabled']:
        OFF_FWD = canon[('beak', 'head')] * gz['eye_fwd_frac']
        OFF_LAT = 0.5*(canon[('back', 'left_wing')] + canon[('back', 'right_wing')]) * gz['eye_lat_frac']
        L_BIRD = float(np.median(np.hypot(X['beak']-X['tail'], Y['beak']-Y['tail'])))
        CONE_R = gz['cone_mult'] * L_BIRD
        BINO_HALF = np.deg2rad(gz['bino_half'])
        MONO_FAR = np.deg2rad(gz['mono_fov']) - BINO_HALF
        hbx, hby = X['beak']-X['head'], Y['beak']-Y['head']
        nrm = np.hypot(hbx, hby); nrm[nrm < 1e-6] = 1e-6
        fwx, fwy = hbx/nrm, hby/nrm
        px, py = -fwy, fwx
        TH = np.arctan2(fwy, fwx)
        EBX, EBY = X['head'] + fwx*OFF_FWD, Y['head'] + fwy*OFF_FWD
        E1X, E1Y = EBX + px*OFF_LAT, EBY + py*OFF_LAT
        E2X, E2Y = EBX - px*OFF_LAT, EBY - py*OFF_LAT
        log('gaze: eyes fwd=%.1f lat=%.1f px, cone R=%.0f px (mono %.0f deg, bino +/-%.0f deg)'
            % (OFF_FWD, OFF_LAT, CONE_R, gz['mono_fov'], gz['bino_half']))

    # ---- figure + background ----
    plate = np.asarray(Image.open(cfg['background']).convert('RGB'))
    log('background %s (%dx%d) -> stretched to %dx%d' %
        (os.path.basename(cfg['background']), plate.shape[1], plate.shape[0], VW, VH))
    fig, ax = plt.subplots(figsize=(7.6, 7.6 * VH / VW), dpi=dpi, facecolor='black')
    ax.imshow(plate, extent=[0, VW, VH, 0], zorder=0, interpolation='bilinear')
    ax.set_xlim(0, VW); ax.set_ylim(VH, 0)
    ax.set_aspect('equal'); ax.axis('off')
    fig.subplots_adjust(left=0, right=1, top=1, bottom=0)

    cones = []
    if gz['enabled']:
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
    for a_, b_ in CONNS:
        ln, = ax.plot([], [], color='#fafafa', lw=1.8, alpha=0.92, solid_capstyle='round',
                      path_effects=[pe.withStroke(linewidth=3.4, foreground='black')], zorder=5)
        bone_lines[(a_, b_)] = ln

    scat = {}
    for name in order:
        col, mk, sz = STYLE[name]
        scat[name] = ax.scatter([], [], s=sz, marker=mk, facecolor=col,
                                edgecolor='#101010', linewidths=0.9, zorder=10)

    frame_txt = ax.text(0.012, 0.985, '', transform=ax.transAxes, ha='left', va='top',
                        color='white', fontsize=11, family='monospace', zorder=20,
                        bbox=dict(boxstyle='round,pad=0.3', facecolor='black', edgecolor='none', alpha=0.55))

    handles = [Line2D([0], [0], marker=STYLE[k][1], linestyle='none',
                      markerfacecolor=STYLE[k][0], markeredgecolor='white',
                      markersize=11 if k != 'beak' else 14,
                      label={'left_wing': 'wing L', 'right_wing': 'wing R'}.get(k, k))
               for k in ['beak', 'head', 'back', 'left_wing', 'right_wing', 'tail']]
    if gz['enabled']:
        handles += [Patch(facecolor=(1, 0.18, 0.18, 0.45), edgecolor='none',
                          label='monocular %d°' % round(gz['mono_fov'])),
                    Patch(facecolor=(0.18, 0.48, 1, 0.6), edgecolor='none',
                          label='binocular ±%d°' % round(gz['bino_half']))]
    ax.legend(handles=handles, loc='upper right', ncol=2, fontsize=8.5,
              facecolor='black', edgecolor='#777', labelcolor='white',
              framealpha=1.0, handletextpad=0.45, columnspacing=1.0, borderpad=0.6)

    def update(i):
        if gz['enabled']:
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
        log('saved preview frame %d -> %s' % (int(cfg['preview_frame']), os.path.basename(cfg['output'])))
        return cfg['output']

    ani = animation.FuncAnimation(fig, update, frames=nmax, blit=False, interval=1000/out_fps)
    w = animation.FFMpegWriter(fps=out_fps, codec='libx264',
                               extra_args=['-pix_fmt', 'yuv420p',
                                           '-crf', str(cfg.get('crf', 24)),
                                           '-preset', cfg.get('preset', 'veryfast')],
                               metadata={'title': 'Jackdaw gaze'})
    log('rendering %d frames @ %.2f fps -> %s' % (nmax, out_fps, os.path.basename(cfg['output'])))
    step_p = max(1, nmax // 20)

    def _progress(i, nn):
        if vb and (i % step_p == 0 or i == nmax - 1):
            print('[render]   frame %d/%d (%d%%)' % (i + 1, nmax, round(100.0*(i+1)/nmax)), flush=True)
    ani.save(cfg['output'], writer=w, dpi=dpi, progress_callback=_progress)
    plt.close(fig)
    log('saved %s' % cfg['output'])
    return cfg['output']
