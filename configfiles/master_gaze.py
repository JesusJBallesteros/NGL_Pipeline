"""Master pipeline: clean a DeepLabCut pose table and render the gaze animation.

Usage:   python master_gaze.py <config.json>

The config is normally written by the MATLAB wrapper (functions/video/process_gaze.m)
but can be hand-written. It locates the python tools in ../toolboxes (or an explicit
"toolboxes" path in the config), runs pose_clean -> pose_render, and prints a single
RESULT json line with the output path and a quality summary.
"""
import os
import sys
import json


def _emit(result, cfg):
    """Print the machine-readable result and, when asked, also drop it in a file.
    The MATLAB wrapper reads that file because it deliberately does NOT capture
    stdout any more - that is what lets the progress messages stream live."""
    print('RESULT ' + json.dumps(result), flush=True)
    path = (cfg or {}).get('result_json')
    if path:
        try:
            with open(path, 'w') as fh:
                json.dump(result, fh)
        except OSError:
            pass
    return result


def main(config_path):
    with open(config_path, 'r') as f:
        cfg = json.load(f)
    verbose = cfg.get('verbose', True)

    def log(msg):
        if verbose:
            print('[master] ' + msg, flush=True)
    log('config: %s' % os.path.basename(config_path))

    here = os.path.dirname(os.path.abspath(__file__))
    root = os.path.dirname(here)
    toolboxes = cfg.get('toolboxes') or os.path.join(root, 'toolboxes')
    sys.path.insert(0, toolboxes)
    import pose_clean
    import pose_render

    if not cfg.get('csv') or not os.path.isfile(cfg['csv']):
        raise FileNotFoundError('csv not found: %r' % cfg.get('csv'))

    if cfg.get('list_parts'):          # just report the labels available in the table
        parts = pose_clean.available_parts(cfg['csv'])
        log('parts in %s: %s' % (os.path.basename(cfg['csv']), ', '.join(parts)))
        return _emit(dict(status='ok', parts=parts), cfg)

    # a preview still is a render too, so honour it even when video is false
    want_video = bool(cfg.get('video', True)) or cfg.get('preview_frame') is not None
    if want_video and not os.path.isfile(cfg.get('background', '')):
        raise FileNotFoundError('background not found: %r' % cfg.get('background'))

    # clean
    clean_cfg = dict(cfg.get('clean') or {})
    for key in ('fps_in', 'downsample_step', 'target_fps', 'start_time', 'end_time',
                'parts', 'roles', 'body_px'):
        if cfg.get(key) is not None:
            clean_cfg[key] = cfg[key]
    clean_cfg['verbose'] = verbose
    log('stage 1/2 - cleaning pose')
    cleaned = pose_clean.clean(cfg['csv'], clean_cfg)
    info = cleaned['info']

    # optional extra output: derived features (.mat) + summary figure
    feat_path, fig_path = '', ''
    want_feat = cfg.get('features', cfg.get('head_direction'))     # old name still works
    if want_feat:
        import features
        target = want_feat if isinstance(want_feat, str) else features.default_path(cfg.get('output'))
        data = features.compute(cleaned)
        feat_path = os.path.abspath(features.save_mat(data, target))
        log('features -> %s' % feat_path)
        want_fig = cfg.get('features_figure', True)
        if want_fig:
            fig_target = want_fig if isinstance(want_fig, str) else os.path.splitext(target)[0] + '.png'
            fig_path = os.path.abspath(features.plot(
                data, fig_target, video_w=cfg.get('video_w'), video_h=cfg.get('video_h'),
                min_likelihood=(cfg.get('clean') or {}).get('p_cut', 0.5)))
            log('figure -> %s' % fig_path)

    # render
    out = ''
    if want_video:
        render_cfg = dict(
            background=cfg['background'], output=cfg['output'],
            video_w=cfg.get('video_w', 1250), video_h=cfg.get('video_h', 1160),
            out_width=cfg.get('out_width'),
            dpi=cfg.get('dpi', 120), crf=cfg.get('crf', 24), preset=cfg.get('preset', 'veryfast'),
            max_frames=cfg.get('max_frames'), preview_frame=cfg.get('preview_frame'),
            gaze=cfg.get('gaze'), verbose=verbose,
        )
        os.makedirs(os.path.dirname(os.path.abspath(cfg['output'])), exist_ok=True)
        log('stage 2/2 - rendering')
        out = os.path.abspath(pose_render.render(cleaned, render_cfg))
        log('done -> %s' % out)
    else:
        log('video rendering skipped (video = false)')

    flags = {k: int(info['FLAG'][k].sum()) for k in info['order']}
    result = dict(status='ok', output=out, features=feat_path, figure=fig_path,
                  frames=info['n'],
                  out_fps=round(info['out_fps'], 3), body_px=round(info['L_body'], 1),
                  downsample_step=info['step'], flagged=flags,
                  parts=info['order'], roles=info['roles'], available=info['available'],
                  canon={('%s-%s' % k): round(v, 1) for k, v in info['canon'].items()})
    return _emit(result, cfg)


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('RESULT ' + json.dumps({'status': 'error', 'message': 'no config path given'}))
        sys.exit(2)
    try:
        main(sys.argv[1])
    except Exception as exc:
        import traceback
        traceback.print_exc()
        cfg = {}
        try:
            with open(sys.argv[1], 'r') as fh:
                cfg = json.load(fh)
        except Exception:
            pass
        _emit({'status': 'error', 'message': str(exc)}, cfg)   # so MATLAB still sees why
        sys.exit(1)
