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
    if not os.path.isfile(cfg.get('background', '')):
        raise FileNotFoundError('background not found: %r' % cfg.get('background'))

    # clean
    clean_cfg = dict(cfg.get('clean') or {})
    for key in ('fps_in', 'downsample_step', 'target_fps', 'start_time', 'end_time'):
        if cfg.get(key) is not None:
            clean_cfg[key] = cfg[key]
    clean_cfg['verbose'] = verbose
    log('stage 1/2 - cleaning pose')
    cleaned = pose_clean.clean(cfg['csv'], clean_cfg)
    info = cleaned['info']

    # render
    render_cfg = dict(
        background=cfg['background'], output=cfg['output'],
        video_w=cfg.get('video_w', 1250), video_h=cfg.get('video_h', 1160),
        dpi=cfg.get('dpi', 120), crf=cfg.get('crf', 24), preset=cfg.get('preset', 'veryfast'),
        max_frames=cfg.get('max_frames'), preview_frame=cfg.get('preview_frame'),
        gaze=cfg.get('gaze'), verbose=verbose,
    )
    os.makedirs(os.path.dirname(os.path.abspath(cfg['output'])), exist_ok=True)
    log('stage 2/2 - rendering')
    out = pose_render.render(cleaned, render_cfg)
    log('done -> %s' % out)

    flags = {k: int(info['FLAG'][k].sum()) for k in info['order']}
    result = dict(status='ok', output=os.path.abspath(out), frames=info['n'],
                  out_fps=round(info['out_fps'], 3), body_px=round(info['L_body'], 1),
                  downsample_step=info['step'], flagged=flags,
                  canon={('%s-%s' % k): round(v, 1) for k, v in info['canon'].items()})
    print('RESULT ' + json.dumps(result))
    return result


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('RESULT ' + json.dumps({'status': 'error', 'message': 'no config path given'}))
        sys.exit(2)
    try:
        main(sys.argv[1])
    except Exception as exc:
        import traceback
        traceback.print_exc()
        print('RESULT ' + json.dumps({'status': 'error', 'message': str(exc)}))
        sys.exit(1)
