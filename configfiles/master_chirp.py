"""Master pipeline: survey a session's raw INTAN channels with CHIRP.

Usage:   python master_chirp.py <config.json>

The config is normally written by the MATLAB wrapper (functions/intan/process_chirp.m)
but can be hand-written. It locates the CHIRP engine in ../toolboxes/CHIRP (or an
explicit "toolboxes" path in the config), runs one statistics survey per channel
group, writes a combined CSV and a text report, and prints a single RESULT json
line with the output paths and a short summary.

Channels are surveyed in groups because the cross-channel sharing measure is
relative to the size of the group it is computed over: CHIRP calls a waveform
"noise" when it appears on more than max(8, 0.20 x N) of the N channels it saw.
With a multi-area probe the pipeline passes one group per area, so an area is
judged against its own channel count rather than against the whole headstage.
"""
import os
import sys
import json
from pathlib import Path


def _emit(result, cfg):
    """Print the machine-readable result and, when asked, also drop it in a file.
    The MATLAB wrapper reads that file because it deliberately does NOT capture
    stdout - that is what lets the scan progress stream live."""
    print('RESULT ' + json.dumps(result), flush=True)
    path = (cfg or {}).get('result_json')
    if path:
        try:
            with open(path, 'w') as fh:
                json.dump(result, fh)
        except OSError:
            pass
    return result


def _write_csv(rows, path, fields):
    import csv
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, 'w', newline='', encoding='utf-8') as fh:
        wr = csv.DictWriter(fh, fieldnames=fields, extrasaction='ignore')
        wr.writeheader()
        for r in rows:
            wr.writerow({k: (round(v, 4) if isinstance(v, float) else v)
                         for k, v in r.items()})
    return os.path.abspath(path)


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
    toolboxes = cfg.get('toolboxes') or os.path.join(root, 'toolboxes', 'CHIRP')
    sys.path.insert(0, toolboxes)
    import dat_to_stats as eng_stats
    import dat_to_video as eng_video

    # The wrapper depends on run_stats(), which only exists from CHIRP 1.3.0.
    # Fail here with a clear message rather than on a missing attribute.
    version = getattr(eng_stats, '__version__', '0')
    if tuple(int(p) for p in version.split('.')[:2]) < (1, 3):
        raise RuntimeError('CHIRP %s in %s is too old; needs 1.3.0 or newer '
                           '(the vendored copy under toolboxes/CHIRP is stale)'
                           % (version, toolboxes))
    log('CHIRP %s from %s' % (version, toolboxes))

    groups = cfg.get('groups') or []
    if not groups:
        raise ValueError('no channel groups given')
    out_dir = cfg.get('out_dir') or os.getcwd()
    os.makedirs(out_dir, exist_ok=True)

    band = tuple(cfg.get('band', (450.0, 8000.0)))
    params = eng_stats.StatsParams(
        duration=cfg.get('duration', 10.0),
        start=cfg.get('start'),
        band=band,
        fs=int(cfg.get('fs', 30000)),
        neg_k=cfg.get('neg_k', 5.0),
        pos_k=cfg.get('pos_k', 8.0),
        artifact_k=cfg.get('artifact_k', 18.0),
        max_k=int(cfg.get('max_k', 3)),
        segments=int(cfg.get('segments', eng_stats.STAT_SEGMENTS)),
        scan_step=cfg.get('scan_step', 60.0),
        scan_range=tuple(cfg['scan_range']) if cfg.get('scan_range') else None,
    ).validate()

    all_rows, sections, results = [], [], {}
    for gi, group in enumerate(groups, start=1):
        name = group.get('name') or ''
        files = [f for f in group.get('files', [])]
        missing = [f for f in files if not os.path.isfile(f)]
        if missing:
            raise FileNotFoundError('missing .dat file(s): %s'
                                    % ', '.join(missing[:5]))
        if not files:
            log('group %r has no files, skipped' % name)
            continue
        log('group %d/%d %r: %d channel(s)'
            % (gi, len(groups), name or 'all', len(files)))

        result = eng_stats.run_stats(files, params, verbose=verbose)
        results[name] = result
        for r in result.rows:
            r['area'] = name
        all_rows += result.rows

        header = 'all channels' if name in ('', 'all') else 'area: %s' % name
        sections.append(header + '\n' + '-' * len(header) + '\n'
                        + eng_stats.build_report(result))

    if not all_rows:
        raise RuntimeError('no clusters found on any channel')

    csv_path = _write_csv(all_rows, os.path.join(out_dir, cfg.get(
        'csv_name', 'chirp_cluster_stats.csv')),
        ['area'] + eng_video.STAT_FIELDS)
    log('table -> %s' % csv_path)

    report = ('CHIRP survey - %s\n%s\n\n' % (cfg.get('label', ''), '=' * 60)
              + '\n\n'.join(sections)
              + '\n\nwritten:\n  %s\n' % csv_path)
    report_path = os.path.abspath(os.path.join(
        out_dir, cfg.get('report_name', 'chirp_report.txt')))
    with open(report_path, 'w', encoding='utf-8') as fh:
        fh.write(report)
    log('report -> %s' % report_path)

    # Optional render, of the channels the survey liked best. Rendering every
    # channel of a 64-site probe is minutes of ffmpeg for a quick check, so the
    # default is none and video_top caps it.
    videos = []
    if cfg.get('video'):
        top = int(cfg.get('video_top', 0))
        by_stem = {Path(f).stem: Path(f)
                   for g in groups for f in g.get('files', [])}
        for name, result in results.items():
            ranked = _rank_channels(result)
            for stem in (ranked[:top] if top > 0 else ranked):
                path = by_stem.get(stem)
                if path is None:
                    continue
                log('rendering %s' % stem)
                # The survey already chose this channel's best window, so the
                # render reuses it rather than scanning the recording again.
                mp4 = eng_video.render(
                    path, Path(out_dir), params.duration,
                    result.segments[stem][0], params.band, params.fs,
                    cfg.get('fps', 30), tuple(cfg.get('size', (1600, 900))),
                    100, 1.0, cfg.get('crf', 20), 16, 1.0,
                    cfg.get('ylim', 100.0), params.neg_k, params.pos_k,
                    params.pre_ms, params.post_ms, params.refractory_ms,
                    0.26, False, max_k=params.max_k)
                videos.append(os.path.abspath(str(mp4)))

    summary = _summarise(all_rows)
    return _emit(dict(status='ok', csv=csv_path, report=report_path,
                      videos=videos, n_rows=len(all_rows),
                      n_channels=summary['n_channels'],
                      groups=[g.get('name') or '' for g in groups],
                      isolated=summary['isolated'], noise=summary['noise'],
                      chirp_version=version), cfg)


def _rank_channels(result):
    """Channel stems, best first: lowest auto_quality tag, then highest SNR."""
    best = {}
    for r in result.rows:
        tag = r.get('auto_quality') or 9      # 0 = not assessed, sorts last
        key = (tag, -r.get('snr', 0.0))
        if r['channel'] not in best or key < best[r['channel']]:
            best[r['channel']] = key
    return [c for c, _ in sorted(best.items(), key=lambda kv: kv[1])]


def _summarise(rows):
    """Counts the wrapper prints, so MATLAB does not have to parse the CSV."""
    best = {}
    for r in rows:
        tag = r.get('auto_quality') or 9
        best[r['channel']] = min(best.get(r['channel'], 9), tag)
    return dict(n_channels=len(best),
                isolated=sum(1 for v in best.values() if v == 1),
                noise=sum(1 for v in best.values() if v == 3))


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('RESULT ' + json.dumps({'status': 'error',
                                      'message': 'no config path given'}))
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
        _emit({'status': 'error', 'message': str(exc)}, cfg)   # so MATLAB sees why
        sys.exit(1)
