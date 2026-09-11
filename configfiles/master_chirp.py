"""Master pipeline: survey a session's raw INTAN channels with CHIRP.

Usage:   python master_chirp.py <config.json>

The config is normally written by the MATLAB wrapper (functions/intan/process_chirp.m)
but can be hand-written. It locates the CHIRP engine in ../toolboxes/CHIRP (or an
explicit "toolboxes" path in the config), surveys every channel group, writes the
tables and a text report, and hands back a small RESULT json with the output paths
and a summary.

Two modes ("mode" in the config):
  fast  the best few windows per channel at one threshold (the default)
  deep  the whole recording, a -4/-5/-6 sigma threshold sweep on the same
        windows, per-channel pooled clusters, and a suggested threshold per area

Channels are surveyed in groups because the cross-channel sharing measure is
relative to the size of the group it is computed over: CHIRP calls a waveform
"noise" when it appears on more than max(8, 0.20 x N) of the N channels it saw.
With a multi-area probe the pipeline passes one group per area, so an area is
judged against its own channel count rather than against the whole headstage.
"""
import os
import sys
import json
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

# config key -> StatsParams field. Absent keys fall back to the mode's preset.
_PARAM_KEYS = ('duration', 'start', 'fs', 'neg_k', 'pos_k', 'artifact_k',
               'max_k', 'segments', 'scan_step', 'sampling', 'share_windows')


def _emit(result, cfg):
    """Hand the machine-readable result back. The MATLAB wrapper reads it from
    the result file, so it is only printed when there is no file to read -
    a hand-run - where it would otherwise be lost."""
    path = (cfg or {}).get('result_json')
    if path:
        try:
            with open(path, 'w') as fh:
                json.dump(result, fh)
            return result
        except OSError:
            pass
    print('RESULT ' + json.dumps(result), flush=True)
    return result


def main(config_path):
    with open(config_path, 'r') as f:
        cfg = json.load(f)
    verbose = cfg.get('verbose', True)

    def log(msg):
        if verbose:
            print('[chirp] ' + msg, flush=True)

    here = os.path.dirname(os.path.abspath(__file__))
    root = os.path.dirname(here)
    toolboxes = cfg.get('toolboxes') or os.path.join(root, 'toolboxes', 'CHIRP')
    sys.path.insert(0, toolboxes)
    import dat_to_stats as eng_stats
    import dat_to_video as eng_video

    # The wrapper needs the NGL fork (parallel runs, deep mode); an upstream
    # copy vendored over it lacks run_deep. Fail here with a clear message
    # rather than on a missing attribute halfway through a session.
    version = getattr(eng_stats, '__version__', '0')
    if not hasattr(eng_stats, 'run_deep'):
        raise RuntimeError('CHIRP %s in %s is not the NGL fork (1.3.0+ngl.2 or '
                           'newer); it was probably re-vendored from upstream'
                           % (version, toolboxes))

    groups = [g for g in (cfg.get('groups') or []) if g.get('files')]
    if not groups:
        raise ValueError('no channel groups given')
    for g in groups:
        missing = [f for f in g['files'] if not os.path.isfile(f)]
        if missing:
            raise FileNotFoundError('missing .dat file(s): %s' % ', '.join(missing[:5]))
    out_dir = cfg.get('out_dir') or os.getcwd()

    mode = cfg.get('mode') or 'fast'
    given = {k: cfg.get(k) for k in _PARAM_KEYS}
    given['band'] = tuple(cfg['band']) if cfg.get('band') else None
    given['scan_range'] = tuple(cfg['scan_range']) if cfg.get('scan_range') else None
    params, ignored = eng_stats.params_for_mode(mode, **given)

    # One pool for the whole session: every area group reuses it, so the
    # workers' start-up (~2 s, mostly each one importing scipy) is paid once.
    all_files = [f for g in groups for f in g['files']]
    n_total = len(all_files)
    n_windows = eng_stats.planned_windows(all_files, params)
    workers = max(1, min(int(cfg.get('workers') or
                             eng_stats.default_workers(n_total, n_windows)),
                         n_total, eng_stats.MAX_WORKERS))
    n_work, per_s, serial_s = eng_stats.estimate_work(
        [(g.get('name') or '', g['files']) for g in groups], params, mode)
    estimate = dict(analyses=n_work, window_ms=round(per_s * 1000, 1),
                    serial_s=round(serial_s, 1))
    log('CHIRP %s, %s mode: %d channel(s) in %d group(s), %d worker(s); '
        '~%s of work if run serially'
        % (version, mode, n_total, len(groups), workers, eng_stats.fmt_seconds(serial_s)))
    for name in ignored:
        log('note: %s is set but ignored in deep mode (the threshold is swept)' % name)
    if cfg.get('estimate_only'):
        return _emit(dict(status='ok', mode=mode, estimate_only=True,
                          estimate=estimate, workers=workers,
                          n_channels=n_total, chirp_version=version), cfg)
    os.makedirs(out_dir, exist_ok=True)

    label = cfg.get('label', '')
    t_survey = time.perf_counter()
    pool = ProcessPoolExecutor(max_workers=workers) if workers > 1 else None
    try:
        if mode == 'deep':
            deep = eng_stats.run_deep(
                [(g.get('name') or '', g['files']) for g in groups], params,
                executor=pool, verbose=verbose)
            refs = {}
            for dr in deep:
                for r in dr.window_rows + dr.clusters + dr.sweep:
                    r['area'] = dr.name
                refs[dr.name] = eng_stats.reference_result(dr)
            work_s = deep[0].timings['work_s'] if deep else 0.0
        else:
            refs, work_s = {}, 0.0
            for g in groups:
                name = g.get('name') or ''
                result = eng_stats.run_stats(
                    g['files'], params, verbose=verbose, executor=pool,
                    label=name if len(groups) > 1 else '')
                for r in result.rows:
                    r.update(area=name, neg_k=params.neg_k)
                refs[name] = result
                work_s += result.timings.get('work_s', 0.0)
    finally:
        if pool is not None:
            pool.shutdown(wait=True)
    survey_s = time.perf_counter() - t_survey
    # Average workers busy, not speed-up over serial: see run_stats timings.
    busy = work_s / survey_s if survey_s > 0 else 1.0

    all_rows = [r for res in refs.values() for r in res.rows]
    if not all_rows:
        raise RuntimeError('no clusters found on any channel')

    if mode == 'deep':
        written = eng_stats.write_deep_outputs(deep, out_dir, label)
        sections = [_header(dr.name) + eng_stats.build_report(refs[dr.name])
                    + '\n\n' + eng_stats.build_deep_report(dr) for dr in deep]
    else:
        written = dict(windows=eng_stats.write_rows_csv(
            all_rows, os.path.join(out_dir, eng_stats.DEFAULT_CSV_NAME),
            eng_stats.WINDOW_FIELDS))
        sections = [_header(name) + eng_stats.build_report(res)
                    for name, res in refs.items()]

    report = ('CHIRP survey - %s\n%s\n' % (label, '=' * 60)
              + 'CHIRP %s, %s mode, %d channel(s): %.1f s on %d worker(s), '
                '%.1f busy on average\n\n' % (version, mode, n_total, survey_s,
                                              workers, busy)
              + '\n\n'.join(sections)
              + '\n\nwritten:\n' + ''.join('  %s\n' % p for p in written.values()))
    report_path = os.path.abspath(os.path.join(
        out_dir, cfg.get('report_name', eng_stats.DEFAULT_REPORT_NAME)))
    with open(report_path, 'w', encoding='utf-8') as fh:
        fh.write(report)

    # Optional render, of the channels the survey liked best. Rendering every
    # channel of a 64-site probe is minutes of ffmpeg for a quick check, so the
    # default is none and video_top caps it.
    videos = []
    if cfg.get('video'):
        top = int(cfg.get('video_top', 0))
        by_stem = {Path(f).stem: Path(f) for f in all_files}
        for name, result in refs.items():
            ranked = _rank_channels(result)
            for stem in (ranked[:top] if top > 0 else ranked):
                path = by_stem.get(stem)
                if path is None:
                    continue
                # The survey already chose this channel's best window, so the
                # render reuses it rather than scanning the recording again.
                # It reports itself on one self-updating line.
                mp4 = eng_video.render(
                    path, Path(out_dir), params.duration,
                    result.segments[stem][0], params.band, params.fs,
                    cfg.get('fps', 30), tuple(cfg.get('size', (1600, 900))),
                    100, 1.0, cfg.get('crf', 20), 16, 1.0,
                    cfg.get('ylim', 100.0), params.neg_k, params.pos_k,
                    params.pre_ms, params.post_ms, params.refractory_ms,
                    0.26, False, max_k=params.max_k, verbose=verbose)
                videos.append(os.path.abspath(str(mp4)))

    # Deep mode summarises the pooled clusters at the reference threshold;
    # fast mode, its per-window rows.
    summary = _summarise([r for dr in deep for r in dr.clusters
                          if r['neg_k'] == eng_stats.REFERENCE_K]
                         if mode == 'deep' else all_rows)
    result = dict(status='ok', mode=mode, report=report_path,
                  csv=str(written['windows']), videos=videos,
                  n_rows=len(all_rows), n_channels=summary['n_channels'],
                  groups=[g.get('name') or '' for g in groups],
                  isolated=summary['isolated'], noise=summary['noise'],
                  chirp_version=version, estimate=estimate,
                  timings=dict(survey_s=round(survey_s, 2), work_s=round(work_s, 2),
                               workers=workers, busy=round(busy, 2)))
    if mode == 'deep':
        result.update(clusters_csv=str(written['clusters']),
                      sweep_csv=str(written['sweep']),
                      suggestion_json=str(written['suggestion']),
                      suggested_neg_k={dr.name: dr.suggestion['suggested_neg_k']
                                       for dr in deep})
    return _emit(result, cfg)


def _header(name):
    title = 'all channels' if name in ('', 'all') else 'area: %s' % name
    return title + '\n' + '-' * len(title) + '\n'


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
