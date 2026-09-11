#!/usr/bin/env python3
"""
Cluster statistics for a whole recording, with nothing rendered.

Sizing a session up is CHIRP's main use: run every channel, read the report,
then choose the sorter's settings. That survey is what this module does, and
it is deliberately separate from dat_to_video.py, which exists to draw a single
excerpt. Nothing here opens a window, writes a WAV or needs ffmpeg.

The detection, clustering and per-cluster measurements all live in
dat_to_video.py; this module owns the run, in one of two modes:

  fast  one cross-channel pass, then each channel's best few windows
        (run_stats). The survey as it was in CHIRP 1.3.0.
  deep  the whole recording scanned, one window per equal slice of it, and the
        detection threshold swept over -4/-5/-6 sigma on those same windows;
        each channel is re-clustered on its spikes pooled across windows, and
        the sweep yields a suggested threshold per area (run_deep).

Usage:
    python dat_to_stats.py --folder /path/to/recording --all
    python dat_to_stats.py --folder /path/to/recording --all --mode deep

Importable too, which is how the NGL pipeline drives it (configfiles/
master_chirp.py):

    import dat_to_stats as eng_stats
    params, _ = eng_stats.params_for_mode("fast", fs=30000)
    result = eng_stats.run_stats(paths, params,
                                 workers=eng_stats.default_workers(len(paths)))
    print(eng_stats.build_report(result, []))

NGL fork: this copy has diverged from upstream CHIRP 1.3.0 (see _version.py).
Channels can run on a process pool; rows are identical to a serial run, in the
same order.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import sys
import time
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict, dataclass, field, replace
from pathlib import Path

import numpy as np

# Same nudge dat_to_video.py and intan_gui.py use: the modules import each
# other by bare name, so their own directory has to be importable even when the
# script is launched from a recording folder somewhere else entirely.
sys.path.insert(0, str(Path(__file__).resolve().parent))

import dat_to_audio as eng_audio
import dat_to_video as eng_video
from dat_to_audio import DEFAULT_BAND, SAMPLE_RATE, __version__
from dat_to_video import AUTO_LABEL, Cancelled, Progress, fmt_seconds

# Statistics are averaged over this many of the best-scoring windows per
# channel rather than a single one - one excerpt is a thin basis for a firing
# rate. When a video is rendered as well, it uses the first of them.
STAT_SEGMENTS = 3

DEFAULT_CSV_NAME = "chirp_cluster_stats.csv"
DEFAULT_REPORT_NAME = "chirp_report.txt"
CLUSTERS_CSV_NAME = "chirp_channel_clusters.csv"
SWEEP_CSV_NAME = "chirp_threshold_sweep.csv"
SUGGESTION_NAME = "chirp_suggestion.json"

# Deep mode sweeps the detection threshold over fixed values rather than taking
# one from the user. Windows are chosen once, at the reference, so every
# threshold is judged on identical data.
SWEEP_K = (4.0, 5.0, 6.0)
REFERENCE_K = 5.0

# Most of its spikes a channel's smallest real unit may lose below threshold
# before a threshold counts as too strict. Matches maxPercSpikesMissing (25%)
# in configfiles/bombcellConfig.m, the QC the sorted units face later.
MAX_MISSING = 0.25

# Settings each mode fills in when the caller leaves them unset. A scan_step
# of None means "the window length": the scan then covers every second of the
# recording exactly once, with no two candidates overlapping.
PRESETS = {
    "fast": dict(segments=STAT_SEGMENTS, scan_step=60.0, sampling="best",
                 share_windows=1),
    "deep": dict(segments=12, scan_step=None, sampling="stratified",
                 share_windows=5),
}

# ProcessPoolExecutor refuses more than 61 workers on Windows.
MAX_WORKERS = 61 if os.name == "nt" else 256

# Below this many windows a pool costs more than it saves. Measured: start-up
# ~2 s (every worker imports scipy), then ~1/3 of the serial time, so parallel
# breaks even near 3 s of serial work - about 200 windows at ~17 ms each.
MIN_PARALLEL_WINDOWS = 250


# ------------------------------------------------------------------- params --
@dataclass
class StatsParams:
    """
    Everything that changes what the numbers mean, in the engine's own
    vocabulary. Callers build one of these rather than passing a dozen loose
    arguments, so the report can state the settings the rows were measured
    under without being handed a second copy of them.
    """
    duration: float = 10.0              # excerpt length, s
    start: float | None = None          # fixed start, s; None = scan for windows
    band: tuple[float, float] = DEFAULT_BAND
    fs: int = SAMPLE_RATE               # acquisition rate, Hz
    neg_k: float = 5.0                  # detection threshold, sigma below zero
    pos_k: float = 8.0                  # artifact rejection, sigma above zero
    artifact_k: float = 18.0            # window is "clean" below this x sigma
    max_k: int = 3                      # most amplitude clusters to consider
    segments: int = STAT_SEGMENTS       # windows per channel
    scan_step: float = 60.0             # spacing of candidate windows, s
    scan_range: tuple[float, float] | None = None
    pre_ms: float = 1.0                 # waveform window before the trough
    post_ms: float = 2.0                # waveform window after the trough
    refractory_ms: float = 1.0          # minimum spacing between detections
    sampling: str = "best"              # "best" top-N, or "stratified" by slice
    share_windows: int = 1              # common windows for the sharing pass

    def validate(self) -> "StatsParams":
        """Reject the combinations that would otherwise fail deep in a scan."""
        if self.duration <= 0:
            raise ValueError("duration must be positive")
        if not 1 <= self.max_k <= 3:
            raise ValueError("max_k must be 1, 2 or 3")
        if self.segments < 1:
            raise ValueError("segments must be at least 1")
        if self.artifact_k <= 0:
            raise ValueError("artifact_k must be positive")
        if self.scan_step <= 0:
            raise ValueError("scan_step must be positive")
        if self.sampling not in ("best", "stratified"):
            raise ValueError("sampling must be 'best' or 'stratified'")
        if self.share_windows < 1:
            raise ValueError("share_windows must be at least 1")
        lo, hi = self.band
        if not 0 < lo < hi < self.fs / 2:
            raise ValueError(f"band must satisfy 0 < low < high < {self.fs / 2:.0f} Hz "
                             f"(Nyquist), got {lo:g}-{hi:g} Hz")
        return self


def params_for_mode(mode, **given):
    """
    StatsParams for a survey mode. Every value in `given` that is not None
    overrides the mode's preset. Deep mode sweeps the threshold itself and
    chooses its own windows, so a given neg_k or start is dropped there.

    Returns (params, names of the given settings that were ignored).
    """
    if mode not in PRESETS:
        raise ValueError(f"mode must be one of {sorted(PRESETS)}, got {mode!r}")
    settings = {k: v for k, v in given.items() if v is not None}
    ignored = []
    if mode == "deep":
        for name in ("neg_k", "start"):
            if settings.pop(name, None) is not None:
                ignored.append(name)
        settings["neg_k"] = REFERENCE_K
    merged = {**PRESETS[mode], **settings}
    if merged.get("scan_step") is None:
        merged["scan_step"] = merged.get("duration", StatsParams.duration)
    return StatsParams(**merged).validate(), ignored


def planned_windows(paths, params):
    """Windows a survey will analyse: every scan candidate, plus the measured ones."""
    n = 0
    for p in paths:
        if params.start is None or params.sampling == "stratified":
            total = eng_audio.file_duration_s(Path(p), params.fs)
            lo, hi = params.scan_range or (0.0, total - params.duration)
            hi = min(hi, total - params.duration)
            n += max(1, int(max(0.0, hi - lo) // params.scan_step) + 1)
        n += params.segments
    return n


def estimate_work(groups, params, mode, probe=5):
    """
    Serial work a survey will take on this machine and this data drive: a few
    windows of the first channel are timed and multiplied by the analyses the
    run will make. Parallel wall time is lower by an amount that depends on
    memory bandwidth and the drive (measured ~3x on 8 cores), so only the
    serial figure is claimed; the live ETA covers the rest.

    groups: list of (name, paths). Returns (analyses, s per analysis, serial s).
    """
    files = [Path(f) for _, paths in groups for f in paths]
    if not files:
        return 0, 0.0, 0.0
    total = eng_audio.file_duration_s(files[0], params.fs)
    t0 = time.perf_counter()
    for start in np.linspace(0.0, max(0.0, total - params.duration), probe):
        x = eng_audio.bandpass(
            eng_audio.read_segment(files[0], float(start), params.duration, params.fs),
            params.band[0], params.band[1], params.fs)
        eng_video.analyse(x, params.fs, params.neg_k, params.pos_k, params.pre_ms,
                          params.post_ms, params.refractory_ms, max_k=params.max_k)
    per = (time.perf_counter() - t0) / probe
    per_k = len(SWEEP_K) if mode == "deep" else 1
    n = planned_windows(files, params) + params.segments * (per_k - 1) * len(files)
    n += sum(len(p) * params.share_windows * per_k for _, p in groups if len(p) >= 2)
    return n, per, n * per


def default_workers(n_channels, n_windows=None):
    """
    Worker processes to use when the caller does not say.

    The per-window work is large numpy array operations, which are bound by
    memory bandwidth rather than arithmetic: measured on an 8-core / 16-thread
    machine, speed-up plateaued at the physical core count (~3x) and more
    workers bought nothing. Half the logical CPUs approximates the physical
    cores without a psutil dependency. Small jobs (see MIN_PARALLEL_WINDOWS)
    run serially, where a pool would only add its start-up.
    """
    if n_windows is not None and n_windows < MIN_PARALLEL_WINDOWS:
        return 1
    return max(1, min((os.cpu_count() or 2) // 2, int(n_channels), MAX_WORKERS))


# ------------------------------------------------------------------ results --
@dataclass
class StatsResult:
    """A finished fast run: the rows, plus the context needed to describe them."""
    rows: list = field(default_factory=list)
    segments: dict = field(default_factory=dict)   # channel stem -> [start_s, ...]
    shares: dict = field(default_factory=dict)     # channel stem -> co-active count
    n_channels: int = 0                            # channels in the sharing pass
    params: StatsParams = field(default_factory=StatsParams)
    timings: dict = field(default_factory=dict)    # wall/work seconds, workers


@dataclass
class DeepResult:
    """A finished deep run for one channel group (one area)."""
    name: str = ""
    params: StatsParams = field(default_factory=StatsParams)
    window_rows: list = field(default_factory=list)  # window x threshold x cluster
    clusters: list = field(default_factory=list)     # pooled: channel x threshold x cluster
    sweep: list = field(default_factory=list)        # channel x threshold
    segments: dict = field(default_factory=dict)     # channel stem -> window starts
    shares: dict = field(default_factory=dict)       # threshold -> stem -> count
    n_channels: int = 0
    suggestion: dict = field(default_factory=dict)
    notes: list = field(default_factory=list)
    timings: dict = field(default_factory=dict)      # the whole run, all groups


# ---------------------------------------------------------------- the pool --
def _pool_size(workers, executor, n_jobs):
    if executor is not None:
        return getattr(executor, "_max_workers", 1)
    if workers and workers > 1 and n_jobs > 1:
        return max(1, min(int(workers), n_jobs, MAX_WORKERS))
    return 1


def _run_jobs(pool, fn, arglist, on_done=None, cancel=None):
    """
    fn(*args) for every args in `arglist`, serially when `pool` is None, else
    on the pool. Results come back in input order; on_done(i, result) fires as
    each finishes, which in parallel is completion order. cancel is polled
    between jobs: a worker cannot be interrupted mid-job, so whatever has not
    started is dropped and in-flight jobs finish on their own.
    """
    out = [None] * len(arglist)
    if pool is None:
        for i, args in enumerate(arglist):
            if cancel is not None and cancel():
                raise Cancelled("cancelled")
            out[i] = fn(*args)
            if on_done is not None:
                on_done(i, out[i])
        return out
    futures = {}
    try:
        futures = {pool.submit(fn, *args): i for i, args in enumerate(arglist)}
        for fut in as_completed(futures):
            i = futures[fut]
            out[i] = fut.result()
            if on_done is not None:
                on_done(i, out[i])
            if cancel is not None and cancel():
                raise Cancelled("cancelled")
    except BaseException:
        for fut in futures:
            fut.cancel()
        raise
    return out


def _counter(prog, progress=None):
    """on_done callback that advances a Progress (and an optional 0..1 hook)."""
    done = [0]

    def on_done(_i, _result):
        done[0] += 1
        prog.update(done[0])
        if progress is not None:
            progress(done[0] / prog.total)
    return on_done


# ---------------------------------------------------------------- windows --
def _choose_windows(path, params):
    """
    Start times of the windows to measure on one channel.

    "best": the top `segments` windows by score, anywhere in the recording
    (the 1.3.0 behaviour). "stratified": the recording cut into `segments`
    equal slices and the best window taken inside each, so the windows cover
    the whole session instead of bunching around its cleanest stretch.
    """
    path = Path(path)
    p = params

    def best(search_range, top_n):
        return eng_video.find_best_windows(
            path, p.duration, p.band, p.fs, p.scan_step, p.artifact_k,
            p.neg_k, p.pos_k, p.pre_ms, p.post_ms, p.refractory_ms,
            search_range, True, False, max_k=p.max_k, top_n=top_n)

    if p.sampling == "best":
        segments = [p.start] if p.start is not None else best(p.scan_range, p.segments)
    else:
        total = eng_audio.file_duration_s(path, p.fs)
        edges = np.linspace(0.0, total, p.segments + 1)
        segments = []
        for lo, hi in zip(edges[:-1], edges[1:]):
            if hi - lo >= p.duration:
                # find_best_windows takes the latest allowed START as its upper
                # bound, and excludes it; the epsilon keeps a window flush
                # with the slice's end as a candidate.
                segments += best((float(lo), float(hi) - p.duration + 1e-6), 1)
    if not segments:
        segments = [(eng_audio.file_duration_s(path, p.fs) - p.duration) / 2]
    return segments


def _sharing_starts(shortest, params):
    """Common windows for the sharing pass: the middle one, or one per slice."""
    d = params.duration
    if params.share_windows <= 1:
        return [params.start if params.start is not None
                else max(0.0, (shortest - d) / 2)]
    k = params.share_windows
    return [min(max(0.0, (i + 0.5) * shortest / k - d / 2), shortest - d)
            for i in range(k)]


def _share_job(paths, start, params, neg_k):
    """One sharing window at one threshold. Module-level so a pool can run it."""
    shares, _ = eng_video.session_sharing(
        [Path(p) for p in paths], start, params.duration, params.band, params.fs,
        neg_k, params.pos_k, pre_ms=params.pre_ms, post_ms=params.post_ms,
        refractory_ms=params.refractory_ms, verbose=False)
    return shares


def _median_shares(per_window):
    """Each channel's co-active count, as the median over sharing windows."""
    return {stem: float(np.median([w[stem] for w in per_window]))
            for stem in per_window[0]}


# ------------------------------------------------------------ fast: channel --
def _channel_job(path, params, share_count, n_channels):
    """
    One channel's share of a fast survey: choose its windows, measure each.
    Module-level so a process pool can pickle it by reference.

    Returns (segments, rows per segment, wall seconds).
    """
    t0 = time.perf_counter()
    path = Path(path)
    segments = _choose_windows(path, params)
    per_segment = []
    for rank, seg_start in enumerate(segments, start=1):
        rows, _ = eng_video.channel_stats(
            path, seg_start, params.duration, params.band, params.fs,
            params.neg_k, params.pos_k,
            pre_ms=params.pre_ms, post_ms=params.post_ms,
            refractory_ms=params.refractory_ms, max_k=params.max_k,
            segment=rank, share_count=share_count, n_channels=n_channels)
        per_segment.append(rows)
    return segments, per_segment, time.perf_counter() - t0


# --------------------------------------------------------------- fast: run --
def run_stats(paths, params, cancel=None, progress=None, on_rows=None,
              on_channel=None, verbose=True, workers=1, executor=None,
              label="") -> StatsResult:
    """
    Fast survey of every channel in `paths`; returns the rows.

    paths       list of amp-*.dat Paths, all from one session.
    params      StatsParams (see params_for_mode).
    cancel      callable() -> bool, polled between channels; True raises Cancelled.
    progress    callable(fraction 0..1) as the run advances.
    on_rows     callable(rows) once per channel x segment.
    on_channel  callable(path, segments, rows) once per finished channel. The
                chosen windows are handed over so a caller that also renders
                does not have to scan for them a second time.
    verbose     one self-updating progress line per phase on stdout.
    workers     processes for the per-channel pass; 1 runs serially in-process.
    executor    an existing ProcessPoolExecutor to reuse instead of starting
                one, so a caller surveying several groups pays start-up once.
    label       prefix for the progress lines, e.g. the area name.

    The cross-channel pass runs first over a window common to every channel:
    sharing is a property of a channel within its session, so each channel's
    own best window would not be comparable. It needs two or more channels
    and is skipped below that, leaving auto_quality on waveform and rate alone.

    In parallel, on_rows / on_channel fire in completion order; result.rows is
    always in input channel order, identical to a serial run.
    """
    t_run = time.perf_counter()
    params = params.validate()
    paths = [Path(p) for p in paths]
    if not paths:
        raise ValueError("no channels given")
    prefix = f"{label}: " if label else ""
    result = StatsResult(params=params)
    share_weight = 0.05 if len(paths) >= 2 else 0.0

    def advance(frac):
        if progress is not None:
            progress(min(1.0, max(0.0, frac)))

    n_workers = _pool_size(workers, executor, len(paths))
    pool = executor or (ProcessPoolExecutor(max_workers=n_workers)
                        if n_workers > 1 else None)
    try:
        # --- cross-channel pass ------------------------------------------
        if len(paths) >= 2:
            shortest = min(eng_audio.file_duration_s(p, params.fs) for p in paths)
            starts = _sharing_starts(shortest, params)
            if len(starts) == 1:
                prog = Progress(prefix + "sharing", len(paths), verbose)

                def share_progress(f):
                    prog.update(round(f * len(paths)))
                    advance(share_weight * f)
                result.shares, result.n_channels = eng_video.session_sharing(
                    paths, starts[0], params.duration, params.band, params.fs,
                    params.neg_k, params.pos_k,
                    pre_ms=params.pre_ms, post_ms=params.post_ms,
                    refractory_ms=params.refractory_ms,
                    cancel=cancel, progress=share_progress, verbose=False)
            else:
                prog = Progress(prefix + "sharing", len(starts), verbose)
                per_window = _run_jobs(
                    pool, _share_job,
                    [([str(p) for p in paths], s, params, params.neg_k)
                     for s in starts],
                    _counter(prog), cancel)
                result.shares, result.n_channels = (_median_shares(per_window),
                                                    len(paths))
            shared = sum(eng_video.shared_across_batch(s, result.n_channels)
                         for s in result.shares.values())
            prog.close(f", {shared}/{result.n_channels} channel(s) shared "
                       f"across a large batch")
        t_share = time.perf_counter() - t_run

        # --- per-channel pass --------------------------------------------
        prog = Progress(prefix + "channels", len(paths), verbose)
        tick = _counter(prog, lambda f: advance(share_weight
                                                + (1 - share_weight) * f))

        def finish(i, out):
            tick(i, out)
            segments, per_segment, _ = out
            result.segments[paths[i].stem] = segments
            if on_rows is not None:
                for rows in per_segment:
                    on_rows(rows)
            if on_channel is not None:
                on_channel(paths[i], segments,
                           [r for rows in per_segment for r in rows])

        t_chan = time.perf_counter()
        outputs = _run_jobs(
            pool, _channel_job,
            [(str(p), params, result.shares.get(p.stem), result.n_channels)
             for p in paths],
            finish, cancel)
        prog.close()
    finally:
        if pool is not None and executor is None:
            pool.shutdown(wait=True)

    result.segments = {paths[i].stem: out[0] for i, out in enumerate(outputs)}
    result.rows = [r for out in outputs for rows in out[1] for r in rows]
    result.timings = _timings(t_run, t_share, t_chan, [o[2] for o in outputs],
                              n_workers)
    advance(1.0)
    return result


def _timings(t_run, t_share, t_chan, job_seconds, n_workers):
    # work_s sums each channel's own wall time. Under parallel load every
    # channel runs slower (memory bandwidth is shared), so work_s / wall is how
    # many workers were busy on average - NOT the speed-up over a serial run,
    # which only a serial run can measure.
    t_end = time.perf_counter()
    work = float(sum(job_seconds))
    return dict(total_s=t_end - t_run, sharing_s=t_share,
                channels_s=t_end - t_chan, work_s=work, workers=max(1, n_workers),
                busy=work / (t_end - t_chan) if t_end > t_chan else 1.0)


# ------------------------------------------------------------ deep: channel --
def _none_if_nan(x):
    return None if x is None or (isinstance(x, float) and math.isnan(x)) else x


def _pool_clusters(channel, k, acc, n_windows, p, share_count, n_channels):
    """
    Re-cluster one channel's spikes pooled over all its windows, at threshold
    k: one clustering per channel instead of one per window, so a unit keeps
    one identity across the session. Drift smears a unit's amplitudes across
    windows; that risk is accepted, and amplitude_cv_windows records it.

    Returns (pooled cluster rows, the channel's sweep row at k).
    """
    n_pre = int(round(p.pre_ms * p.fs / 1000))
    width = n_pre + int(round(p.post_ms * p.fs / 1000))
    amp = np.concatenate(acc["amp"]) if acc["amp"] else np.zeros(0)
    filled = [w for w in acc["waves"] if len(w)]
    waves = np.vstack(filled) if filled else np.zeros((0, width))
    win = np.concatenate(acc["win"]) if acc["win"] else np.zeros(0, dtype=int)
    sig_w = np.asarray(acc["sigma"], dtype=float)
    # Every spike is judged against its own window's threshold, k x sigma.
    thr = k * sig_w[win - 1] if len(win) else np.zeros(0)
    sigma = float(np.median(sig_w)) if len(sig_w) else 0.0
    total = n_windows * p.duration

    if len(amp):
        labels, centres, _ = eng_video.cluster_amplitudes(amp, max_k=p.max_k)
    else:
        labels, centres = np.zeros(0, dtype=int), np.array([0.0])
    rows = eng_video.cluster_stats(
        dict(waves=waves, amp=amp, labels=labels, centres=centres, sigma=sigma),
        p.fs, total, p.pre_ms)

    share_frac = (share_count / n_channels
                  if share_count is not None and n_channels else None)
    for r in rows:
        sel = labels == r["cluster"] - 1
        tag = eng_video.auto_quality(share_count, n_channels, r["wf_residual"],
                                     r["firing_rate_sp_s"])
        rates, amps, tags = [], [], []
        for w in range(1, n_windows + 1):
            m = sel & (win == w)
            n_w = int(m.sum())
            rates.append(n_w / p.duration)
            if n_w:
                amps.append(float(amp[m].mean()))
                t_w = eng_video.auto_quality(
                    share_count, n_channels,
                    eng_video.waveform_residual(waves[m]), n_w / p.duration)
                if t_w:                                  # 0 = too few to judge
                    tags.append(t_w)
        mean_amp = float(np.mean(amps)) if amps else 0.0
        r.update(channel=channel, neg_k=k, n_clusters=len(centres),
                 n_windows=n_windows, n_windows_present=len(amps),
                 rate_sd_sp_s=float(np.std(rates)),
                 amplitude_cv_windows=(float(np.std(amps) / abs(mean_amp))
                                       if len(amps) >= 2 and mean_amp else None),
                 share_count=share_count, share_frac=share_frac,
                 auto_quality=tag,
                 tag_agreement=(sum(t == tag for t in tags) / len(tags)
                                if tags else None),
                 missing_frac=_none_if_nan(
                     eng_video.missing_fraction(amp[sel], thr[sel])))

    # "Real" units are the ones tagged isolated or multi-unit. Isolation alone
    # cannot drive the sweep: the residual cut-off needs a trough ~6.7x the
    # noise, so a unit near any of the swept thresholds never counts as
    # isolated however well it is detected.
    iso = [r for r in rows if r["auto_quality"] == 1]
    real = [r for r in rows if r["auto_quality"] in (1, 2)]
    tags = [r["auto_quality"] for r in rows]
    smallest = min(real, key=lambda r: abs(r["mean_amplitude_uV"])) if real else None
    sweep = dict(
        channel=channel, neg_k=k, n_windows=n_windows, n_clusters=len(rows),
        n_real=len(real), n_isolated=len(iso), n_multiunit=tags.count(2),
        n_noise=tags.count(3),
        noise_flag=int(eng_video.shared_across_batch(share_count, n_channels)),
        total_rate_sp_s=len(amp) / total if total else 0.0,
        n_rejected=acc["n_rej"], sigma_uV=sigma,
        share_count=share_count, share_frac=share_frac,
        missing_frac=smallest["missing_frac"] if smallest else None,
        residual_isolated=(float(np.median([r["wf_residual"] for r in iso]))
                           if iso else None))
    return rows, sweep


def _deep_channel_job(path, params, shares, n_channels):
    """
    One channel's share of a deep survey. Windows are chosen once, at the
    reference threshold; each window is read and filtered once and analysed at
    every threshold of the sweep. Module-level so a process pool can run it.

    Returns (segments, window rows, pooled cluster rows, sweep rows, seconds).
    """
    t0 = time.perf_counter()
    path = Path(path)
    p = params
    segments = _choose_windows(path, replace(p, neg_k=REFERENCE_K))
    window_rows = []
    pooled = {k: dict(amp=[], waves=[], win=[], sigma=[], n_rej=0) for k in SWEEP_K}
    for w, start in enumerate(segments, start=1):
        x = eng_video.load_excerpt(path, start, p.duration, p.band, p.fs)
        for k in SWEEP_K:
            a = eng_video.analyse(x, p.fs, k, p.pos_k, p.pre_ms, p.post_ms,
                                  p.refractory_ms, max_k=p.max_k)
            rows = eng_video.rows_from_analysis(
                a, path.stem, start, p.duration, p.fs, p.pre_ms, segment=w,
                share_count=shares.get(k), n_channels=n_channels)
            for r in rows:
                r["neg_k"] = k
            window_rows += rows
            acc = pooled[k]
            acc["amp"].append(a["amp"])
            acc["waves"].append(a["waves"])
            acc["win"].append(np.full(len(a["amp"]), w))
            acc["sigma"].append(a["sigma"])
            acc["n_rej"] += a["n_rej"]

    clusters, sweep = [], []
    for k in SWEEP_K:
        rows, row = _pool_clusters(path.stem, k, pooled[k], len(segments), p,
                                   shares.get(k), n_channels)
        clusters += rows
        sweep.append(row)

    # The channel's pick: the strictest threshold that keeps every real unit
    # the channel shows at its best, while its smallest unit loses at most
    # MAX_MISSING of its spikes below threshold. A stricter threshold admits
    # less noise, so it wins unless it costs a unit or too many spikes. If even
    # the loosest loses too many, the loosest (the least loss) is picked.
    # A threshold at which the channel trips the sharing test tags every
    # cluster 3 (no real units), so it can never be picked - the "no unit
    # turns into noise" condition holds by construction. An unknown loss
    # (under 20 spikes to fit) does not count against a threshold.
    top = max(row["n_real"] for row in sweep)
    keeps = [row for row in sweep if top and row["n_real"] == top]
    for row in keeps:
        row["loss_ok"] = int(row["missing_frac"] is None
                             or row["missing_frac"] <= MAX_MISSING)
    fine = [row["neg_k"] for row in keeps if row["loss_ok"]]
    preferred = (max(fine) if fine else min(row["neg_k"] for row in keeps)
                 if keeps else None)
    for row in sweep:
        row.setdefault("loss_ok", 0)
        row["preferred_neg_k"] = preferred
        row["is_preferred"] = int(row["neg_k"] == preferred)
    return segments, window_rows, clusters, sweep, time.perf_counter() - t0


def _kkey(k):
    return f"{k:g}"


def suggest_threshold(sweep):
    """
    The area's suggested threshold: the most common channel pick (ties go to
    the stricter threshold), with the vote and the evidence behind it.
    Channels with no real unit at any threshold do not vote.
    """
    picks = {r["channel"]: r["preferred_neg_k"] for r in sweep}
    votes = Counter(k for k in picks.values() if k is not None)
    n_voting = sum(votes.values())
    suggested, confidence = None, None
    if n_voting:
        top = max(votes.values())
        suggested = max(k for k, v in votes.items() if v == top)
        confidence = votes[suggested] / n_voting

    def median(vals):
        vals = [v for v in vals if v is not None]
        return float(np.median(vals)) if vals else None

    evidence = {}
    for k in SWEEP_K:
        at_k = [r for r in sweep if r["neg_k"] == k]
        evidence[_kkey(k)] = dict(
            channels_with_unit=sum(r["n_real"] > 0 for r in at_k),
            real_units=sum(r["n_real"] for r in at_k),
            isolated_clusters=sum(r["n_isolated"] for r in at_k),
            noise_channels=sum(r["noise_flag"] for r in at_k),
            median_missing_frac=median(r["missing_frac"] for r in at_k),
            median_residual=median(r["residual_isolated"] for r in at_k))
    return dict(suggested_neg_k=suggested, confidence=confidence,
                max_missing=MAX_MISSING,
                votes={_kkey(k): votes.get(k, 0) for k in SWEEP_K},
                n_voting=n_voting, n_channels=len(picks), evidence=evidence)


# --------------------------------------------------------------- deep: run --
def run_deep(groups, params, workers=1, executor=None, cancel=None,
             verbose=True):
    """
    Deep survey of one or more channel groups (one per area).

    groups      list of (name, paths). Sharing is measured within each group,
                so each area is judged against its own channel count.
    params      StatsParams from params_for_mode("deep", ...).

    Every group's jobs go into one queue - sharing first, then all channels -
    so no worker idles between areas. Returns a DeepResult per group, in order.
    """
    t_run = time.perf_counter()
    params = params.validate()
    plans = []
    for name, paths in groups:
        paths = [Path(p) for p in paths]
        if not paths:
            continue
        shortest = min(eng_audio.file_duration_s(p, params.fs) for p in paths)
        fit = int(shortest // params.duration)
        p, notes = params, []
        if fit < params.segments:
            p = replace(params, segments=max(1, fit))
            notes.append(f"recording holds only {fit} window(s) of "
                         f"{params.duration:g} s; measured {p.segments} instead "
                         f"of {params.segments}")
        plans.append(dict(name=name, paths=paths, params=p, notes=notes,
                          n=len(paths) if len(paths) >= 2 else 0,
                          starts=(_sharing_starts(shortest, p)
                                  if len(paths) >= 2 else [])))

    n_chan = sum(len(pl["paths"]) for pl in plans)
    n_workers = _pool_size(workers, executor, n_chan)
    pool = executor or (ProcessPoolExecutor(max_workers=n_workers)
                        if n_workers > 1 else None)
    try:
        # --- sharing: one job per (group, threshold, window) --------------
        keys, args = [], []
        for gi, pl in enumerate(plans):
            for k in SWEEP_K:
                for s in pl["starts"]:
                    keys.append((gi, k))
                    args.append(([str(x) for x in pl["paths"]], s, pl["params"], k))
        prog = Progress("sharing", len(args), verbose)
        per_job = _run_jobs(pool, _share_job, args, _counter(prog), cancel)
        grouped = {}
        for key, res in zip(keys, per_job):
            grouped.setdefault(key, []).append(res)
        shares = [{k: {} for k in SWEEP_K} for _ in plans]
        for (gi, k), windows in grouped.items():
            shares[gi][k] = _median_shares(windows)
        flagged = sum(eng_video.shared_across_batch(s, pl["n"])
                      for gi, pl in enumerate(plans)
                      for s in shares[gi][REFERENCE_K].values())
        prog.close(f", {flagged}/{n_chan} channel(s) shared across a large "
                   f"batch at -{REFERENCE_K:g} sigma")
        t_share = time.perf_counter() - t_run

        # --- channels: every channel of every group, one queue ------------
        owners, chan_args = [], []
        for gi, pl in enumerate(plans):
            for path in pl["paths"]:
                owners.append((gi, path))
                chan_args.append((str(path), pl["params"],
                                  {k: shares[gi][k].get(path.stem) for k in SWEEP_K},
                                  pl["n"]))
        prog = Progress("channels", len(chan_args), verbose)
        t_chan = time.perf_counter()
        outs = _run_jobs(pool, _deep_channel_job, chan_args, _counter(prog), cancel)
        prog.close()
    finally:
        if pool is not None and executor is None:
            pool.shutdown(wait=True)

    timings = _timings(t_run, t_share, t_chan, [o[4] for o in outs], n_workers)
    results = []
    for gi, pl in enumerate(plans):
        dr = DeepResult(name=pl["name"], params=pl["params"], shares=shares[gi],
                        n_channels=pl["n"], notes=pl["notes"], timings=timings)
        mine = [(path, out) for (g, path), out in zip(owners, outs) if g == gi]
        for path, (segments, windows, clusters, sweep, _) in mine:
            dr.segments[path.stem] = segments
            dr.window_rows += windows
            dr.clusters += clusters
            dr.sweep += sweep
        dr.suggestion = suggest_threshold(dr.sweep)
        results.append(dr)
    return results


# ------------------------------------------------------------------ report --
def build_report(result: StatsResult, files=()) -> str:
    """Plain-text summary of a whole run. `files` are the paths it wrote."""
    import statistics as st

    rows, p = result.rows, result.params
    if not rows:
        return "CHIRP - cluster statistics report\nno clusters found"

    # Cluster count is a property of one excerpt, so tally it per
    # channel x segment rather than collapsing it onto the channel.
    per_seg = {}
    for r in rows:
        per_seg[(r["channel"], r["segment"])] = r["n_clusters"]
    kdist = Counter(per_seg.values())
    channels = {r["channel"] for r in rows}

    def pm(key, fmt):
        vals = [r[key] for r in rows]
        m = st.fmean(vals)
        sd = st.pstdev(vals) if len(vals) > 1 else 0.0
        return f"{fmt.format(m)} +/- {fmt.format(sd)}"

    if p.sampling == "stratified":
        how = f", {p.segments} windows, one per equal slice of the recording"
    elif p.start is not None:
        how = f" from {p.start:.0f} s"
    else:
        how = f", best {p.segments} non-overlapping windows per channel"

    L = []
    L.append("CHIRP - cluster statistics report")
    L.append("=" * 60)
    L.append(f"channels analysed   : {len(channels)}")
    L.append(f"segments analysed   : {len(per_seg)} "
             f"({len(per_seg) / max(1, len(channels)):.1f} per channel)")
    L.append(f"excerpt             : {p.duration:.0f} s" + how)
    L.append(f"band-pass           : {p.band[0]:.0f}-{p.band[1]:.0f} Hz")
    L.append(f"detect / reject     : -{p.neg_k:g} sigma / +{p.pos_k:g} sigma")
    L.append(f"artifact scan       : {p.artifact_k:g} sigma")
    L.append(f"max clusters        : {p.max_k}")
    t = result.timings
    if t:
        L.append(f"time                : {t['total_s']:.1f} s (sharing "
                 f"{t['sharing_s']:.1f} s, channels {t['channels_s']:.1f} s on "
                 f"{t['workers']} worker(s), {t['busy']:.1f} busy on average)")
    L.append("")
    L.append("clusters per segment : " + ", ".join(
        f"{k} cluster(s) x {v} segment(s)" for k, v in sorted(kdist.items())))
    L.append(f"cluster rows in total: {len(rows)}")
    L.append(f"spikes in total      : {sum(r['n_spikes'] for r in rows)}")
    L.append("")
    L.append("across all cluster rows (mean +/- sd):")
    L.append(f"  mean amplitude     : {pm('mean_amplitude_uV', '{:.1f}')} uV")
    L.append(f"  firing rate        : {pm('firing_rate_sp_s', '{:.1f}')} sp/s")
    L.append(f"  SNR                : {pm('snr', '{:.1f}')}")
    L.append(f"  half-width         : {pm('half_width_ms', '{:.3f}')} ms")
    L.append(f"  trough-to-peak     : {pm('trough_to_peak_ms', '{:.2f}')} ms")
    L.append(f"  noise sigma        : {pm('sigma_uV', '{:.2f}')} uV")
    L.append("")
    tagged = [r for r in rows if r.get("auto_quality")]
    if tagged:
        ac = Counter(r["auto_quality"] for r in tagged)
        L.append("first-pass tag, per cluster:")
        for v in (1, 2, 3):
            if ac.get(v):
                L.append(f"  {AUTO_LABEL[v]:<14s} {ac[v]:4d} / {len(tagged)}")
        if len(tagged) < len(rows):
            L.append(f"  {'not assessed':<14s} {len(rows) - len(tagged):4d}"
                     f" / {len(rows)}  (too few spikes to judge)")
        # A channel is worth opening if any of its clusters looks isolated.
        best = {}
        for r in tagged:
            b = best.get(r["channel"], 9)
            best[r["channel"]] = min(b, r["auto_quality"])
        bc = Counter(best.values())
        L.append("")
        L.append("first-pass tag, per channel (best cluster on it):")
        for v in (1, 2, 3):
            if bc.get(v):
                L.append(f"  {AUTO_LABEL[v]:<14s} {bc[v]:4d} / {len(best)}")
        worth = sorted(c for c, v in best.items() if v == 1)
        if worth:
            L.append("")
            L.append("channels with a cluster that looks isolated:")
            for i in range(0, len(worth), 6):
                L.append("  " + ", ".join(worth[i:i + 6]))
        L.append("")
        L.append("This tag is a first pass meant to be reviewed, not a")
        L.append("verdict. Against one hand-tagged session it agreed on")
        L.append("95% of channels and 91% of clusters.")
    if files:
        L.append("")
        L.append("written:")
        for f in files:
            L.append(f"  {f}")
    L.append("")
    L.append("note: half-width and trough-to-peak shift with the")
    L.append("      band-pass, so compare them only between recordings")
    L.append("      filtered the same way.")
    return "\n".join(L)


def reference_result(dr: DeepResult) -> StatsResult:
    """A deep run seen as a fast one at the reference threshold, for build_report."""
    return StatsResult(rows=[r for r in dr.window_rows if r["neg_k"] == REFERENCE_K],
                       segments=dr.segments,
                       shares=dr.shares.get(REFERENCE_K, {}),
                       n_channels=dr.n_channels, params=dr.params)


def build_deep_report(dr: DeepResult) -> str:
    """The threshold sweep and its suggestion, for one group."""
    p, s = dr.params, dr.suggestion

    def fmt(x, spec):
        return "-" if x is None else format(x, spec)

    L = ["threshold sweep", "-" * 60,
         f"windows per channel : {p.segments} x {p.duration:g} s, one per equal "
         f"slice of the recording",
         f"window choice       : whole recording scanned every {p.scan_step:g} s "
         f"at -{REFERENCE_K:g} sigma",
         "thresholds          : " + ", ".join(f"-{k:g}" for k in SWEEP_K)
         + " sigma, each on the same windows",
         f"sharing             : median over {p.share_windows} common windows, "
         f"per threshold",
         "clusters            : re-clustered per channel on spikes pooled "
         "across its windows",
         "",
         "  threshold   channels     real    isolated   noise      median",
         "              with a unit  units   clusters   channels   missing"]
    for k in SWEEP_K:
        e = s["evidence"][_kkey(k)]
        L.append(f"  -{k:g} sigma    {e['channels_with_unit']:4d}/{s['n_channels']:<4d}"
                 f"  {e['real_units']:5d}    {e['isolated_clusters']:5d}      "
                 f"{e['noise_channels']:5d}      {fmt(e['median_missing_frac'], '6.1%')}")
    L.append("")
    if s["suggested_neg_k"] is None:
        L.append("suggested threshold : none - no channel has a unit at any "
                 "threshold")
    else:
        k = s["suggested_neg_k"]
        L.append(f"suggested threshold : -{k:g} sigma, picked by "
                 f"{s['votes'][_kkey(k)]} of {s['n_voting']} voting channels "
                 f"({s['confidence']:.0%})")
        L.append("votes               : " + ", ".join(
            f"-{kk} sigma: {v}" for kk, v in s["votes"].items()))
    L += ["",
          "rule: each channel picks the strictest threshold that keeps all its",
          "      real units (tagged isolated or multi-unit) and loses at most",
          f"      {MAX_MISSING:.0%} of its smallest unit's spikes below threshold"
          " (as",
          "      Bombcell's maxPercSpikesMissing); the area takes the most",
          "      common pick (ties go to the stricter). Channels with no unit",
          "      at any threshold do not vote.",
          "caveats: 'missing' is estimated from a truncated-Gaussian fit to the",
          "      smallest unit's trough depths, and is unknown below 20 spikes.",
          "      The tag's cut-offs were fitted at CHIRP's default threshold and",
          "      are untested at -4 sigma. CHIRP's sigma (MAD of the band-passed",
          "      raw trace) is not interchangeable with Kilosort's Th_*",
          "      thresholds, which apply after whitening."]
    for n in dr.notes:
        L.append(f"note: {n}")
    return "\n".join(L)


# ------------------------------------------------------------------ writers --
WINDOW_FIELDS = ["area", "neg_k"] + eng_video.STAT_FIELDS
CLUSTER_FIELDS = ["area", "channel", "neg_k", "cluster", "n_clusters",
                  "n_windows", "n_windows_present", "n_spikes",
                  "mean_amplitude_uV", "amplitude_sd_uV", "amplitude_cv_windows",
                  "firing_rate_sp_s", "rate_sd_sp_s", "snr", "half_width_ms",
                  "trough_to_peak_ms", "peak_uV", "sigma_uV", "wf_residual",
                  "share_count", "share_frac", "auto_quality", "tag_agreement",
                  "missing_frac"]
SWEEP_FIELDS = ["area", "channel", "neg_k", "n_windows", "n_clusters",
                "n_real", "n_isolated", "n_multiunit", "n_noise", "noise_flag",
                "total_rate_sp_s", "n_rejected", "sigma_uV", "share_count",
                "share_frac", "missing_frac", "residual_isolated", "loss_ok",
                "preferred_neg_k", "is_preferred"]


def write_rows_csv(rows, path, fields):
    """One CSV; floats rounded to 4 places, None written as an empty cell."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        wr = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        wr.writeheader()
        for r in rows:
            wr.writerow({k: (round(v, 4) if isinstance(v, float) else v)
                         for k, v in r.items()})
    return path.resolve()


def write_report(text, path):
    """Save the report next to the CSV, so a run leaves both behind."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    return path


def write_deep_outputs(results, out_dir, label=""):
    """
    The three deep tables and the suggestion file, for all groups of a run.
    Rows must already carry their 'area'. Returns {kind: path}.
    """
    out_dir = Path(out_dir)
    paths = dict(
        windows=write_rows_csv([r for dr in results for r in dr.window_rows],
                               out_dir / DEFAULT_CSV_NAME, WINDOW_FIELDS),
        clusters=write_rows_csv([r for dr in results for r in dr.clusters],
                                out_dir / CLUSTERS_CSV_NAME, CLUSTER_FIELDS),
        sweep=write_rows_csv([r for dr in results for r in dr.sweep],
                             out_dir / SWEEP_CSV_NAME, SWEEP_FIELDS))
    settings = asdict(results[0].params) if results else {}
    settings.pop("neg_k", None)
    payload = dict(
        chirp_version=__version__, label=label, mode="deep",
        detection="negative threshold at -neg_k x sigma (MAD of the band-passed "
                  "raw trace); suggested_neg_k uses the same convention as "
                  "opt.chirp.negK",
        thresholds=list(SWEEP_K), reference_neg_k=REFERENCE_K, settings=settings,
        areas={dr.name: dict(dr.suggestion, notes=dr.notes) for dr in results})
    paths["suggestion"] = out_dir / SUGGESTION_NAME
    paths["suggestion"].write_text(json.dumps(payload, indent=2), encoding="utf-8")
    return paths


# --------------------------------------------------------------------- main --
def main(argv=None) -> int:
    p = argparse.ArgumentParser(
        description="Survey Intan .dat channels and write cluster statistics. "
                    "Renders nothing; see dat_to_video.py for that.")
    p.add_argument("files", nargs="*", type=Path,
                   help="specific .dat files (default: --all)")
    p.add_argument("--all", action="store_true",
                   help="every amp-*.dat in the folder")
    p.add_argument("--folder", type=Path, default=None,
                   help="folder holding the amp-*.dat files (default: here)")
    p.add_argument("-o", "--out-dir", type=Path, default=None,
                   help="output directory (default: <folder>/stats)")
    p.add_argument("--mode", choices=sorted(PRESETS), default="fast",
                   help="fast: best windows at one threshold; deep: whole "
                        "recording, threshold sweep, suggested threshold")
    p.add_argument("-d", "--duration", type=float, default=10.0,
                   help="seconds of signal per window (default 10)")
    p.add_argument("-s", "--start", type=float, default=None,
                   help="fixed start in s, fast mode only (default: scan)")
    p.add_argument("-b", "--band", type=float, nargs=2, metavar=("LOW", "HIGH"),
                   default=list(DEFAULT_BAND), help="band-pass limits in Hz")
    p.add_argument("--fs", type=int, default=SAMPLE_RATE,
                   help="acquisition sample rate of the .dat files")
    p.add_argument("--neg-k", type=float, default=None,
                   help="detection threshold in sigma, fast mode only (default 5)")
    p.add_argument("--pos-k", type=float, default=8.0,
                   help="artifact rejection threshold in sigma (default 8)")
    p.add_argument("--artifact-k", type=float, default=18.0,
                   help="window is 'clean' if no sample exceeds this x sigma")
    p.add_argument("--max-clusters", type=int, default=3, choices=(1, 2, 3))
    p.add_argument("--segments", type=int, default=None,
                   help="windows measured per channel (default: 3 fast, 12 deep)")
    p.add_argument("--scan-step", type=float, default=None,
                   help="candidate window spacing in s (default: 60 fast, "
                        "the window length deep)")
    p.add_argument("--scan-range", type=float, nargs=2, default=None,
                   metavar=("FROM", "TO"), help="restrict the scan, in s")
    p.add_argument("--sampling", choices=("best", "stratified"), default=None,
                   help="window choice (default: best fast, stratified deep)")
    p.add_argument("--share-windows", type=int, default=None,
                   help="common windows for sharing (default: 1 fast, 5 deep)")
    p.add_argument("-w", "--workers", type=int, default=None,
                   help="processes (default: ~physical cores, serial if small)")
    p.add_argument("-q", "--quiet", action="store_true",
                   help="only print the final report")
    p.add_argument("--version", action="version", version=f"CHIRP {__version__}")
    args = p.parse_args(argv)

    folder = (args.folder or Path.cwd()).resolve()
    out_dir = args.out_dir or folder / "stats"
    if args.files:
        targets = [f if f.is_absolute() else folder / f for f in args.files]
        missing = [t for t in targets if not t.is_file()]
        if missing:
            print("Not found: " + ", ".join(map(str, missing)), file=sys.stderr)
            return 1
    elif args.all:
        targets = sorted(folder.glob("amp-*.dat"))
    else:
        print("Give one or more .dat files, or --all", file=sys.stderr)
        return 1
    if not targets:
        print(f"No amp-*.dat files found in {folder}", file=sys.stderr)
        return 1

    try:
        params, ignored = params_for_mode(
            args.mode, duration=args.duration, start=args.start,
            band=(args.band[0], args.band[1]), fs=args.fs, neg_k=args.neg_k,
            pos_k=args.pos_k, artifact_k=args.artifact_k,
            max_k=args.max_clusters, segments=args.segments,
            scan_step=args.scan_step,
            scan_range=tuple(args.scan_range) if args.scan_range else None,
            sampling=args.sampling, share_windows=args.share_windows)
    except ValueError as exc:
        print(f"bad parameters: {exc}", file=sys.stderr)
        return 1
    for name in ignored:
        print(f"note: {name} is ignored in deep mode", file=sys.stderr)
    workers = args.workers or default_workers(len(targets),
                                              planned_windows(targets, params))

    try:
        if args.mode == "deep":
            dr = run_deep([("all", targets)], params, workers=workers,
                          verbose=not args.quiet)[0]
            for r in dr.window_rows + dr.clusters + dr.sweep:
                r["area"] = "all"
            written = write_deep_outputs([dr], out_dir)
            report = (build_report(reference_result(dr), list(written.values()))
                      + "\n\n" + build_deep_report(dr))
        else:
            result = run_stats(targets, params, verbose=not args.quiet,
                               workers=workers)
            if not result.rows:
                print("no clusters found in any channel", file=sys.stderr)
                return 1
            for r in result.rows:
                r.update(area="all", neg_k=params.neg_k)
            written = dict(windows=write_rows_csv(
                result.rows, out_dir / DEFAULT_CSV_NAME, WINDOW_FIELDS))
            report = build_report(result, list(written.values()))
    except Cancelled:
        print("cancelled", file=sys.stderr)
        return 1

    report_path = write_report(report, out_dir / DEFAULT_REPORT_NAME)
    if not args.quiet:
        print()
    print(report)
    print(f"\n  -> {report_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
