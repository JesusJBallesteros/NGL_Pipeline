"""Synthesise a 'one file per channel' Intan recording for testing CHIRP.

Writes headerless int16 little-endian at 0.195 uV/bit, the format CHIRP reads.

Channels cycle through five archetypes (injected trough amplitudes, white
noise of 6 uV before band-pass, which leaves ~4.3 uV after 450-8000 Hz):

  0  two units, 90 and 45 uV
  1  two units, 110 and 50 uV
  2  one unit, 75 uV
  3  noise only
  4  one near-threshold unit, 30 uV (~5 sigma after band-pass): detected at
     -4 and -5 sigma, largely lost at -6 - a known answer for a threshold sweep

A common-mode transient is injected into every channel so the cross-channel
sharing pass has something to find. The defaults (4 channels, 60 s, seed 7)
reproduce the original fixture.

    python tests/make_fake_intan.py <out_dir> [--duration S] [--channels N] [--seed N]
"""
import argparse
from pathlib import Path

import numpy as np

UV_PER_BIT = 0.195
FS = 30000
ARCHETYPES = [
    [(12, 90.0), (25, 45.0)],
    [(8, 110.0), (30, 50.0)],
    [(15, 75.0)],
    [],
    [(20, 30.0)],
]


def spike_waveform(fs, amp_uv, width_ms=0.4):
    """A trough followed by a positive rebound, the shape the detector expects."""
    n = int(round(width_ms * 3 * fs / 1000))
    t = np.arange(n) / fs * 1000.0
    trough = -amp_uv * np.exp(-0.5 * ((t - width_ms) / (width_ms / 2.5)) ** 2)
    rebound = 0.35 * amp_uv * np.exp(-0.5 * ((t - 2.1 * width_ms) / (width_ms / 1.6)) ** 2)
    return trough + rebound


def build_channel(rng, n, duration, rates_amps, noise_uv=6.0):
    x = rng.normal(0.0, noise_uv, n)
    for rate, amp in rates_amps:
        idx = rng.integers(200, n - 200, int(rate * duration))
        w = spike_waveform(FS, amp)
        for i in idx:
            x[i:i + len(w)] += w
    return x


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("out_dir", nargs="?", default="fake_recording")
    ap.add_argument("--duration", type=float, default=60.0, help="seconds")
    ap.add_argument("--channels", type=int, default=4)
    ap.add_argument("--seed", type=int, default=7)
    args = ap.parse_args(argv)

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)
    n = int(FS * args.duration)

    # Common-mode transient shared by every channel, so session_sharing has a
    # signal that is present everywhere at the same instants.
    common = np.zeros(n)
    cw = spike_waveform(FS, 55.0)
    for i in rng.integers(200, n - 200, int(3 * args.duration)):
        common[i:i + len(cw)] += cw

    # Written one channel at a time: a long recording held whole would need
    # several hundred MB per channel in float64.
    for c in range(args.channels):
        x = build_channel(rng, n, args.duration, ARCHETYPES[c % len(ARCHETYPES)])
        raw = np.clip(np.round((x + common) / UV_PER_BIT), -32768, 32767).astype("<i2")
        (out / f"amp-A-{c:03d}.dat").write_bytes(raw.tobytes())
    print(f"wrote {args.channels} channels, {args.duration:g} s at {FS} Hz -> {out}")


if __name__ == "__main__":
    main()
