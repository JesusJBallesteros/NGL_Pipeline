"""Synthesise a small 'one file per channel' Intan recording for testing.

Writes headerless int16 little-endian at 0.195 uV/bit, the format CHIRP reads.
Four channels: two carry two amplitude classes of spikes, one carries a single
class, one is noise only. A common-mode transient is injected into every
channel so the cross-channel sharing pass has something to find.
"""
import sys
from pathlib import Path
import numpy as np

UV_PER_BIT = 0.195
FS = 30000
DUR = 60.0
N = int(FS * DUR)
rng = np.random.default_rng(7)


def spike_waveform(fs, amp_uv, width_ms=0.4):
    """A trough followed by a positive rebound, the shape the detector expects."""
    n = int(round(width_ms * 3 * fs / 1000))
    t = np.arange(n) / fs * 1000.0
    trough = -amp_uv * np.exp(-0.5 * ((t - width_ms) / (width_ms / 2.5)) ** 2)
    rebound = 0.35 * amp_uv * np.exp(-0.5 * ((t - 2.1 * width_ms) / (width_ms / 1.6)) ** 2)
    return trough + rebound


def build_channel(rates_amps, noise_uv=6.0):
    x = rng.normal(0.0, noise_uv, N)
    for rate, amp in rates_amps:
        n_spk = int(rate * DUR)
        idx = rng.integers(200, N - 200, n_spk)
        w = spike_waveform(FS, amp)
        for i in idx:
            x[i:i + len(w)] += w
    return x


def main(out_dir):
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)

    # Common-mode transient shared by every channel, so session_sharing has a
    # signal that is present everywhere at the same instants.
    common_idx = rng.integers(200, N - 200, int(3 * DUR))
    common = np.zeros(N)
    cw = spike_waveform(FS, 55.0)
    for i in common_idx:
        common[i:i + len(cw)] += cw

    channels = {
        "amp-A-000": build_channel([(12, 90.0), (25, 45.0)]),
        "amp-A-001": build_channel([(8, 110.0), (30, 50.0)]),
        "amp-A-002": build_channel([(15, 75.0)]),
        "amp-A-003": build_channel([]),
    }
    for name, x in channels.items():
        y = x + common
        raw = np.clip(np.round(y / UV_PER_BIT), -32768, 32767).astype("<i2")
        (out / f"{name}.dat").write_bytes(raw.tobytes())
        print(f"{name}.dat  {raw.nbytes / 1e6:.1f} MB")
    print(f"wrote {len(channels)} channels, {DUR:g} s at {FS} Hz -> {out}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "fake_recording")
