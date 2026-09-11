"""Single place the version is written down, imported by every entry point.

Kept in its own module so nothing has to import a heavy one just to ask.

This is the NGL pipeline's fork of CHIRP. Upstream stays at 1.3.0; the local
suffix counts NGL-only changes, so reports and result files say which code
produced them. Do not re-vendor from upstream over this directory: it would
silently drop the pipeline's changes.

  +ngl.1  parallel per-channel pass (process pool)
  +ngl.2  deep mode: stratified windows, threshold sweep, pooled clusters,
          suggested threshold (lost-spikes rule); one self-updating progress
          line per phase; serial-work estimate from a short probe
"""

__version__ = "1.3.0+ngl.2"
