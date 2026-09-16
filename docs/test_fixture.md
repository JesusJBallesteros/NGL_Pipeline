# Test fixture — subject MRX

A synthetic study whose answers are known in advance. Every analysis in the
pipeline can be checked against **what was planted** rather than against what
looks plausible, and everyone checks against the same thing: the generator is
seeded, so the same call produces the same bytes on any machine.

The repository carries the **generator, not the data**. One session is ~270 MB,
which does not belong in git, and a fixed seed gives the same bytes anyway.

Files: `tests/fixture/`

| File | What it does |
|---|---|
| `makeFakeStudy.m` | builds the three shipped sessions and a README beside them |
| `makeFakeSession.m` | one session: blocks, trials, events, and the files the pipeline expects |
| `makeFixtureSignals.m` | the LFP, spike trains and IMU, and the record of what was planted |
| `makeFixtureFT.m` | packages LFP as FieldTrip data, continuous or per alignment |
| `checkFixture.m` | runs the real analyses and compares them to the planted answers |
| `makeFixtureExamples.m` | regenerates `examples/` — the committed run log and figures |
| `NGL_RunFixture.m` | the user-level script: runs the pipeline over MRX as its own study |
| `fixtureNFTblocks.m` | tagging windows valid in every session, for `opt.nft.blocks` |
| `examples/` | what a correct run looks like, kept in the repo |
| `../make_fake_intan.py` | separate: raw INTAN `.dat` generator (pre-sorting tests) |

## Use it

```matlab
addpath('tests/fixture')
makeFakeStudy('D:\TESTSTUDY')            % three sessions, ~10 min each
checkFixture('D:\TESTSTUDY')             % errors if any check fails
```

`makeFakeStudy(root, 'length', 'full')` gives ~20-minute sessions;
`'sessions', {'19850214'}` builds just one.

The dates are **pre-1990**, in the normal `YYYYMMDD` format, so a test session
can never be mistaken for a real recording. The subject `MRX` exists for the
same reason.

### At the user-script level

`tests/fixture/NGL_RunFixture.m` is the `NGL_SetAndRunMe` for the fixture: same
sections, same `opt`, same stages. Copy it beside your own script, point
`datadrive` / `studyname` at a location of your choice, and run it whenever you
want to know whether the pipeline still returns the planted answers. It runs
`NGL02_LFP`, `NGL07_LFPanalysis` and `NGL08_NFT` (the earlier stages have
nothing to do — the generator already wrote what they produce), then calls
`checkFixture`.

One wrinkle it has to work around, which a real study shares:
`opt.nft.blocks` is **one set of windows for the whole run**, while each
session's tagging block starts at a slightly different second. The script sets
it from `fixtureNFTblocks`, which returns the window common to every session
requested — narrower than any single session's block, and correct in all of
them. A study whose blocks drift more than this one's would need per-session
windows, which `NGL08_NFT` does not currently accept.

**It is a separate study, not a subject added to yours.** Study-level outputs
are named by analysis rather than by subject, so a fixture run inside your
study would write over your aggregates — and synthetic results have no business
in the same tree as real ones. Run it in parallel with your own work, in its
own `studyname`.

### Why it cannot contaminate a real analysis

The separation above is the actual protection. Underneath it there is a
backstop, because "keep them apart" is a rule that people and scripts break:
the generator writes a `SYNTHETIC_DATA.txt` marker into the subject's folder in
every data tree, and `set_default` — the single point where **every** stage
resolves its subject list, `NGL03_aggregate` and `NGL04_fireRate` included —
reads it through `isSyntheticSubject` and applies three rules
(`guardSyntheticSubjects`):

| You ask for | What happens |
|---|---|
| `subjects = 'all'` | synthetic subjects are skipped, with one line saying so |
| `subjects = {'MRX'}`, in its own study | runs, printing a `*** SYNTHETIC DATA RUN ***` banner |
| `{'MRX','A42'}`, or MRX inside a study holding real subjects | **error** — refused outright |

The `isolation` check tests all of this rather than trusting it: it asserts the
marker is on disk in every tree, that `'all'` drops the subject, that a
fixture-only run is flagged, and that both mixing cases raise their specific
error. Deleting the marker removes the backstop — the file says so.

## What a session contains

Three task blocks in order, each with its own trial anatomy — this is one
session containing several tasks, which is itself something the pipeline has to
survive:

| Block | Trials | Anatomy |
|---|---|---|
| `dms` | 40 | perched, pecking a screen: ITI → sample → delay → test → response |
| `arena` | 24 | freely moving, 3 of 6 screens alternating: ITI → cue → travel → response → return |
| `nft` | 234 | passive rapid stream, 1.3 Hz then 2.6 Hz, no response |

Outcomes are `correct` / `incorrect` / `omission` / `aborted` (`passive` for the
tagging stream), each with its own timing consequence: an omission runs the
response window out, an abort ends early, a response comes at a latency drawn
per trial between a floor (nobody reacts instantly) and the window's ceiling.
**Trial lengths therefore vary**, from 0.4 s to over 10 s, and code that assumes
a fixed trial length fails here — which is the point.

Files written, in the layout `set_default` builds:

```
data\preprocessing\MRX\<sess>\  <sess>_FTcont.mat  EventRecord.mat
                                MotionData_raw.mat  fixture_truth.mat
data\trialSorted\MRX\<sess>\    events.mat  trialdef.mat  condition.mat  blocks.mat
                                <sess>_itiOn.mat  <sess>_stimOn1.mat  <sess>_stimOn2.mat
data\spikeSorted\MRX\<sess>\    spike.mat
data\analysis\MRX\<sess>\       (empty; analyses write here)
data\*\MRX\                     SYNTHETIC_DATA.txt   (the isolation marker)
```

`fixture_truth.mat` holds the planted answers in machine-readable form, so a
test asserts against the design and not against a previous run.

## What is planted

**LFP** (32 channels, 2 shanks, the lab's ATLAS map; NCL on shank 1, STR on shank 2)

- a 20 Hz burst 0.35 s after `stimOn2` on **correct** dms trials, **NCL only**
- 6 Hz theta while walking in arena trials, both areas
- a stimulus-locked response at the tagging rates, **NCL only**
- a current **sink at 325 µm on shank 1**, ~60 ms after `stimOn1`, built by
  double-integrating a known CSD so the analysis recovers what was planted
- trouble on purpose: one dead channel, one noisy channel, movement artifacts

**Spikes** (8 units, in the shape `NGL02_postPhy` writes)

- units driven 40–70 ms after `stimOn1` — in **every** task, including the
  passive stream, because a stimulus drives a stimulus-responsive unit whether
  or not there is anything to do about it
- one unit driven by reward, one locked to the tagging stream, one whose
  amplitude drifts ~45 % across the session, one cluster left labelled `noise`
  as Phy curation would leave it

**IMU** (3-axis accelerometer, m/s² at 1 kHz — the shape `GetMotionSensors`
writes for INTAN)

- a sharp 30 ms jerk at every peck, ~50× the noise MAD on the jerk magnitude
- slow sway while walking, ramped in and out (a step onset would itself look
  like a peck to a jerk detector)

**Events** — the digital-input record `NGL01` would extract:

- `preIni`, `itiOn`, `stimOn1`, `stimOn2`, then the response `bhv`;
- `rwd` (correct) or `pun` (incorrect) **120 ms after the peck** — the outcome
  is the response's consequence, so it follows it and cannot float free of it;
- an `end1`/`end2`/`end3` code closing **every** trial, at a **fixed maximum
  length**: the same distance from the last stimulus whatever the animal did.
  A fast response and a slow one therefore give the same trial length and the
  same analysis window, which is what lets the LFP window reach ±1 s beyond the
  effects (`opt.toi = [-2 1.8]`) on every trial rather than on some of them;
- block markers (`16`/`17` plus a task code) around each block — see below;
- tagging stimuli carry project codes 8001–8006 after their `stimOn1`.

**Blocks** — a block marker is read **once per session, not per trial**
(`sessionBlocks`), and gives a table of stage boundaries: label, time range,
and the trials inside it. That is what an extinction-style paradigm needs to
draw its stage boundaries, and here it does a second job: not every task
carries every analysis, so each analysis asks the block table which trials it
may use. A tagging spectrum over a delay-match block is not a weaker result,
it is a meaningless one. `blocks.mat` holds the recovered table; the gating is
checked, including the case of an analysis nobody declared, which gets no
trials at all.

## The checks

`checkFixture` runs the real functions — not copies — and fails when an analysis
**misses** a planted effect *or* **finds one where nothing was planted**. The
second half is the one a test usually forgets, and it is why each check has a
negative control.

| Check | What it asserts |
|---|---|
| `layout` | the files a session must have, where the pipeline looks for them |
| `isolation` | the marker is in every tree, and the guard skips / flags / refuses as it should |
| `trials` | counts agree, outcomes partition the trials, latencies inside their windows, the end event fixed, lengths vary |
| `events` | times monotonic, every trial paired `itiOn`→end code, outcomes follow the peck, tagging codes present |
| `blocks` | the recovered blocks match the ones that ran, they tile the session, and the analysis gating holds |
| `channels` | the dead channel is dead (and answers no tagging), the noisy one is noisy |
| `spikes` | driven units beat their own baseline and undriven ones do not; drift and curation labels survive |
| `psth` | the same through the pipeline's own path (`sort2trials` → `calcFireRate`), per task, against the pooled undriven units |
| `imu` | every peck recoverable from the jerk magnitude, with nothing else above threshold |
| `csd` | the strongest sink lands within a contact of 325 µm, ~120 ms after `stimOn1` |
| `contrast` | a cluster at ~20 Hz in NCL (p ≤ 0.01) and **nothing comparable in STR** |
| `nft` | driven channels reach every harmonic, undriven ones at most one |

The full run over all three sessions is committed under
`tests/fixture/examples/expected_output.txt`, alongside one figure per
analysis.

## The examples

`examples/` is what a correct run looks like, so nobody has to run anything to
find out. It is produced by `makeFixtureExamples` and never edited by hand:

| File | What it shows |
|---|---|
| `expected_output.txt` | the whole run: generation, then every check on every session |
| `example_contrast.png` | correct vs incorrect in NCL: the planted 20 Hz burst and its cluster |
| `example_csd.png` | the sink at ~325 µm on shank 1, with the LFP traces beside it |
| `example_tagging.png` | the 1.3 Hz stream: harmonics, SNR and summed response |
| `example_psth.png` | one PSTH per task, the driven unit against an undriven one |
| `example_pecks.png` | peck detection in both pecking tasks: one interaction, then all of them averaged |

The figures are made by running the real functions, so a change that breaks an
API breaks this script — which is the point of keeping it runnable rather than
pasting screenshots. Regenerate with:

```matlab
makeFixtureExamples                            % all three sessions, ~15 min
makeFixtureExamples('sessions', {'19850214'})  % just the first
```

Two things worth knowing when reading them. In `example_psth.png` the tagging
panel shows a second peak about 385 ms before the alignment — that is the
previous stimulus of the 2.6 Hz stream, not an artifact — and the rate falls to
zero past +220 ms because the trial itself has ended there. In
`example_pecks.png` the two tasks look alike because the peck is the same beak
movement in both; the walking that distinguishes them is slow, and the 20 Hz
high-pass the detector runs on removes it by design.

## What building it exposed

Bugs the fixture found in the pipeline, **now fixed**:

1. **`sort2trials` matched alignments by position, not by name.** It walked
   `opt.alignto` and read `trialdef{2,a}` at the same index, so asking for a
   subset (`opt.alignto = {'stimOn1'}`) returned the *first* alignment's
   windows under the name `stimOn1` — silently, and looking entirely plausible:
   a PSTH with no response in it. It now looks the column up by name and keeps
   position only as a fallback.
2. **Every signed map crashed on schema defaults.** `opt.lfp.plot.divergingColormap`
   defaults to `''`, documented as "use the built-in blue-white-red", but
   `lfpStyle` passed the empty string to `colormap()` — `Invalid function name ''`
   for any contrast or CSD figure made from unmodified defaults.
3. **Long TFR cache keys silently failed to cache.** The key encodes every
   setting, which under a deep study folder pushed the path past the Windows
   260-character limit; MATLAB reports that as *"the file appears to be
   corrupt"*, so it reads as a broken cache rather than a long name. Keys over
   120 characters now keep their readable head and fold the rest into a hash.

Things to know before real mixed-task sessions arrive:

4. **A trial without an end code is not a trial.** `trialdefGen` pairs `itiOn`
   with `end1`/`end2`/`end3`; an outcome marker (`rwd`, `pun`) is not enough.
5. **An alignment can be missing from a trial.** The tagging stream has no
   `stimOn2`, so those rows carry a NaN offset and the trial parser drops them.
   The surviving trials then no longer line up with the per-trial `condition`
   vectors — the fixture records `FT_data.trialinfo` (original trial indices) so
   a consumer can subset. **NGL07's contrast and CSD paths do not do this yet**;
   they assume every trial has the alignment, which holds only when one session
   is one task.
6. **CSD needs a common time axis.** Real trials end when the animal responds,
   so they differ in length; `computeCSD` requires them squared off first. The
   check cuts a fixed window around the alignment before calling it, and
   selects trials that span the whole window — not the most common length,
   because the most numerous trials here are the shortest ones.
7. **A dead contact invents a sink.** Zeroing a channel punches a hole in the
   depth profile, and a second spatial derivative across that hole produces a
   sink exactly where the hole is. The fixture keeps its bad channels on the
   shank the CSD does not use and says so; interpolating a dead contact before
   a CSD is a real requirement the pipeline does not meet yet.

Still open: **`detect_jerkEvents` cannot be used on a jerk magnitude.** Its
Pass 2 learns templates that keep a DC component, so on a strictly positive
signal the matched-filter output sits above threshold everywhere and
`find_peaks` returns one event for the whole recording; on a mean-centred or
signed trace it saturates the refractory limit instead. In both cases the
result does not respond to `ThresholdSD2` at all — 2 or 24 gives the same
count. The IMU check therefore uses a plain threshold on the jerk magnitude,
which is what the fixture actually guarantees. The likely one-line cause is
that learned templates are normalised but never mean-removed.

## Two fixture bugs worth remembering

Both were caught by looking at the figures rather than at the numbers, which is
the argument for keeping the examples in the repo:

- The **dead channel was zeroed before** the planted responses were added, so it
  still carried them — and because its noise floor was the lowest on the probe,
  every SNR measure ranked it first and the tagging figure picked it to show
  (z ≈ 2000 on a channel documented as dead). The `channels` check now asserts
  that the dead channel is dead *and* answers no tagging.
- The **arena's walking sway started as a step**, and a step in acceleration is
  exactly what a jerk detector should find: 23 spurious "pecks" per session. It
  is ramped in and out now.

Last modified 16.09.2026 (Jesus) — new (#F test fixture).
