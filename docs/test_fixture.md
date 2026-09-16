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
| `../make_fake_intan.py` | separate: raw INTAN `.dat` generator (pre-sorting tests) |

## Use it

```matlab
addpath('tests/fixture')
makeFakeStudy('D:\TESTSTUDY')            % three sessions, ~8 min each
checkFixture('D:\TESTSTUDY')             % errors if any check fails
```

Then point a study at it: `subjects = {'MRX'}`, `dates = {'19850214',
'19860704', '19870321'}` in `NGL_SetAndRunMe`, with the study root set to the
folder above. `makeFakeStudy(root, 'length', 'full')` gives ~20-minute
sessions; `'sessions', {'19850214'}` builds just one.

The dates are **pre-1990**, in the normal `YYYYMMDD` format, so a test session
can never be mistaken for a real recording. The subject `MRX` exists for the
same reason.

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
data\trialSorted\MRX\<sess>\    events.mat  trialdef.mat  condition.mat
                                <sess>_itiOn.mat  <sess>_stimOn1.mat  <sess>_stimOn2.mat
data\spikeSorted\MRX\<sess>\    spike.mat
data\analysis\MRX\<sess>\       (empty; analyses write here)
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

- units driven 40–70 ms after `stimOn1`, one driven by reward, one locked to the
  tagging stream, one whose amplitude drifts ~45 % across the session, one
  cluster left labelled `noise` as Phy curation would leave it

**IMU** (3-axis accelerometer, m/s² at 1 kHz — the shape `GetMotionSensors`
writes for INTAN)

- a sharp 30 ms jerk at every peck, ~50× the noise MAD on the jerk magnitude
- slow sway while walking, ramped in and out (a step onset would itself look
  like a peck to a jerk detector)

**Events** — the digital-input record `NGL01` would extract: `preIni`, `itiOn`,
`stimOn1`, `stimOn2`, `bhv`, `rwd`/`pun`, and an `end1`/`end2`/`end3` code
closing **every** trial, which is the pairing `trialdefGen` needs to find trials
at all. Tagging stimuli carry project codes 8001–8006 after their `stimOn1`.

## The checks

`checkFixture` runs the real functions — not copies — and fails when an analysis
**misses** a planted effect *or* **finds one where nothing was planted**. The
second half is the one a test usually forgets, and it is why each check has a
negative control.

| Check | What it asserts |
|---|---|
| `layout` | the files a session must have, where the pipeline looks for them |
| `trials` | counts agree, outcomes partition the trials, latencies inside their windows, lengths vary |
| `events` | times monotonic, every trial paired `itiOn`→end code, tagging codes present |
| `spikes` | driven units beat their own baseline and undriven ones do not; drift and curation labels survive |
| `imu` | every peck recoverable from the jerk magnitude, with nothing else above threshold |
| `csd` | the strongest sink lands within a contact of 325 µm, ~120 ms after `stimOn1` |
| `contrast` | a cluster at ~20 Hz in NCL (p ≈ 0.002) and **no cluster in STR** |
| `nft` | driven channels reach every harmonic, undriven ones at most one |

Typical run (session `19850214`):

```
  [PASS] trials    298 trials, 42 correct / 13 incorrect / 7 omitted / 2 aborted / 234 passive
  [PASS] spikes    8 units; stimOn1 drive 5.6x (others 1.2x), drift -25% over the session
  [PASS] imu       recall 100% of 80 pecks, precision 100%; peck jerk is 53x the noise MAD
  [PASS] csd       sink at 300 um (planted 325, spacing 50) at +121 ms
  [PASS] contrast  NCL cluster p=0.002 at 23 Hz (planted 20, band 15-30); STR no cluster
  [PASS] nft       1.3 Hz: NCL z 4.9, 8.0 harmonics; STR z 0.1, 0.0
```

## What building it exposed

Three things the fixture surfaced that are worth knowing before real mixed-task
sessions arrive:

1. **A trial without an end code is not a trial.** `trialdefGen` pairs `itiOn`
   with `end1`/`end2`/`end3`; an outcome marker (`rwd`, `pun`) is not enough.
2. **An alignment can be missing from a trial.** The tagging stream has no
   `stimOn2`, so those rows carry a NaN offset and the trial parser drops them.
   The surviving trials then no longer line up with the per-trial `condition`
   vectors — the fixture records `FT_data.trialinfo` (original trial indices) so
   a consumer can subset. **NGL07's contrast and CSD paths do not do this yet**;
   they assume every trial has the alignment, which holds only when one session
   is one task.
3. **CSD needs a common time axis.** Real trials end when the animal responds,
   so they differ in length; `computeCSD` requires them squared off first. The
   check cuts a fixed window around the alignment before calling it, and
   selects trials that span the whole window — not the most common length,
   because the most numerous trials here are the shortest ones.

Separately, `detect_jerkEvents` could not be used for the IMU check: its Pass 2
learns templates that keep a DC component, so on a strictly positive signal
(a jerk *magnitude*) the matched-filter output sits above threshold everywhere
and `find_peaks` returns one event for the whole recording; on a mean-centred or
signed trace it saturates the refractory limit instead, and in both cases the
result does not respond to `ThresholdSD2` at all. The IMU check therefore uses a
plain threshold on the jerk magnitude, which is what the fixture actually
guarantees. The detector itself is worth a separate look.

Last modified 16.09.2026 (Jesus) — new (#F test fixture).
