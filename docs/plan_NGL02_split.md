# Plan: split NGL02 into spike-only and LFP-only pipelines

**Status:** proposal, not yet executed. Awaiting user sign-off on the open decisions in §6 before any code changes.
**Author:** planning pass, 29.05.2026.

---

## 1. Motivation

`NGL02_postPhy.m` today does two things whose only shared dependency is the NGL01 output:

1. **Spike work** — read curated KS/Phy clusters, sort into trials, compute firing rate, run population dynamics. **Requires human Phy curation.**
2. **LFP work** — load FieldTrip data, do artifact rejection and time-frequency analysis. **Has nothing to do with Phy.**

Because LFP is gated behind a script whose name says "postPhy", users are forced to wait for curation before they can look at LFP. That's an artificial constraint — every input the LFP branch needs (`trialdef.mat`, `events.mat`, `*_FTcont.mat` / `*_stimOn2.mat`) is produced by NGL01 directly. The current single-script design conflates two timelines.

Goal: separate them into two scripts so the LFP path can run as soon as NGL01 finishes, in parallel with manual curation.

---

## 2. What's in NGL02 today (section-by-section)

| Section | Lines | Belongs to | Notes |
|---|---|---|---|
| 00. NGL00_Prep | 11-13 | shared | scaffolding |
| 00b. Areas recovery from preprocInfo | 15-40 | shared | scaffolding |
| `set_default`, `findSessions` | 42-47 | shared | scaffolding |
| 2.01 prepforsession | 52-54 | shared per-session | scaffolding |
| 2.01a checkNGL01Outputs | 56-60 | shared per-session | already gated by opt flags |
| 2.01b applyPreprocInfo | 62-98 | shared per-session | scaffolding |
| 2.02 offline video (`opt.offlineTrack`) | 100-105 | **spike** | produces `blob.mat`, only consumer is the social-tracking sub-step in 2.03 |
| 2.03 SPIKE branch (`opt.doSpikething`) | 107-276 | **spike** | spike.mat, neurons.mat, fireRate.mat, neuralDynamics.mat, social tracking |
| 2.04 LFP branch (`opt.doLFPthing`) | 278-349 | **LFP** | FT_data load, artifact rejection, spectrograms |
| cleanup clear | 351-356 | shared | per-iteration cleanup |

So the body splits cleanly: 2.02 + 2.03 go to the spike script; 2.04 goes to the LFP script. Everything else is shared scaffolding.

---

## 3. What moves where

### Spike script (renamed-or-kept `NGL02_postPhy.m`)

Body lives entirely in the existing 2.02-2.03 blocks. No code changes to the *content* — the work is to remove the LFP section and tidy comments.

Files produced (unchanged): `spike.mat`, `neurons.mat`, `fireRate.mat`, `neuralDynamics.mat`. Saved to `opt.spikeSorted` (spike.mat) and `opt.analysis` (everything else). Multi-area aware.

Functions consumed: `loadSpikes`, `sort2trials`, `calculate_fireRate_general`, `calculate_population_dynamics` + the population-dynamics family, plus social-tracking helpers (`processAndTrack_video`, `getSocialEvents`).

### LFP script (new `NGL02_LFP.m`)

Body comes from the current 2.04 block (lines 278-349). Cleanups recommended *as part of the move*:

- Drop the inline `data_all.mat` / `allconditions{x,y}` block (lines 312-315). It's a cross-session aggregation that doesn't belong inside a per-session loop. Slated for task #7 (NGL03_acrossSession); easiest to remove now and let task #7 reintroduce it in its proper home.
- Inline project-specific bits behind opt.proj_* flags as a passing courtesy. The `'ASL_Clean_Final_Correct'` testname (line 335), and the `opt.analysisCode = input.analysisCode` hack (line 338) are the obvious candidates. Strictly speaking that's task #8 territory, but they're small and isolated enough to handle inline.
- The cross-session save (`save(fullfile(input.analysis, 'allTFR_continuous.mat'), ...)` at line 331, similar at 344) belongs in NGL03_acrossSession. Same reason as `data_all.mat`. Strip during the move; task #7 reintroduces.

Files produced: `*_artifact_rejected.mat` (if `opt.artifdet`), `allTFR_continuous.mat` / `allTFR_chbych_trialparsed.mat` (if `opt.spectrogram` — TO BE MOVED to NGL03_acrossSession). Multi-area awareness for LFP is not in scope (LFP-side functions currently take all-channel FT_data; if needed later it's a separate task).

Functions consumed: `artifact_detRej_lfp`, `continous_MTspectrogram`, `trialparsed_MTspectrogram`, optionally `MAT2FieldTrip`, optionally `vFLIP_NGL` (commented today).

### Shared infrastructure

Two viable handlings:

**Option A — duplicate.** Each script carries its own copy of:
- Sections 00, 00b, set_default, findSessions before the loop.
- Sections 2.01, 2.01a, 2.01b inside the loop.

~85 lines of duplication. Each script reads top-to-bottom without jumps.

**Option B — extract.** New helper `prepSession(input, opt)` (functions/pipeline/) wraps the per-session preamble (`prepforsession` + `checkNGL01Outputs` + `applyPreprocInfo`). Pre-loop sections still live inline in each script. Per-session preamble shrinks to one call.

Less duplication; one extra hop when reading.

My recommendation: **B**, because (i) the per-session preamble is already 47 lines, almost all of it the `applyPreprocInfo` try/catch tower; (ii) when we update one of those steps later it's good to have one place to edit; (iii) the pre-loop block stays inline because it interacts with `input`/`opt` lexically and is short enough not to obscure flow.

---

## 4. SetAndRunMe wiring

### Current sections (NGL_SetAndRunMe.m template)

```
3.1 NGL01_Main
3.2 NGL02_postPhy
3.3 NGL03_plotting   (TODO)
3.4 NGL04_aggregate  (SocialLearning only)
```

### Proposed sections

```
3.1 NGL01_Main             — preprocessing; produces NGL01 outputs
3.2 NGL02_LFP              — LFP analysis; can run any time after 3.1
3.3 NGL02_postPhy          — spike analysis; run after Phy curation
3.4 NGL03_plotting         — (TODO)
3.5 NGL03_acrossSession    — cross-session aggregation (planned, task #7)
```

LFP goes before postPhy in the user-visible ordering because:
- Chronologically natural (no Phy wait).
- User can review LFP results while curating in Phy.
- Users running spike-only work just skip 3.2; users running LFP-only just skip 3.3.

Both 3.2 and 3.3 are independent — running one doesn't require the other to have completed.

---

## 5. `opt.doSpikething` / `opt.doLFPthing` after split

Two possible roles:

**(a) Master gate per script.** Set the corresponding flag true at the top of each script; the script's user-visible existence is the actual gate.

**(b) Per-session gate inside each script.** Keep them as inner guards so users can skip certain sessions while looping.

(b) is more flexible and is what the current code does. Recommend keeping (b) — minor cost (one `if opt.doX` per script) for meaningful flexibility.

`opt.useTrack`, `opt.offlineTrack`, `opt.popDyn.do`, `opt.neurDyn.do`, `opt.artifdet`, `opt.spectrogram` all stay where they are; they're sub-step flags inside the bigger spike or LFP guard.

---

## 6. Open decisions (need your sign-off before implementation)

1. **Naming.** Confirm `NGL02_postPhy.m` keeps its name (informative: "after Phy") and the new file is `NGL02_LFP.m`. Alternatives:
   - `NGL02_LFP.m` + rename existing to `NGL02_spike.m` (clearer symmetry, breaks SetAndRunMe references everywhere).
   - `NGL02b_LFP.m` (kind-of-hierarchical).
   - Keep `NGL02_postPhy.m` + add `NGL02_LFP.m` (**my recommendation**).

2. **SetAndRunMe ordering.** LFP before postPhy (chronological)? Or LFP after postPhy (current alphabetical-ish)? My recommendation: **LFP first** (3.2), postPhy second (3.3).

3. **Scaffolding factoring.** Duplicate (Option A) or extract `prepSession` helper (Option B)? **My recommendation: B.**

4. **Section 2.02 video tracker.** Confirm it stays with the spike pipeline (its consumer `opt.useTrack` lives in the spike branch).

5. **`data_all.mat` / `allconditions` block.** Strip during the split (it'll be re-added properly by task #7), or leave in the LFP script for now and let task #7 remove it later? My recommendation: **strip now** — it's broken-ish anyway and the LFP script is cleaner without it. Task #7 reintroduces in NGL03_acrossSession.

6. **Cross-session save inside LFP block.** Strip (the `allTFR_*{x,y}` accumulation and the end-of-loop save) and let NGL03_acrossSession own it, or keep for now? My recommendation: **strip now**. Same reasoning as 5.

7. **Project-specific code in LFP block** (`'ASL_Clean_Final_Correct'` testname, `opt.analysisCode = input.analysisCode` hack). Gate behind a new `opt.proj_ASL` flag during the move (slight overlap with task #8), or leave commented-with-warning and let task #8 own it? My recommendation: **leave for task #8**; tag in code as a `% TODO proj_ASL gate (task #8)` comment.

8. **Section numbering inside the new script.** Restart at `01`, `02`, … in `NGL02_LFP.m`, or keep some shared scheme (`2.04` would be odd after extraction)? My recommendation: **restart at `01`** — each script has its own internal numbering.

---

## 7. Risks and edge cases

- **Cached files.** Splitting doesn't change output filenames or locations, so no cache invalidation.
- **NGL_SetAndRunMe templates in user projects.** Each project's own `NGL_SetAndRunMe.m` (e.g. colorWheel's, S3_ExperimentArena's) will continue to work — they call `NGL02_postPhy` which still exists. To run LFP separately, users add the new section 3.2 call. Backwards compatible.
- **`opt.doLFPthing` set true in old SetAndRunMe.** If a user's existing SetAndRunMe has `opt.doLFPthing=true` but still calls only `NGL02_postPhy`, the LFP branch will silently no-op after the split (because the LFP block has been removed from postPhy). Need to add a one-line warning when `opt.doLFPthing=true` reaches the end of NGL02_postPhy without LFP being run — point user at the new `NGL02_LFP.m` script.
- **In-memory leakage between scripts.** If a user runs both 3.2 and 3.3 in the same MATLAB session, the existing `clear` at the end of each per-session iteration handles the per-session cleanup. Between scripts, the workspace cleanly contains the last script's outputs only.
- **Path setup.** `set_default` adds `functions/pipeline/` and `functions/analysis/` etc. to the path. Both scripts call it, so no order dependency. Good.
- **`checkNGL01Outputs` already script-aware** via `opt.doSpikething` / `opt.doLFPthing`. After split, each script sets exactly the right gate, and the checker reports only the relevant missing files. No changes needed to the checker.

---

## 8. Concrete implementation plan (once decisions are locked)

1. Create `functions/pipeline/prepSession.m` (Option B for scaffolding) wrapping per-session prep.
2. Create `NGL02_LFP.m`. Body = current section 2.04, modulo stripping per §6 decisions 5, 6, 7. Pre-session preamble = `prepSession` call. Pre-loop preamble = NGL00_Prep + Areas recovery + set_default + findSessions.
3. Edit `NGL02_postPhy.m`. Remove section 2.04 entirely. Add a tail-of-iteration warning when `opt.doLFPthing=true` (§7 mitigation). Replace per-session preamble with `prepSession` call.
4. Update `NGL_SetAndRunMe.m` (toolbox template) with new section 3.2 → `NGL02_LFP`. Bump postPhy to 3.3. Bump plotting to 3.4. Add comment block explaining the chronology.
5. Update `README.md` Pipeline Stages table and section comments.
6. Update existing project SetAndRunMe files (colorWheel, S3_ExperimentArena) as a one-time follow-up — or leave it to the user.

Each step is small and can be code-reviewed independently.

---

## 9. Estimated effort

| Step | Effort |
|---|---|
| `prepSession` helper | trivial (15 min) |
| Create NGL02_LFP from current section 2.04 | small (30 min, plus deciding which clean-ups to apply during the move) |
| Trim NGL02_postPhy | trivial (15 min) |
| Update SetAndRunMe template + README | small (15 min) |
| Update existing project SetAndRunMe files | trivial (5 min each, two projects) |
| Smoke test on colorWheel | depends on you |

Total toolbox work: **~1.5 hours**, low risk because no analysis-side code changes are required; this is structural.
