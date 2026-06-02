# NGL02_postPhy — Audit (pre-refactor)

**Scope:** read-through of `NGL02_postPhy.m` and every function it calls, against the conventions established during the NGL01 refactor.
**Priority lens:** multi-area / structural consistency first.
**Author:** audit pass, 27.05.2026.

This is a findings document, not a redesign. Each item is tagged with a severity (🔴 bug, 🟠 inconsistency, 🟡 cleanup, 🔵 design question) so we can prioritise edits in a follow-up pass.

---

## 1. Genuine bugs in `NGL02_postPhy.m`

### 1.1 🔴 `prepforsession` call signature is wrong
NGL02 currently does:
```matlab
[input.sessions(input.run(1)).info, opt] = prepforsession(input, opt);
```
But `prepforsession` returns `[input, opt]` — the **full** input struct. So this assigns the entire input struct into `input.sessions(x).info`, silently corrupting that field. NGL01 does it correctly:
```matlab
[input, opt] = prepforsession(input, opt);
```
The reason nothing visibly breaks today is that the rest of NGL02's loop never reads `input.sessions(x).info` after the call. But it's wrong and will bite us when multi-area or any new prepforsession side effect needs to flow out.

### 1.2 🔴 `conditions` vs `condition` typo when calling fireRate
NGL02 loads `condition` (singular, from `condition.mat`) but passes `conditions` (plural, undefined) to `calculate_fireRate_general`:
```matlab
if ~exist('condition','var'), load(fullfile(opt.trialSorted, "condition.mat")); end
% ...
fireRate = calculate_fireRate_general(neurons, events, conditions, opt, param);
```
`calculate_fireRate_general`'s signature is `(neurons, events, condition, opt, param)` — singular. This will throw "Undefined function or variable 'conditions'" on the first session.

### 1.3 🔴 `artifact_detRej_lfp` call signature mismatch
NGL02 calls:
```matlab
FT_data = artifact_detRej_lfp(FT_data, condition, opt, param);
```
but the function is `artifact_detRej_lfp(FT_data, opt)` (only two args). The extra args get silently dropped because of MATLAB's positional binding — `condition` becomes `opt`, which then has no `artZvalue` field, and the actual `opt` is ignored. Same kind of bug, masked.

### 1.4 🔴 `continous_MTspectrogram` / `trialparsed_MTspectrogram` argument naming
Both functions are declared `(FT_data, conditions, param, opt)` — plural `conditions`. NGL02 passes `condition` (singular). Positional, so it still works, but the inconsistency with `calculate_fireRate_general` (singular) is a footgun. Pick one and stick to it.

### 1.5 🟠 `calculate_neural_dynamics` synth branch is broken
Inside the function:
```matlab
nSpikes = poissrnd(firingRate(neuron) * trialDuration);
```
`firingRate` is a struct here (`firingRate.actual`, `firingRate.synth`), so `firingRate(neuron)` errors. Probably wants `firingRate.synth(neuron)`. Only fires if `opt.neurDyn.synth=true`, which is the default in `postPhy_param`, so this *will* fire.

### 1.6 🟡 Dead `clear` list
```matlab
clear events EventRecord trialdef blob neurons spike
```
`EventRecord` is never loaded in NGL02 (left over from NGL01). `condition` and `fireRate` *are* loaded but not cleared — they leak across sessions, which is exactly the kind of bug the `clear` was meant to prevent.

---

## 2. Multi-area inconsistencies (the priority lane)

### 2.1 🔴 No multi-area awareness in the spike branch
NGL01 explicitly loops over `input.areaMap.uniqueAreas` for both Kilosort and Bombcell, creating per-area `optArea` copies. NGL02 ignores this entirely:
- `loadSpikes(opt)` reads `opt.KSfolder` (singular). In multi-area mode `opt.KSfolder` points to the default `kilosort4/` subfolder, which has nothing in it because KS4 was run into per-area folders (`opt.KSfolders.NCL`, etc.).
- All downstream files (`spike.mat`, `neurons.mat`, `fireRate.mat`) are saved to a single location per session, with no area discriminator.

The NGL01 pattern to mirror:
```matlab
if isfield(input, 'areaMap') && ~isempty(input.areaMap)
    for a = 1:numel(input.areaMap.uniqueAreas)
        areaName         = input.areaMap.uniqueAreas{a};
        optArea          = opt;
        optArea.KSfolder = opt.KSfolders.(areaName);
        % ... loadSpikes(optArea), sort2trials, fireRate, save with area suffix ...
    end
else
    % single-area path (current behaviour)
end
```

### 2.2 🔵 Output layout decision needed for multi-area
Two viable patterns:
- **Per-area files:** `spike_NCL.mat`, `neurons_NCL.mat`, `fireRate_NCL.mat`, etc. Mirrors KS output folder structure. Easy to inspect; tedious to load.
- **Nested struct:** single `spike.mat` containing `spike.NCL = ...`, `spike.STR = ...`. Easier to consume downstream; rewrites the entire file every area iteration.

I lean toward per-area files (matches NGL01's `preprocessing/<session>/<Area>/` layout). Needs a sign-off before any code changes.

### 2.3 🟠 `opt.shank` / `opt.mapkey` overlaps with `input.Areas`
`postPhy_param.m` defines:
```matlab
opt.shank  = [1, 2, 3];
opt.mapkey = {'NCL', 'NCL', 'STR'};
```
This is the *old* way of associating shanks with area labels, and it duplicates `input.Areas` (which is the canonical NGL01 mechanism). `loadSpikes` uses `opt.mapkey(spike.shank{cl})` to tag each cluster with an ROI. We should either:
- Drop `opt.shank`/`opt.mapkey` from `postPhy_param.m` and derive the equivalent from `input.Areas` / `input.areaMap`.
- Keep them but have `set_default` validate that they agree with `input.Areas` when both are present.

Either way it shouldn't be two parallel sources of truth.

### 2.4 🟠 `opt.spikeSorted` is per-session, not per-area
`prepforsession` sets `opt.spikeSorted = <spikeSorted>/<subject>/<session>`. For multi-area, do we want `<.../<session>/<Area>/`? Tied to the output-layout decision in 2.2.

---

## 3. `opt.analysis` vs `input.analysis` confusion

NGL02 mixes the two:
- `opt.analysis` = `<study>/data/analysis/<subject>/<session>` (per-session, set by `prepforsession`). Used at lines 33, 53, 67, 71, 79, 90, 94, 97, 111, 116, 119.
- `input.analysis` = `<study>/data/analysis` (study-wide, set by `set_default`). Used at lines 140, 141 only.

Both lines that use `input.analysis` are in the LFP block and load `data_all.mat` / `allconditions` — i.e., an aggregated, study-wide output that was apparently produced by something like `collect_allData`. This is project-specific aggregator output and probably shouldn't be in the per-session loop at all.

**Action:** either move the `data_all.mat` lookup outside the loop (and gate it behind a flag), or drop it entirely and require `condition` to be loaded from the per-session `condition.mat` (which it already is at line 74).

---

## 4. Convention drift from NGL01

### 4.1 🟡 Doc header
NGL02 has only a 1-paragraph comment header. NGL01 has the full `PURPOSE / USAGE / REQUIRED WORKSPACE VARIABLES / PIPELINE / OUTPUTS / DEPENDENCIES / Last modified` template. Bring NGL02 to parity.

### 4.2 🟡 Repetitive load idiom
The pattern
```matlab
if ~exist('events','var'), load(fullfile(opt.trialSorted, "events.mat")); end
```
appears ~8 times. NGL01 had similar boilerplate that we condensed. A `loadIfMissing(varname, path)` helper in `functions/utils/` would clean this up considerably.

### 4.3 🟡 `addtime` semantics
`opt.addtime` defaults to `0` in `default_opt.m`, `500` in the SetAndRunMe template. NGL01 uses it during trialdef generation. NGL02 doesn't touch it but the LFP fallback regenerates `trialdef` via `trialdefGen` — so the value at NGL02 time must equal the value at NGL01 time, or the regenerated trialdef will be slightly off. This is exactly the kind of silent drift our new `preprocInfo` overlay is meant to prevent (`addtime` is in `applyPreprocInfo`'s default field list — good).

### 4.4 🟠 Trialdef regeneration in the LFP block
Lines 145-152:
```matlab
if ~exist('trialdef','var')
    load(fullfile(opt.analysis, "events.mat"));
    if exist('events','var')
        [~, trialdef, ~] = trialdefGen(events, opt, 1);
        save(fullfile(opt.trialSorted, "trialdef.mat"), 'trialdef');
    end
end
```
Two problems:
1. `events.mat` lives in `opt.trialSorted`, not `opt.analysis` (mistake).
2. Regenerating `trialdef` from `events` (rather than from `EventRecord`) is a degraded path — `trialdefGen` was designed to take `EventRecord`. The output may not match what NGL01 produced.

If `trialdef.mat` is missing, that's a real problem in the preprocessing and we should error rather than silently regenerate.

### 4.5 🔵 NGL02 only loops over `input.sessions(x).list`; never validates that NGL01 actually ran for each session
There's no check that `preprocessing/<subj>/<session>/kilosort4/` exists, that `cluster_info.tsv` from Phy is present, etc. If the user accidentally runs NGL02 on a session that wasn't curated, they'll get a confusing crash deep inside `loadSpikes`. A pre-flight check (mirroring `checkAnalysisCode`) would be cheap and helpful.

---

## 5. Project-specific code leakage

NGL02 still has hard-coded project-specific bits that should live in project hooks, not in the toolbox script:

- Line 155 (commented): `trialdef{2,1} = ceil(trialdef{2,1}/32);` — `chgDtctPCue`-specific.
- Line 176: `param.testname = 'ASL_Clean_Final_Correct';` — `SocialLearning`/`ASL`-specific magic string.
- Lines 191-199 (commented): `vFLIP_NGL` call with hard-coded `laminaraxis`, `freqaxis` — Miller-lab specific.
- Line 84-85 (commented): `calculate_fireRate_extintion` — Extintion-specific.
- Comments throughout referencing specific paradigms (`SocialLearning`, `chgDtctPCue`, `S3-Extintion`).

These belong either deleted, or wrapped in a project-hook mechanism (`functions/<project>/hook_postSpike.m` or similar), or moved to user-side scripts that the user opts into.

---

## 6. Downstream function issues

### 6.1 `loadSpikes.m`
- 🟠 Inline defaults for `opt.spparams`, `opt.getwF`, `opt.isibins` — should be moved to `default_opt` (post-Phy section) or `postPhy_param`. Once NGL02 starts calling `set_default`, those defaults should already be present and not duplicated here.
- 🔴 Reads `opt.gwfparams.wfWin`, `opt.gwfparams.nWf` without `isfield` checks — crashes if `postPhy_param` wasn't run.
- 🟠 References `opt.mapkey` (see 2.3) — same fate.
- 🟡 Commented `if ~exist(...)` wrapper around the whole body — dead code.
- 🟠 No docstring; doesn't follow NGL01 function-header convention.

### 6.2 `sort2trials.m`
- 🟠 The non-cell-trialdef branch ("Social interactions") is project-specific and is wedged into a general-purpose function. Should be split into `sort2trials` (general) and `sort2trials_blob` (social), or gated by an explicit mode argument.
- 🟠 `neurons_FT` always returns empty struct. The FT-spike-trial-parsing block is commented out. Either fix it or drop the second output.
- 🟡 `% TODO make FT spike trial parsing work` — open work item.
- 🟠 No docstring; doesn't follow NGL01 function-header convention.

### 6.3 `calculate_fireRate_general.m`
- 🟠 Inline defaults for `param.*` — should live in `postPhy_param`.
- 🟠 `param.levels = 1; ... for p = 1:param.levels` — the "levels" concept is unfinished, with all the project-specific examples commented out.
- 🟠 Calls `plot_multi_fireRate` directly — analysis tied to plotting, opposite of separation in NGL01. Either move plotting to NGL03 or make it `if param.plot, ... end` (it has the flag but plot is true by default).
- 🟠 Hard-codes `'allInitiated'` and `condition.aborted` — implies a condition struct shape that may not be universal.
- 🟠 No docstring.

### 6.4 `calculate_neural_dynamics.m`
- 🔴 Synth branch is broken (see 1.5).
- 🟠 Hardcodes `neurons.itiOn` — assumes `itiOn` is always an alignment. Should iterate `opt.alignto` like the rest.
- 🟠 Inline defaults; no docstring.

### 6.5 `trialdefGen.m`
- 🟡 `useevents` branch duplicates significant logic from the no-`useevents` branch. Could be unified.
- 🟠 References `opt.newEvent` which isn't in `default_opt` — undefined in the standard path.
- 🟠 Hard-coded project FIX (`-1000 is an Extintion arena FIX for 079's #1-11`).

### 6.6 LFP function signatures
- 🔴 `artifact_detRej_lfp` signature mismatch (see 1.3).
- 🟠 Plural `conditions` in `continous_MTspectrogram` / `trialparsed_MTspectrogram` (see 1.4).
- 🟠 No multi-area awareness anywhere in LFP path (probably outside our immediate scope but worth noting).

---

## 7. Recommended next-step task list

Prioritised, derived from above. Each one is a candidate work item.

| # | Item | Severity | Effort |
|---|---|---|---|
| A | Fix `prepforsession` call signature in NGL02 (1.1) | 🔴 | trivial |
| B | Fix `conditions` → `condition` typo (1.2) | 🔴 | trivial |
| C | Fix `artifact_detRej_lfp` call signature (1.3) | 🔴 | trivial |
| D | Decide singular vs plural `condition(s)` convention; align signatures (1.4) | 🟠 | small |
| E | Fix `calculate_neural_dynamics` synth branch (1.5) | 🔴 | trivial |
| F | Clean up `clear` list to match what NGL02 actually uses (1.6) | 🟡 | trivial |
| G | Multi-area design decision: per-area files vs nested struct (2.2) | 🔵 | design |
| H | Implement multi-area loop in spike branch, mirroring NGL01 (2.1, 2.4) | 🔴 | medium |
| I | Reconcile `opt.shank`/`opt.mapkey` with `input.Areas` (2.3) | 🟠 | small |
| J | Drop or move `input.analysis` references (data_all.mat) out of the loop (§3) | 🟠 | small |
| K | NGL01-style doc header on NGL02 (4.1) | 🟡 | small |
| L | `loadIfMissing` helper to dedupe load idiom (4.2) | 🟡 | small |
| M | Replace trialdef regeneration with hard error if missing (4.4) | 🟠 | small |
| N | Pre-flight check that NGL01 outputs exist for each session (4.5) | 🟠 | small |
| O | Strip / hook-ify project-specific code (§5) | 🟡 | medium |
| P | Move inline defaults from analysis functions into `default_opt` / `postPhy_param` (§6) | 🟠 | medium |
| Q | Add NGL01-style docstrings to `loadSpikes`, `sort2trials`, `calculate_*` (§6) | 🟡 | medium |
| R | Split / gate the "Social interactions" branch in `sort2trials` (6.2) | 🟠 | small |
| S | Decouple `calculate_fireRate_general` from `plot_multi_fireRate` (6.3) | 🟠 | small |

Suggested order: **A → B → C → E → F** (zero-risk bug fixes first), then **G** (design decision), then **H → I → J** (multi-area + structural cleanup), then **K → L → M → N** (convention hygiene), then **O → P → Q → R → S** (deeper cleanup).

---

## 8. Open design questions for you

These are the questions I need answered before any code changes:

1. **Multi-area output layout** — per-area files (`spike_NCL.mat`) or nested struct (`spike.NCL`)? (item G)
2. **`opt.shank` / `opt.mapkey` future** — keep with cross-check, or drop in favour of `input.Areas`? (item I)
3. **`data_all.mat` lookup** — kill it, or move it to a separate aggregation stage outside NGL02's loop? (item J)
4. **`condition` vs `conditions` convention** — which one becomes the standard? (item D)
5. **Project-specific code** — delete-on-sight, or formalise a hook mechanism? (item O)

Once those five are answered, the rest is mostly mechanical.
