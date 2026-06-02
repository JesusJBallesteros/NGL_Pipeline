# `calculate_neural_dynamics.m` — Audit and Forward Plan

**Scope:** explain what the function does, what it's *trying* to do, where the implementation falls short of the goal, and what alternative methods are worth supporting next.
**Author:** audit pass, 29.05.2026.

---

## 1. What the function is for

The stated goal (per the original one-liner) is "neural dynamics analysis using PCA, t-SNE, or UMAP." In practice this means:

> *Given a session's curated clusters and the trial-aligned spike times, build a low-dimensional embedding of population activity so the user can see structure across trials — drift, condition clusters, stimulus geometry — without staring at hundreds of per-cluster heat-maps.*

That's a legitimate and common goal in systems-neuroscience analysis. It's the same goal that motivates jPCA, dPCA, GPFA, LFADS and the visual tools built around them. The function is a stepping-stone toward that capability, currently implemented at the "exploratory prototype" level.

## 2. What it actually does today

Walking the file top-to-bottom:

1. Reads `numNeurons = size(neurons.itiOn, 1)`, `numTrials = size(neurons.itiOn{1,1}, 1)`, and a `trialDuration` = median trial length in seconds (from `trialdef{2,1}`).
2. Builds `timeBins = 0 : opt.neurDyn.binSize : trialDuration + binSize`. `binSize` defaults to **0.1 s**.
3. **Branches on `opt.neurDyn.synth`.**
   - **synth = true (the default!):** generates a homogeneous-rate Poisson surrogate for every (neuron, trial) using a per-neuron mean rate drawn uniformly between `min(fireRate.sps, [], 'all')` and `max(...)`. Spikes are uniformly distributed in `[0, trialDuration]`. This becomes `spTimes.synth`.
   - **synth = false:** `spTimes = neurons` (i.e., the input struct).
4. Bins into a `(numNeurons × numTimeBins × numTrials)` array via `histcounts`.
5. Reshapes to `(numTrials × (numNeurons*numTimeBins))` — **each trial becomes one row.**
6. Runs PCA / t-SNE / UMAP on that matrix and keeps the first 3 components.
7. Draws two figures: a `scatter3` of the 3-D embedding coloured by trial index, and a "trajectories" figure that calls `plot3` with one point per trial.
8. Saves a `.fig` to the current working directory.

The output argument `neuralDynamics` is declared but never assigned.

## 3. Fitness to purpose — what's broken or misaligned

### 3.1 Severity-ranked findings

| # | Severity | Issue |
|---|---|---|
| A | 🔴 critical | `opt.neurDyn.synth` defaults to **true**, so the function embeds a Poisson surrogate by default. Whatever the user sees is the dim-reduction of noise, not of the real session. |
| B | 🔴 critical | The `else` branch (real data) is dead — `spTimes = neurons` is then indexed as `spTimes{neuron, trial}` inside the binning loop, but `neurons` is a struct with `.itiOn`, not a 2-D cell. This will error the moment `synth=false`. |
| C | 🔴 critical | Hardcodes `neurons.itiOn`. In multi-area runs (post-recent refactor), `neurons` has shape `neurons.<Area>.itiOn`. The function will fail in multi-area mode. |
| D | 🟠 logic | Units inconsistent. `neurons` stores spike times in **ms**; the synth surrogate generates them in **seconds**; `timeBins` is in **seconds**. The dead "real-data" path would also bin ms spikes with second-scaled edges, producing all zeros. |
| E | 🟠 logic | The "trajectory" figure does not draw trajectories. `plot3(x, y, z)` with three scalars draws a single point, not a line — and `reducedData` has no time dimension to draw a line through. To get trajectories you have to keep time as a separate axis (don't reshape it away). |
| F | 🟠 logic | `reducedData(:, 1:3)` from PCA can blow up if `pca` returns fewer than 3 components, which happens whenever `numTrials < numNeurons*numTimeBins` — i.e., almost always. No guard. |
| G | 🟠 logic | `fireRate` is taken as input but only used to bound the surrogate. The real `fireRate.sps` is already a binned, rate-normalised `[Ntrials × Nbins]` matrix per cluster — exactly the input dim-reduction expects. The function does an independent rebinning from raw spike times, with a different bin size, in inconsistent units. |
| H | 🟠 logic | `histcounts(spikes, [timeBins, trialDuration])` appends a duplicate edge (`timeBins` already ends at `≥ trialDuration`), giving an off-by-one bin. Compensates for a count-vs-edge mistake elsewhere; cleaner to use `histcounts(spikes, timeBins)` and accept `numTimeBins-1` bins. |
| I | 🟡 product | "Trial-state PCA" is a valid analysis, but it's not what *neural dynamics* usually means in the literature (which is time-resolved trajectories). The function name promises something stronger than the implementation delivers. |
| J | 🟡 product | Colouring by trial index reveals drift over trials but hides condition structure (correct/incorrect, stimulus type, block phase). For most paradigms condition colouring is the more interesting view. |
| K | 🟡 hygiene | t-SNE perplexity = 10 is hard-coded; UMAP path silently `addpath('umap')` (relative). |
| L | 🟡 hygiene | `savefig([opt.neurDyn.method '_neural_dynamics.fig'])` lands in `pwd`, not in `opt.analysis/plots/`. |
| M | 🟡 hygiene | The declared output variable `neuralDynamics` is never assigned. NGL02 then saves it to `.mat`, which currently relies on undefined behaviour. |
| N | 🟡 hygiene | No defaults are validated; `set_default` doesn't know about `opt.neurDyn.binSize`/`method`/`synth`. |

### 3.2 Conceptual misalignment

Beyond the bugs, there's a deeper mismatch. The standard "neural population dynamics" pipeline in the field is:

1. Smooth each neuron's spike train with a kernel (Gaussian, ~50 ms σ is typical).
2. Trial-average per condition (or keep single-trial if using GPFA / LFADS).
3. Stack neurons into a `(N_neurons × N_timepoints × N_conditions)` tensor.
4. Reduce along the neuron axis using a method appropriate to the question:
   - **PCA** → variance-explaining axes; trajectory through time per condition.
   - **jPCA** → planes with strong rotational dynamics; great for motor / decision tasks.
   - **dPCA** → axes demixed by task variable (stimulus / decision / time / interaction).
   - **GPFA** → smooth single-trial trajectories.
   - **LFADS** → nonlinear single-trial inference of latent dynamics.
5. Plot trajectories *through time*, coloured by condition.

The current function collapses time into a feature vector and reduces over trials. That's a "trial-similarity embedding," useful for spotting clusters of trials with unusual population states, but it doesn't show neural dynamics in the traditional sense. The two are complementary; the function should make the distinction clear or pick one.

## 4. Research summary — what the field uses

### 4.1 Methods worth supporting

#### PCA on trial-averaged smoothed rates (current standard for "first look")
- Smooth spikes with Gaussian kernel (σ ≈ 50 ms).
- Average across trials within each condition.
- Stack into `(N_neurons × (N_time × N_conditions))` and run `pca`.
- Reshape PCs back to `(N_PCs × N_time × N_conditions)`, plot one line per condition through time.
- **MATLAB:** trivial with built-in `pca`.
- **Why:** still the most common first-pass analysis in 2024-2025 population-dynamics papers.

#### jPCA (Churchland, Cunningham, et al., 2012)
- Finds linear combinations of PCs that show *rotational* dynamics in 2-D planes.
- Great for motor planning / preparation / perceptual-decision tasks where trajectories rotate.
- **MATLAB:** original code from the Churchland Lab is the reference; community Python port at [bantin/jPCA](https://github.com/bantin/jPCA).
- **Why add it:** small, self-contained, MATLAB-native, well-suited to NGL paradigms where time-locked dynamics are the question.

#### dPCA — demixed PCA (Kobak & Brendel, 2016)
- Decomposes population variance into axes attributable to specific task variables (stimulus, decision, time, and their interactions).
- Requires a balanced condition structure.
- **MATLAB:** Kobak Lab reference code (their original repo is the canonical MATLAB implementation).
- **Why add it:** the most interpretable dim-reduction for structured paradigms (colour wheel, S3 extinction). Each axis answers a specific question like "is this neuron tuned to stimulus identity?"

#### GPFA — Gaussian-Process Factor Analysis (Yu, Cunningham, Churchland, Sahani, Shenoy, 2009)
- One-stage method: smoothing and dim-reduction done jointly via a GP model on factors.
- **Single-trial** trajectories — no need to average across trials.
- **MATLAB:** Byron Yu's `DataHigh` GUI ([DataHigh — Yu Lab software](https://users.ece.cmu.edu/~byronyu/software/DataHigh/)) and the standalone GPFA codepack from the same page. Also community implementation at [wrongu/gpfa](https://github.com/wrongu/gpfa) and [aecker/gpfa](https://github.com/aecker/gpfa).
- **Why add it:** the right tool when single-trial variability matters and you don't want to throw it away by averaging.

#### LFADS — Latent Factor Analysis via Dynamical Systems (Pandarinath et al., 2018)
- RNN-based, fully nonlinear single-trial inference.
- Strong results but heavy: Python + TensorFlow, GPU strongly recommended.
- Reference: [Pandarinath et al., *Nat. Methods* 2018](https://www.nature.com/articles/s41592-018-0109-9). Generalised follow-up: [Keshtkaran et al., 2022, PMC9825111](https://pmc.ncbi.nlm.nih.gov/articles/PMC9825111/).
- **Why we'd defer it:** Python pipeline, not a quick add to a MATLAB toolbox. Worth pointing users at if they need single-trial inference quality beyond GPFA.

### 4.2 What stays the same regardless of method

- **Smoothing.** Replace `histcounts` with a Gaussian kernel convolution before any dim-reduction. ~50 ms σ is a defensible default; this is what every neural-dynamics paper does. Reference: [Kernel methods on spike train space for neuroscience, Park & Brockmeier, 2013, arxiv:1302.5964](https://arxiv.org/pdf/1302.5964); [Yu et al. 2009, *J. Neurophysiol.*](https://journals.physiology.org/doi/full/10.1152/jn.90941.2008).
- **Trial averaging by condition.** Even for single-trial methods, the comparison view is per-condition PSTH trajectories.
- **Time stays a separate axis.** Reduction is across neurons, *not* across time.
- **Colour by condition, not trial index** (with an optional second view coloured by time for drift checks).

### 4.3 What we already have that helps

- `fireRate.sps{c, 1}` is already `(N_trials × N_bins)` per cluster, with consistent binning across the toolbox (param.binSize / stepSz / interval). This is the natural input to a population-dynamics function. Stacking across clusters gives the `(N_neurons × N_bins × N_trials)` tensor every method above needs as a starting point.
- `param.ROI`, `param.HumanLabel`, `param.bc_unitType` are now available per cluster (from the recent labels pass), so the dim-reduction can pre-filter to "good" units only, or weight by `bc_unitType`.

## 5. Recommended fixes for the existing function (small, near-term)

These keep the function's current scope ("quick dim-reduction view of trial-level population state") but make it correct and useful.

1. **Default `opt.neurDyn.synth = false`.** Real data is the primary use case. Move the synth branch behind a flag whose name makes clear it's a null/control: e.g. `opt.neurDyn.runNullSurrogate`, off by default.
2. **Repair the real-data path.** Build the per-trial-per-cluster matrix directly from `fireRate.sps` — it's already binned, in consistent units, in the right shape. Drop the manual `histcounts` block.
3. **Multi-area + alignment-generic.** Take an explicit `alignment` argument (default `opt.alignto{1}`). Make NGL02 call this once per area in multi-area mode, exactly the way it does for `fireRate.<Area>` now.
4. **Compute `numComponents` safely.** Cap at `min(3, size(reducedData, 2))`; warn if fewer than 3.
5. **Real trajectory plot.** Optionally produce a *time-resolved* trajectory: reshape `fireRate.sps` to keep time as an axis, average across trials per condition, project per time-bin into the same 3-D space, and `plot3` a line per condition.
6. **Save to `opt.analysis/plots/neural_dynamics/`** with a session-scoped filename. Don't write to `pwd`.
7. **Assign `neuralDynamics`.** Return at least `.reducedData`, `.method`, `.explainedVariance` (for PCA), `.timeAxis` (when applicable), so the `.mat` file has something usable.
8. **Move all `opt.neurDyn.*` defaults into `default_opt.m`.** Validate in `set_default`.
9. **Drop the hard-coded UMAP `addpath('umap')`.** Check `exist('run_umap','file')` and error with a clear "install UMAP from File Exchange #71902" message.

These are an afternoon's work. They make the function honest about what it does and remove the silent-Poisson-surrogate trap.

## 6. Forward plan — what to add next

Concrete proposal, in order of effort-to-payoff.

| Tier | Method | New function | Effort | Pre-req |
|---|---|---|---|---|
| 1a | PCA on trial-averaged smoothed rates | `calculate_neural_pca.m` | small | smoother helper (see below) |
| 1b | jPCA | `calculate_neural_jpca.m` | small-medium | import Churchland Lab codepack into `toolboxes/` |
| 2  | dPCA | `calculate_neural_dpca.m` | medium | import Kobak codepack; needs balanced conditions per project |
| 3  | GPFA | `calculate_neural_gpfa.m` | medium-large | import Yu Lab GPFA codepack |
| 4  | LFADS hook | — | large | separate Python pipeline; defer |

Supporting infrastructure to add once and reuse:

- `functions/analysis/smooth_spikes.m` — Gaussian-kernel rate estimate, takes `neurons.<align>{c,t}` or `fireRate.sps{c}` and a σ; returns smoothed `(N_neurons × N_time × N_trials)` tensor.
- A small `param.dimRed` block in `default_opt` covering common knobs (`smoothSigma`, `nComponents`, `condition` for colouring), validated in `set_default`.

After this, the current `calculate_neural_dynamics.m` could be either:
- **Renamed** to `calculate_neural_trialEmbedding.m` to honestly reflect what it does (per-trial state embedding), kept as a sibling of the trajectory functions; or
- **Retired** in favour of the new family, with the trial-state view becoming an optional mode of `calculate_neural_pca.m`.

I'd recommend the rename — the trial-state view is genuinely useful for drift detection and trial-clustering questions; it just shouldn't be called "dynamics."

## 7. Sources

- Yu, Cunningham, Santhanam, Ryu, Shenoy, Sahani (2009). *Gaussian-Process Factor Analysis for Low-Dimensional Single-Trial Analysis of Neural Population Activity.* [J. Neurophysiol., PMC2712272](https://pmc.ncbi.nlm.nih.gov/articles/PMC2712272/) · [Journal of Neurophysiology version](https://journals.physiology.org/doi/full/10.1152/jn.90941.2008).
- Churchland et al. (2012). *Neural population dynamics during reaching.* [PLOS Comp Biol, 1005175](https://journals.plos.org/ploscompbiol/article?id=10.1371%2Fjournal.pcbi.1005175). jPCA reference.
- Pandarinath et al. (2018). *Inferring single-trial neural population dynamics using sequential auto-encoders (LFADS).* [Nature Methods](https://www.nature.com/articles/s41592-018-0109-9) · [bioRxiv preprint](https://www.biorxiv.org/content/10.1101/152884.full.pdf) · [PMC6380887](https://pmc.ncbi.nlm.nih.gov/articles/PMC6380887/).
- Keshtkaran et al. (2022). *A large-scale neural network training framework for generalized estimation of single-trial population dynamics.* [PMC9825111](https://pmc.ncbi.nlm.nih.gov/articles/PMC9825111/).
- Park & Brockmeier (2013). *Kernel Methods on Spike Train Space for Neuroscience: A Tutorial.* [arxiv:1302.5964](https://arxiv.org/pdf/1302.5964).
- Unnikrishnan & Jacob (2015). *Kernel bandwidth optimization in spike rate estimation.* [PMC2940025](https://ncbi.nlm.nih.gov/pmc/articles/PMC2940025).
- *Comparing Open-Source Toolboxes for Processing and Analysis of Spike and Local Field Potentials Data.* [Frontiers in Neuroinformatics](https://www.frontiersin.org/journals/neuroinformatics/articles/10.3389/fninf.2019.00057/full) · [PMC6682703](https://pmc.ncbi.nlm.nih.gov/articles/PMC6682703/).
- **DataHigh** (Cowley, Smith, Yu): [Yu Lab software page](https://users.ece.cmu.edu/~byronyu/software.shtml). MATLAB GUI for population dim-reduction.
- **Elephant** (Python): [GPFA tutorial](https://elephant.readthedocs.io/en/latest/tutorials/gpfa.html). Useful reference for binning / smoothing conventions.
- **jPCA implementations**: [bantin/jPCA (Python)](https://github.com/bantin/jPCA), [alexmorley/jPCA (Julia)](https://github.com/alexmorley/jPCA).
- **GPFA implementations (MATLAB)**: [wrongu/gpfa](https://github.com/wrongu/gpfa), [aecker/gpfa](https://github.com/aecker/gpfa), [karthikcl2590/NeuralTraj](https://github.com/karthikcl2590/NeuralTraj).
