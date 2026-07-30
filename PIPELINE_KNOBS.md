# Pipeline knobs

A running inventory of every tunable parameter in the classydms pipeline, from raw file → classified prediction. Intended as reference for:

1. **Ad-hoc tuning** — quickly finding what to try when a step underperforms.
2. **Future comprehensive optimization** — designing an end-to-end sweep (probably a small factorial or Bayesian-optimization run over the case/control simulation harness) to pick reasonable defaults for a new cohort.

For each knob:
- **Where**: function + arg name
- **Default**: current value in the code (as of the most recent commit)
- **Type**: numeric / integer / categorical / boolean
- **Reasonable range**: sweep bounds worth exploring
- **Sweep priority**: how likely it is to matter (H = high, M = medium, L = low)
- **Cost**: relative expense of one evaluation (cheap / moderate / expensive)
- **Notes**: rationale, interactions with other knobs, references

The "reasonable range" is a starting bracket for sweeps — not a hard limit. Priorities reflect our best current guess about impact on downstream classification AUC; they should be updated as we accumulate sweep evidence.

---

## Stage 1 — Raw file → intensity matrix

| Knob | Where | Default | Type | Range | Priority | Cost | Notes |
|---|---|---|---|---|---|---|---|
| `pattern` | `load_dms_directory` | `POS.*\.(xls\|txt\|tsv\|csv)$` | regex | POS/NEG variants | L | cheap | Selects channel; POS is standard for classydms |

---

## Stage 2 — Preprocessing (`process_one_sample`, `baseline_basement`, trim)

| Knob | Where | Default | Type | Range | Priority | Cost | Notes |
|---|---|---|---|---|---|---|---|
| `sg_p` | `process_one_sample` | 3 | int | 2–5 | M | cheap | Savitzky-Golay polynomial order. Higher = better preservation of sharp peaks but more noise. Never formally tuned. |
| `sg_n` | `process_one_sample` | 21 | int (odd) | 11–51 | M | cheap | SG window length. `data_load_process_pipeline.Rmd` uses 51; `resnet_script.R` uses 21. **Inconsistent across pipelines — standardize before sweeping.** |
| `asls_lambda` | `process_one_sample` | 1e6 | numeric (log-scale) | 1e4–1e9 | L | cheap | AsLS baseline smoothness. Higher = smoother baseline. Rarely need to tune. |
| `asls_p` | `process_one_sample` | 0.01 | numeric | 0.001–0.05 | L | cheap | AsLS asymmetry. Controls how much we assume the baseline sits below the signal. |
| `asls_niter` | `process_one_sample` | 10 | int | 5–20 | L | cheap | AsLS iterations. Converges quickly. |
| `basement_thr` | `baseline_basement` | 0.005 | numeric | 0.001–0.02 | H | cheap | Dust threshold. Zeros pixels below this value. Directly affects sparsity of both real and synthetic samples. `synthetic_quality_check` sparsity metric depends on this matching. |
| `q_lo`, `q_hi` | `pooled_trim_bounds` | 0.05, 0.95 | numeric | 0.01–0.10 / 0.90–0.99 | M | cheap | Occupancy quantiles for cross-sample trim bounds. Aggressive trim → smaller images → faster training but risk losing peaks at margins. |
| `smooth_k_rt`, `smooth_k_cv` | `trim_bounds_from_occupancy` | 31, 9 | int | 15–51 / 5–15 | L | cheap | Occupancy smoothing window before threshold. Hand-chosen from visual inspection; never formally tuned. See `todos.R`. |
| `thr_rt`, `thr_cv` | `trim_bounds_from_occupancy` | 0.005, 0.01 | numeric | 0.001–0.02 | L | cheap | Occupancy threshold defining "populated" RT/CV pixels. Same status as above. |

---

## Stage 3 — Peak detection (`pick_peak_centers`, `estimate_peak_params`)

| Knob | Where | Default | Type | Range | Priority | Cost | Notes |
|---|---|---|---|---|---|---|---|
| `top_k` | `estimate_peak_params` | 150 (function); 2500 (vignette) | int | 100–5000 | M | moderate | Max peaks per sample. Too low biases intensity distribution toward brightest peaks; vignette bumps to 2500 to capture dim tail. |
| `frac_height` | `estimate_peak_params` | 0.25 (function); 0.5 (vignette) | numeric | 0.2–0.7 | M | moderate | Fraction of peak height defining "the peak". 0.5 = FWHM. Directly sets the observed sigma distributions. |
| `min_sep_rt` | `pick_peak_centers` | 8 | int | 4–15 | H | moderate | NMS separation in RT pixels. Loose → over-merge co-eluting peaks; tight → fragment single peaks into duplicates. See `todos.R` for the min_sep sweep task. |
| `min_sep_cv` | `pick_peak_centers` | 2 | int | 1–5 | H | moderate | Same for CV. Aggressive at 2 given CV pixel resolution. |
| `eps` | `pick_peak_centers` / `robust_eps` | `robust_eps(Z)` | numeric | see `robust_eps` params | H | cheap | Intensity floor. Current `robust_eps` = 0.01 × Q95(positive). Not tied to noise directly. Follow-up: hybrid `max(3 × noise_sd, 0.01 × Q95)`. |
| `q`, `frac` | `robust_eps` | 0.95, 0.01 | numeric | 0.9–0.99 / 0.005–0.05 | M | cheap | Determine the adaptive eps. |
| `max_radius_rt`, `max_radius_cv` | `peak_width_at` | `Inf`, `Inf` | int | 15–50 / 5–15 | L | cheap | Bound on width-measurement walk. Prevents accidental jumps into neighboring peaks; usually not binding. |

---

## Stage 4 — Synthetic data generation (`generate_one_synthetic`, `generate_synthetic_dataset`)

| Knob | Where | Default | Type | Range | Priority | Cost | Notes |
|---|---|---|---|---|---|---|---|
| `location_mode` | `generate_one_synthetic` | `"empirical"` | categorical | empirical / marginal | H | cheap | Bootstrap observed `(rt, cv)` (default) vs independent Normal marginals. Empirical preserves hotspot structure; marginal scatters uniformly. |
| `attribute_mode` | `generate_one_synthetic` | `"joint"` | categorical | joint / marginal | H | cheap | Bootstrap observed `(σ_rt, σ_cv, intensity)` triple (default) vs three independent log-normals. Joint preserves wide↔bright correlation. |
| `location_jitter_rt` | `generate_one_synthetic` | 2 | numeric (pixels) | 0–10 | M | cheap | SD of per-peak location jitter on RT axis. Represents instrument reproducibility. |
| `location_jitter_cv` | `generate_one_synthetic` | 1 | numeric (pixels) | 0–3 | M | cheap | Same, CV axis. |
| `size_jitter` | `generate_one_synthetic` | 0.15 | numeric | 0.05–0.6 | H | cheap | SD of log-normal per-peak size multiplier. Was 0.6 (too loose); tuned down to 0.15 based on QC of synthetic-real match. Also documented in getting-started vignette. |
| `add_noise` | `generate_one_synthetic` | TRUE | bool | TRUE / FALSE | H | cheap | Whether to add background noise. FALSE = reconstruction AE; TRUE = denoising AE. |
| `noise_scale` | `generate_one_synthetic` | 1.0 | numeric | 1–20 | H | cheap | **New** multiplier on `params$noise$sd`. Larger = harder denoising task = stronger regularization. Vincent et al. sweet spot around 3–10 for typical AEs; classydms untuned. |
| `dust_threshold` | `generate_synthetic_dataset` | 0 | numeric | 0 / match `basement_thr` | H | cheap | Zero pixels below this on synthetic samples BEFORE normalization, matching real preprocessing. Should equal `basement_thr` for consistent sparsity. |

---

## Stage 5 — Autoencoder architecture

| Knob | Where | Default | Type | Range | Priority | Cost | Notes |
|---|---|---|---|---|---|---|---|
| `stem_stride_rt` | `dms_encoder` | 4 | int | 2–8 | M | expensive | Initial downsampling on RT axis. Larger = smaller feature map = faster but coarser localization. 4 chosen to handle ~1400px RT dim. |
| `stem_stride_cv` | `dms_encoder` | 1 | int | 1–2 | M | expensive | Same, CV axis. 1 preserves CV resolution since CV has fewer pixels. |
| `target_H`, `target_W` | `compute_pad_targets` / `pad_to_target` | max cohort dims | int | data-dependent | L | moderate | Padded input size. Determined by cohort, not user-chosen. |
| encoder channel widths | `dms_encoder` (baked in) | (32, 64) | int | 16–128 | M | expensive | Channel counts per stage. Currently fixed; would need code change to sweep. |
| bottleneck dim | `dms_encoder` (baked in) | 64 | int | 16–256 | H | expensive | Encoder feature dimension after GAP. Sets capacity of the learned representation. Currently fixed; code change needed to sweep. |
| number of ResNet blocks per stage | `dms_encoder` (baked in) | 2 | int | 1–4 | L | expensive | Depth. |

---

## Stage 6 — Pretraining (`pretrain_denoising_online`, `pretrain_autoencoder`)

| Knob | Where | Default | Type | Range | Priority | Cost | Notes |
|---|---|---|---|---|---|---|---|
| `steps_per_epoch` | `pretrain_denoising_online` | 500 | int | 100–2000 | M | expensive | Batches per epoch. Determines samples-seen per epoch alongside `batch_size`. |
| `epochs` | `pretrain_denoising_online` | 30 | int | 10–100 | M | expensive | Total epochs. Watch validation loss for plateaus. |
| `batch_size` | `pretrain_denoising_online` | 32 | int | 8–128 | M | expensive | Batch size. Larger = smoother gradient but more memory. |
| `lr` | `pretrain_denoising_online` | 1e-3 | numeric (log-scale) | 1e-5 – 1e-2 | H | expensive | Adam learning rate. Standard sweep target. |
| `weight_decay` | `pretrain_denoising_online` | 1e-4 | numeric (log-scale) | 1e-6 – 1e-2 | M | expensive | L2 regularization on encoder weights. |
| `grad_clip` | `pretrain_denoising_online` | 1.0 | numeric | 0.5–5.0 | L | cheap | Gradient norm clip. Rare that this matters after the first epoch. |
| `norm_clamp` | `pretrain_denoising_online` | 10.0 | numeric | 5.0–100.0 | L | cheap | Loss clamp for stability. Belt-and-suspenders. |
| `val_n` | `pretrain_denoising_online` | 200 | int | 100–1000 | L | cheap | Synthetic validation set size for tracking. |
| `val_real` | `pretrain_denoising_online` | NULL | list of Z or NULL | data-dependent | H | cheap | Real held-out samples for real-val MSE. **Critical** — noise_scale sweeps depend on this metric. |
| `checkpoint_every` | `pretrain_denoising_online` | 10 | int | 5–20 | L | cheap | Just controls disk I/O. |
| `num_workers` | `pretrain_denoising_online` | 1 | int | 1–8 | L | moderate | Parallel data-gen workers. Speed only. |
| `device` | `pretrain_denoising_online` | auto | categorical | cpu / cuda / mps | L | cheap | Hardware target. |
| `seed` | `pretrain_denoising_online` | 42 | int | any | L | cheap | Reproducibility. Sweep for variance estimation. |

---

## Stage 7 — Classification (`cv_encoder_elastic_net`, `train_one_fold_pretrained`)

| Knob | Where | Default | Type | Range | Priority | Cost | Notes |
|---|---|---|---|---|---|---|---|
| `k` | `cv_encoder_elastic_net` | 5 | int | 3–10 | L | moderate | CV folds. Standard 5 is fine. |
| `alpha` | `cv_encoder_elastic_net` | 0.5 | numeric | 0–1 | H | cheap | Elastic net mixing: 0 = ridge, 1 = lasso. 0.5 balances L1 sparsity with L2 stability. Sweep for interpretability vs stability. |
| `lambda_sequence` | `cv_encoder_elastic_net` | glmnet default | numeric vector | | L | cheap | glmnet path; usually leave to auto. |
| `epochs_head`, `lr_head` | `train_one_fold_pretrained` | — | int / numeric | | M | expensive | Two-phase fine-tune parameters (head-only phase). |
| `epochs_backbone`, `lr_backbone` | `train_one_fold_pretrained` | — | int / numeric | | M | expensive | Full-model fine-tune. |
| `unfreeze_from_layer` | `train_one_fold_pretrained` | — | categorical | | M | expensive | How deep to unfreeze during phase 2. |

---

## Stage 8 — Augmentation (`augment_batch`)

| Knob | Where | Default | Type | Range | Priority | Cost | Notes |
|---|---|---|---|---|---|---|---|
| `rt_shift` | `augment_batch` | 5 | int (pixels) | 0–25 | H | cheap | Random RT shift augmentation. Larger → better invariance to drift. Sweep prior evidence: 25 outperformed 5 in one experiment. |
| `cv_shift` | `augment_batch` | 2 | int (pixels) | 0–5 | H | cheap | Same, CV. |
| `intensity_scale` | `augment_batch` | (0.9, 1.1) | numeric range | (0.8, 1.2) etc. | M | cheap | Random intensity rescale per batch. |
| `noise_sd` | `augment_batch` | 0.0 | numeric | 0.0–0.05 | M | cheap | Additional Gaussian noise added during training. Distinct from `noise_scale` (which is at synthesis time). |

---

## Stage 9 — Saliency (`compute_saliency_map`, `aggregate_saliency_masked`)

| Knob | Where | Default | Type | Range | Priority | Cost | Notes |
|---|---|---|---|---|---|---|---|
| `layer` | `compute_saliency_map` | encoder$layer3 | categorical | layer1..4 | M | moderate | Which encoder layer to hook for Grad-CAM. Deeper = coarser but semantically richer. |
| aggregation | `aggregate_saliency_masked` | mean | categorical | mean / max | L | cheap | How to combine per-sample saliency into a group-level map. |
| mask threshold | `aggregate_saliency_masked` | — | numeric | | L | cheap | Threshold below which saliency is zeroed. |

---

## Stage 10 — Case/control simulation (`simulate_case_control.R`)

| Knob | Where | Default | Type | Range | Priority | Cost | Notes |
|---|---|---|---|---|---|---|---|
| `n_biomarkers` | `holdout_biomarkers` | — | int | 1–20 | H | cheap | Level 1/2/3 set to 1 / 3 / 6. |
| `protect_radius_rt` | `holdout_biomarkers` | 15 / 1400 | numeric | matching `min_sep_rt` | M | cheap | Radius around biomarkers where pool peaks are also stripped. |
| `protect_radius_cv` | `holdout_biomarkers` | 3 / 40 | numeric | matching `min_sep_cv` | M | cheap | Same, CV. |
| `min_intensity_quantile` | `holdout_biomarkers` | NULL | numeric or NULL | 0.0–0.99 | H | cheap | Restrict biomarker selection to top X% brightest peaks. Level 1 uses 0.90, Level 2 uses 0.50, Level 3 uses NULL. |
| `case_prevalence` / `control_prevalence` | `build_biomarker_spec` | — | numeric (0–1) | 0–1 | H | cheap | Sets Bernoulli probability of biomarker firing per group. Level 1: 0.95/0.05; Level 2: 0.80/0.20; Level 3: 0.60/0.40. |
| `case_intensity_mult` / `control_intensity_mult` | `build_biomarker_spec` | — | numeric | 0.5–10.0 | H | cheap | Amplitude multiplier per group. Level 1: 5.0/1.0; Level 2: 1.5/1.0; Level 3: 1.2/1.0. |
| `n_cases`, `n_controls` | `simulate_case_control_cohort` | — | int | 10–500 | M | moderate | Cohort size. Small cohorts have small-sample variance; sweep to see AUC vs N curve. |

---

## Stage 11 — Biomarker recovery evaluation

| Knob | Where | Default | Type | Range | Priority | Cost | Notes |
|---|---|---|---|---|---|---|---|
| `match_radius_rt`, `match_radius_cv` | `evaluate_biomarker_recovery` | 8, 2 | int | matching `min_sep` | L | cheap | Pixel radius within which a saliency hotspot counts as matching a biomarker. |
| `top_k` | `evaluate_biomarker_recovery` | 10 | int | matching `n_biomarkers` | L | cheap | Number of local maxima to consider for precision/recall. |

---

## Metrics we optimize against

All sweeps ultimately point at these downstream measurements:

1. **Real-validation MSE** (`pretrain_denoising_online` with `val_real`) — encoder generalization from synthetic to real. Cheap; use for pretraining hyperparameter picks.
2. **Case/control AUC** (`cv_encoder_elastic_net$auc`) — classification performance. Primary target.
3. **Biomarker recovery per-pixel AUC** (`evaluate_biomarker_recovery`) — interpretability layer quality.
4. **`synthetic_quality_check` KS distances** (5 axes: total signal, mean positive, sparsity, n_peaks, pixel intensity) — synthetic-real distributional match. Diagnostic; not something to optimize directly, but a red flag if a knob drives these apart.

## Design of a comprehensive sweep (future)

When we tackle systematic optimization:

1. Start on the **case/control simulation harness** since it has ground truth for both AUCs. Real cohorts have only classification AUC and no ground-truth localization to score against.
2. Priority-H knobs first, using the case/control Level 2 setup as a fixed benchmark: `min_sep_rt`, `min_sep_cv`, `basement_thr`, `size_jitter`, `noise_scale`, `lr`, `alpha`, `rt_shift`.
3. Then Priority-M knobs conditional on best H-knob settings.
4. Priority-L knobs only if a specific pathology (unstable training, poor recovery) suggests it.
5. Bayesian optimization (e.g., `mlrMBO`, `ParBayesianOptimization`) is well-suited given the moderate-cost evaluations and high-dimensional space.

The alternative is a small 2-level factorial on the H knobs (~2^7 = 128 configs) — still tractable if each evaluation is under 15 min.
