# Pre-train the encoder via online denoising autoencoder

Generates fresh synthetic batches on-the-fly during training, giving
constant memory usage regardless of total samples seen. Recommended for
scaling to hundreds of thousands or millions of total samples.

## Usage

``` r
pretrain_denoising_online(
  peak_params = NULL,
  H,
  W,
  steps_per_epoch = 500L,
  epochs = 30L,
  batch_size = 32L,
  lr = 0.001,
  weight_decay = 1e-04,
  add_noise = TRUE,
  size_jitter = 0.15,
  grad_clip = 1,
  norm_clamp = 10,
  dust_threshold = 0.01,
  noisy_dust_threshold = NULL,
  location_mode = c("empirical", "marginal"),
  location_jitter_rt = 2,
  location_jitter_cv = 1,
  attribute_mode = c("joint", "marginal"),
  noise_scale = 1,
  catalog = NULL,
  anchor_ids = NULL,
  variable_config = NULL,
  contamination_config = NULL,
  clean_layers = NULL,
  noise_layers = NULL,
  contamination_schedule = NULL,
  mask_config = NULL,
  affine_shift_rt = 0L,
  affine_shift_cv = 0L,
  shift_group_size = 1L,
  baseline_drift_config = NULL,
  val_real = NULL,
  val_n = 200L,
  checkpoint_every = 10L,
  loss_diagnostics = TRUE,
  stem_stride_rt = 4L,
  stem_stride_cv = 1L,
  encoder_channels = c(16L, 16L, 32L, 64L, 64L),
  encoder_dropout = 0,
  anchor_dropout = 0,
  anchor_draw = c("all", "prevalence"),
  anchor_jitter_variability = NULL,
  affine_shift_prob = 1,
  stripe_noise_config = NULL,
  rt_warp_config = NULL,
  variable_source_peak_params = NULL,
  peak_shape = c("gaussian", "emg", "mix"),
  emg_tau_range = c(0.2, 0.8),
  emg_mix_prob = 0.3,
  presence_head = FALSE,
  presence_weight = 1,
  presence_compounds = NULL,
  presence_pos_weight_max = 20,
  presence_d_attn = 32L,
  device = if (torch::cuda_is_available()) "cuda" else "cpu",
  num_threads = NULL,
  num_workers = 1L,
  save_path = NULL,
  resume_from = NULL,
  seed = 42L,
  verbose = TRUE
)
```

## Arguments

- peak_params:

  Output of
  [`estimate_peak_params()`](https://mljaniczek.github.io/classydms/reference/estimate_peak_params.md).

- H, W:

  Synthetic image dimensions.

- steps_per_epoch:

  Number of training batches per epoch.

- epochs:

  Number of epochs.

- batch_size:

  Batch size.

- lr, weight_decay:

  Optimizer parameters.

- add_noise:

  If TRUE, denoising AE; if FALSE, reconstruction AE.

- size_jitter:

  Per-peak random log-normal size scale factor SD.

- grad_clip:

  Maximum L2 norm for gradient clipping (0 disables).

- norm_clamp:

  Maximum normalized pixel value (safety clamp).

- dust_threshold:

  Zero out pixels below this value in the CLEAN target BEFORE
  normalization. Matches the dust thresholding applied to real data via
  [`baseline_basement()`](https://mljaniczek.github.io/classydms/reference/baseline_basement.md)
  so synthetic and real samples have the same sparsity profile entering
  the encoder. **Default `0.01`** — was `0` previously, which produced
  dense synthetic Z matrices that mismatched real (baseline_basement'ed)
  validation samples and pushed the encoder toward feature collapse. Set
  to `0` to disable. If your real preprocessing uses a different
  `basement_thr`, match it here.

- noisy_dust_threshold:

  Zero-out threshold for the NOISY input. `NULL` (default) chooses a
  mode-appropriate value: in split-layer mode (any `clean_layers` /
  `noise_layers` set) the noisy input is left UNTHRESHOLDED so
  contamination and sensor noise remain full-resolution corruption for
  the encoder to strip; in standard mode the same `dust_threshold` is
  applied to noisy as to clean (matches previous behavior). Pass an
  explicit scalar to override either default (e.g., `0` to force no
  thresholding, or a small positive value to keep some background floor
  in the noisy input).

- location_mode:

  Passed to
  [`generate_one_synthetic()`](https://mljaniczek.github.io/classydms/reference/generate_one_synthetic.md).
  Default `"empirical"` samples peak locations from the real cohort's
  observed hotspots (`params$rt_loc_raw` / `cv_loc_raw`), producing
  synthetic samples whose peaks cluster at cohort-typical (RT, CV)
  coordinates rather than scattering across the image. `"marginal"` is
  the previous behavior kept for backward compatibility.

- location_jitter_rt, location_jitter_cv:

  Passed to
  [`generate_one_synthetic()`](https://mljaniczek.github.io/classydms/reference/generate_one_synthetic.md).
  SD of per-peak location jitter.

- attribute_mode:

  Passed to
  [`generate_one_synthetic()`](https://mljaniczek.github.io/classydms/reference/generate_one_synthetic.md).
  Default `"joint"` couples per-peak sigma_rt / sigma_cv / intensity by
  resampling the observed triple from `params$*_raw`, so wider peaks are
  drawn with correspondingly higher intensity (matching the physical
  coupling in real GC-DMS). `"marginal"` reverts to independent
  log-normal draws.

- noise_scale:

  Passed to
  [`generate_one_synthetic()`](https://mljaniczek.github.io/classydms/reference/generate_one_synthetic.md).
  Multiplier on the noise SD used at synthesis time. `1.0` (default)
  uses real-cohort noise levels; larger values increase the corruption
  the denoising autoencoder must learn to remove, which typically
  improves encoder feature quality up to some optimum (Vincent et al.
  2008). Sweep `c(1, 3, 5, 10)` and pick the largest value that doesn't
  degrade real-validation MSE. Accepts a **length-2 vector**
  `c(min, max)` to draw a fresh noise scale per sample uniformly from
  that range — teaches the encoder to handle varying real-instrument
  noise levels rather than assuming a fixed floor. Recommended range for
  the "wild" recipe: `c(1.0, 3.0)`.

- anchor_ids:

  Optional integer vector of catalog `compound_id`s. When provided
  (along with `catalog`), synthetic samples are generated via the
  three-layer noisy scheme (see
  [`generate_noisy_pretrain_sample()`](https://mljaniczek.github.io/classydms/reference/generate_noisy_pretrain_sample.md))
  rather than the default catalog-prevalence firing. The anchor set
  fires unconditionally in every sample; `variable_config` and
  `contamination_config` optionally add per-sample layers on top.

- variable_config:

  Optional list; see
  [`generate_noisy_pretrain_sample()`](https://mljaniczek.github.io/classydms/reference/generate_noisy_pretrain_sample.md).
  Requires `anchor_ids`. Adds a per-sample random subset of
  low-prevalence catalog compounds. `variable_config$draw` picks how:
  `"uniform"` (default) draws `n_per_sample` compounds uniformly from
  the eligible pool; `"prevalence"` fires each eligible compound with
  probability `min(1, prevalence * prevalence_scale)`
  (`prevalence_scale` default `1`), with `n_per_sample` then an optional
  cap (default none). The resolved mode and scale are checked at entry
  and recorded in the manifest as `hyperparams$variable_draw` /
  `hyperparams$variable_prevalence_scale`.

- contamination_config:

  Optional list; see
  [`generate_noisy_pretrain_sample()`](https://mljaniczek.github.io/classydms/reference/generate_noisy_pretrain_sample.md).
  Requires `anchor_ids`. Adds per-sample uniform-random peaks at
  positions outside catalog compound exclusion boxes. The default
  `intensity_range` and exclusion boxes are much brighter / wider than
  real off-catalog peaks; see that function's docs for measured values.

- clean_layers, noise_layers:

  Optional character vectors naming which of the three catalog layers
  (`"anchors"`, `"variable"`, `"contamination"`) go into the CLEAN
  target vs. only into the NOISY input. Default `NULL / NULL` uses the
  returned `pair$clean` and `pair$noisy` directly, which is the current
  behavior where all three layers are in the clean target and only
  sensor noise (from `noise_scale`) is corruption. Setting either
  switches to split-layer denoising: `Z_clean = sum of clean_layers` and
  `Z_noisy = Z_clean + sum of noise_layers + sensor_noise`. Recommended
  recipe for skin GC-DMS pretraining:
  `clean_layers = c("anchors", "variable")`,
  `noise_layers = c("contamination")`, plus a stronger
  `contamination_config` (say 30–100 peaks, intensity 0.3–1.0) and
  `noise_scale` bumped ~2×. This teaches the encoder to strip apparatus
  contamination and unfamiliar peaks while preserving the catalog
  compound structure — the actual denoising objective the training loss
  is measuring. Layers can't appear in both vectors; naming a layer
  whose config is `NULL` produces a warning.

- contamination_schedule:

  Optional named list of regime weights giving stochastic per-sample
  variation in contamination level. Currently uses fixed scale
  multipliers: `clean = 0x`, `light = 0.3x`, `heavy = 1x` — where the
  multiplier scales BOTH `contamination_config$n_per_sample` and
  `intensity_range`. Weights normalize automatically. Example:
  `contamination_schedule = list(clean = 0.3, light = 0.4, heavy = 0.3)`
  makes 30% of samples clean, 40% lightly contaminated, 30% heavily
  contaminated. Prevents the encoder from overfitting to "always messy"
  — real skin GC-DMS ranges from clean to messy per subject, and the
  encoder should generalize to both extremes. Requires
  `contamination_config` (the config to scale) and `anchor_ids` (the
  three-layer path). `NULL` (default) disables scheduling — every sample
  uses the full `contamination_config`.

- mask_config:

  Optional list enabling MAE-style masking augmentation. When provided,
  zeros out random rectangles from the NOISY input only (target
  unchanged), forcing the encoder to reconstruct the masked regions from
  surrounding context. Fields: `n_rects = c(min, max)` (default
  `c(3L, 5L)`) — how many rectangles per sample, drawn uniformly;
  `size_frac = c(min, max)` (default `c(0.05, 0.15)`) — each rectangle's
  area as a fraction of image size, drawn uniformly per rectangle.
  Aspect ratio is randomized per rectangle in `[0.5, 2.0]`. `NULL`
  (default) disables masking. Regularizes against
  synthetic-pixel-pattern memorization; effective on top of split-layer
  denoising.

- affine_shift_rt, affine_shift_cv:

  Coherent per-sample affine shift augmentation, in pixels. Each accepts
  either shape:

  - Scalar `M` (integer, default `0L`): fixed maximum magnitude.
    Per-sample shift `dr ~ U(-M, +M)` (or `dc` for CV).

  - Length-2 vector `c(M_lo, M_hi)`: stochastic magnitude. Each sample
    first draws a per-sample max `M ~ U(M_lo, M_hi)`, then
    `dr ~ U(-M, +M)`. Models real cohorts where drift magnitude itself
    varies session-to-session.

  Applied to BOTH clean and noisy identically (data augmentation, not
  drift correction). Reasonable values on 1400 x 90 images:
  `affine_shift_rt = 5L, affine_shift_cv = 1L` (fixed) or
  `affine_shift_rt = c(5L, 25L), affine_shift_cv = c(1L, 5L)`
  (stochastic).

- shift_group_size:

  Group size for the coherent shift draw. Accepts either shape:

  - Scalar `N` (integer, default `1L`): fixed. The first N samples in
    each batch share ONE shift draw (magnitude, direction, and Bernoulli
    activation all coherent within the group), the next N share another,
    etc. `1L` reproduces the current per-sample-independent behavior.

  - Length-2 vector `c(N_lo, N_hi)`: each successive group draws its own
    size from `Uniform(N_lo, N_hi)`, so group boundaries themselves vary
    within a batch. Matches real acquisitions where sessions differ in
    length.

  Setting scalar to `batch_size` models each batch as one "acquisition
  session" in which every sample drifted the same way. Also gates
  whether stripe patterns are grouped, when
  `stripe_noise_config$grouped = TRUE`.

- baseline_drift_config:

  Optional list controlling smooth low-frequency baseline drift added to
  the NOISY input only (clean target untouched — this is corruption the
  encoder must strip, teaching it to ignore residual baseline that AsLS
  didn't fully remove). Fields:

  - `prob` (default `0.3`) — Bernoulli fraction of samples that get
    baseline drift.

  - `direction` (`"rt"`, `"cv"`, or `"both"`, default `"rt"`) — axis of
    variation. `"both"` adds an outer sum of two 1D drifts. May be a
    subset e.g. `c("rt", "cv")` — one sampled per sample.

  - `shape` (`"polynomial"` or `"sinusoid"`, default `"polynomial"`) —
    smooth 1D function family.

  - `order` (default `c(2L, 4L)`) — polynomial-order range; per sample
    an integer is drawn from this range (ignored for sinusoid).

  - `amplitude` (default `c(0.02, 0.10)`) — peak-to-peak amplitude range
    as a fraction of image intensity units; drawn per sample.

- val_real:

  Optional list of real preprocessed, padded Z matrices (each of size
  `H x W`) held out as a real-data validation set. When provided, an
  additional per-epoch metric `real_val_MSE` is computed as the
  reconstruction MSE of the autoencoder on these real samples (target =
  input, no noise added). Directly measures whether encoder features
  learned from synthetic transfer to real. If `NULL` (default), only the
  synthetic validation metric is computed.

- val_n:

  Number of fixed validation samples (0 disables validation).

- checkpoint_every:

  Save encoder checkpoint every N epochs.

- loss_diagnostics:

  If TRUE, also print median, p95, max batch loss.

- stem_stride_rt, stem_stride_cv:

  Encoder stem strides.

- encoder_channels:

  Length-5 integer vector controlling encoder capacity:
  `c(stem_ch, layer1_ch, layer2_ch, layer3_ch, latent_ch)`. The final
  entry is the GAP-pooled feature dimension (input width for the
  downstream elastic-net classifier). Default
  `c(16L, 16L, 32L, 64L, 64L)` matches previous behavior exactly (no
  projection layer added). Wider examples:
  `c(16L, 32L, 64L, 128L, 128L)` for 128-dim latent;
  `c(16L, 32L, 64L, 128L, 256L)` for 256-dim latent (adds a 1x1
  projection on top). Same 32x/4x RT/CV spatial downsample regardless —
  only channel counts differ. Downstream classifier and feature
  extractor auto-detect the latent dim from the trained encoder, so no
  other args need updating.

- encoder_dropout:

  Channel-dropout probability for the encoder (see
  [`dms_encoder()`](https://mljaniczek.github.io/classydms/reference/dms_encoder.md)).
  Applied via `nn_dropout2d` between residual block groups — whole
  feature-map channels are zeroed independently with probability
  `encoder_dropout`. Standard regularization to reduce overfit to the
  synthetic distribution and improve transferability of encoder features
  to real data. Default `0.0` disables dropout (backward-compatible).
  Recommended try-values `0.1` and `0.2` when real-validation MSE
  plateaus or begins rising past a certain epoch. Decoder is left clean
  — only encoder-side matters for downstream transfer. Dropout is
  training-only; encoder\$eval() disables it for feature extraction and
  inference automatically.

- anchor_dropout:

  Per-sample Bernoulli drop probability applied to the anchor set before
  rendering. Each sample, each anchor compound is independently dropped
  with this probability; surviving anchors fire as usual. Default `0.0`
  — all anchors always fire (backward-compatible). Recommended
  `0.1`–`0.2` to simulate real subject-to-subject compound absence.

- anchor_draw:

  `"all"` (default) or `"prevalence"`; see
  [`generate_noisy_pretrain_sample()`](https://mljaniczek.github.io/classydms/reference/generate_noisy_pretrain_sample.md).
  With `"prevalence"`, each anchor fires with probability
  `prevalence * (1 - anchor_dropout)`. Recorded in the manifest as
  `hyperparams$anchor_draw`.

- anchor_jitter_variability:

  Optional length-2 numeric range for a per-sample scalar multiplier
  applied to `location_jitter_rt` and `location_jitter_cv`. When
  provided, each sample draws `jm ~ Uniform(low, high)` and effective
  jitter for that sample is `location_jitter_rt * jm` on RT and
  `location_jitter_cv * jm` on CV. Whole-sample scalar (per-anchor
  variation deferred as a follow-up refactor). Default `NULL` disables
  (fixed jitter per sample).

- affine_shift_prob:

  Bernoulli probability that any given sample (or group; see
  `shift_group_size`) gets an affine shift. Default `1.0` — every sample
  when `affine_shift_rt`/`cv` \> 0. Set to e.g. `0.6` to have 40% of
  samples receive NO shift (better mimics real acquisition where most
  runs are well-aligned and only some drift).

- stripe_noise_config:

  Optional list controlling row- or column-correlated stripe noise added
  to the NOISY input on a fraction of samples. Simulates thermal /
  electrical acquisition artifacts. Fields:

  - `prob` (default `0.2`) — Bernoulli fraction of samples that get
    stripe noise this batch.

  - `direction` (`"rt"` or `"cv"`, default `"rt"`) — direction along
    which the noise is 1-D correlated. `"rt"` gives horizontal bands
    (one noise value per row, tiled across cols); `"cv"` gives vertical
    bands.

  - `intensity` (default `c(0.05, 0.15)`) — length-2 SD range for the
    1-D noise. Drawn uniformly per sample.

  - `grouped` (default `FALSE`) — when `TRUE`, the stripe amplitude AND
    pattern are drawn once per shift group (see `shift_group_size`) and
    shared across all samples in that group, so a whole "acquisition
    session" gets the same stripe. Applied to noisy only (before dust
    threshold).

- rt_warp_config:

  Optional list controlling smooth per-sample RT warp applied
  identically to clean and noisy. Simulates column-condition drift where
  later peaks shift more than early peaks. Fields:

  - `strength` (default `0.02`) — fractional shift at the LAST row
    (`s * H` pixels). Linear ramp from 0 at row 1. Sign is randomized
    per sample. `|strength|` must be `< 0.5`.

  - `prob` (default `1.0`) — Bernoulli fraction of samples that get
    warped. Nearest-neighbor row indexing (O(H) per sample).

- variable_source_peak_params:

  Optional `peak_params` bundle (output of
  [`estimate_peak_params()`](https://mljaniczek.github.io/classydms/reference/estimate_peak_params.md))
  used to render the VARIABLE layer via
  [`generate_one_synthetic()`](https://mljaniczek.github.io/classydms/reference/generate_one_synthetic.md)
  instead of the default catalog-clustered path. When provided, each
  sample's variable layer is generated by bootstrapping real observed
  peaks from this bundle — real per-peak morphology, empirical location
  distribution, coupled `(sigma_rt, sigma_cv, intensity)` triples. Peak
  count is driven by `variable_source_peak_params$n_peaks$mean/sd` (each
  sample draws a count from that Normal). Bypasses `variable_config`
  entirely; passing both errors. Recommended when you want the variable
  layer to reflect real skin-VOC morphology diversity that clustering
  averages out. `NULL` (default) uses the catalog-based path controlled
  by `variable_config`.

- peak_shape:

  One of `"gaussian"` (default), `"emg"`, or `"mix"`. Selects the
  per-peak shape family:

  - `"gaussian"` — symmetric Gaussian on both axes. Preserves previous
    behavior exactly.

  - `"emg"` — bi-Gaussian along RT (symmetric on the left of the peak
    center, tail-stretched on the right). CV stays symmetric.
    Approximates the asymmetric shape of real GC peaks (sharp leading
    edge, gentle trailing edge). Tau (right-side widening factor) drawn
    per peak from `emg_tau_range`.

  - `"mix"` — per peak, Bernoulli(`emg_mix_prob`) picks EMG; otherwise
    Gaussian. Gives a natural mix of symmetric and asymmetric peaks per
    sample. Recommended when transitioning from pure-Gaussian to EMG
    since real cohorts have both.

- emg_tau_range:

  Length-2 numeric range for the EMG tau parameter (right-side sigma
  widening factor). Default `c(0.2, 0.8)` corresponds to right side
  20-80% wider than left. Drawn per peak uniformly.

- emg_mix_prob:

  Probability of EMG shape per peak in `"mix"` mode. Default `0.3`.

- presence_head:

  If `TRUE`, trains a
  [`dms_presence_head()`](https://mljaniczek.github.io/classydms/reference/dms_presence_head.md)
  jointly with the autoencoder: from the encoder feature map it predicts
  which catalog compounds were rendered in the clean target layers
  (ground truth comes straight from the generator). Total loss is
  `reconstruction MSE + presence_weight * BCE`. Requires the three-layer
  path (`catalog` + `anchor_ids`). The head is saved as
  `*_presence_head.pt` next to the encoder and returned as
  `$presence_head`; use it downstream via
  [`extract_encoder_features()`](https://mljaniczek.github.io/classydms/reference/extract_encoder_features.md)
  /
  [`cv_encoder_elastic_net()`](https://mljaniczek.github.io/classydms/reference/cv_encoder_elastic_net.md)
  with `readout = "presence"`, or
  [`build_pretrained_classifier()`](https://mljaniczek.github.io/classydms/reference/build_pretrained_classifier.md)
  with `readout = "presence"`. Default `FALSE`.

- presence_weight:

  Weight on the presence BCE term. Default `1.0`.

- presence_compounds:

  Optional subset of `compound_id`s the head covers. Default `NULL` =
  every catalog compound with a finite position. Compounds that never
  fire during pretraining (e.g. outside
  `variable_config$prevalence_range` and not anchors) only ever see
  negative labels.

- presence_pos_weight_max:

  Cap on the per-compound positive weight `(1 - f) / f`, where `f` is
  the compound's running firing frequency. Rarity weighting keeps
  rarely-rendered compounds from being ignored by the loss. Default
  `20`.

- presence_d_attn:

  Query/key dimension of the head's attention. Default `32L`.

- device:

  "cpu" or "cuda".

- num_threads:

  If non-NULL, sets the number of CPU threads used by torch for
  forward/backward passes via
  [`torch::torch_set_num_threads()`](https://torch.mlverse.org/docs/reference/threads.html).
  Apple Silicon machines benefit most from setting this to the number of
  performance cores (e.g. 4 on an M-series Pro chip with 4P+10E layout),
  since the efficiency cores contribute less per thread. Default `NULL`
  leaves torch's automatic choice in place.

- num_workers:

  Number of parallel workers for R-side synthetic data generation
  (default `1L` = serial). When \> 1, uses
  [`future.apply::future_lapply`](https://future.apply.futureverse.org/reference/future_lapply.html)
  with a
  [`future::multisession`](https://future.futureverse.org/reference/multisession.html)
  plan — fresh R processes (not fork), so it side-steps the macOS
  libtorch fork crash where `mclapply` children inherit poisoned
  threadpool mutexes. Cross-platform (Linux, macOS, Windows) and handles
  global-variable discovery automatically. Most useful on GPU / MPS
  backends where the GPU is otherwise waiting on R between batches. A
  reasonable choice on a chip with N performance cores is
  `num_workers = N - 2` (leaving a couple cores for the GPU's
  communication threads, BLAS, and the OS). Each worker task gets an
  independent L'Ecuyer-CMRG substream via `future.seed = TRUE` for
  reproducibility. Requires the `future` and `future.apply` packages
  when \> 1; falls back to serial with a message if either is missing.

- save_path:

  If non-NULL, writes encoder, autoencoder, and manifest to disk after
  every epoch (so a mid-run crash leaves recoverable artifacts).

- resume_from:

  If non-NULL, path to a previously-saved encoder (`*.pt`). The function
  looks for the corresponding `_autoencoder.pt` and `_manifest.Rdata`
  next to it and resumes training from `last_epoch_completed + 1`,
  inheriting prior loss history. Useful for a "trained 32 epochs, want
  64" continuation and for picking up after a crash. All
  training-affecting args (`H`, `W`, `batch_size`, `lr`, `weight_decay`,
  `noise_scale`, `size_jitter`, `dust_threshold`, `location_mode`,
  `location_jitter_rt/cv`, `attribute_mode`, `stem_stride_rt/cv`,
  `anchor_ids`, `variable_config`, `contamination_config`,
  `clean_layers`, `noise_layers`, `mask_config`, `affine_shift_rt/cv`,
  `shift_group_size`, `baseline_drift_config`) should match the original
  run. Values that differ trigger a warning on load listing the
  mismatches — training continues with the CURRENT-call values, but the
  model was trained so far under the manifest values. To resume
  identically, read the manifest and copy its `hyperparams` list:

      e <- new.env(); load("path/to/enc_manifest.Rdata", envir = e)
      hp <- e$training_manifest$hyperparams
      hp$epochs <- 64L                     # extend past original target
      hp$noise_scale_range <- NULL         # internal-only field, drop it
      do.call(pretrain_denoising_online,
        c(list(catalog = catalog,
                resume_from = "path/to/enc.pt",
                save_path   = "path/to/enc.pt"), hp))

  Two housekeeping notes on the replay pattern:

  1.  `noise_scale_range` is an internal-derived field on the manifest;
      drop it before splicing (the original `noise_scale` still lives in
      the same hp list, so scalar/range detection re-runs cleanly).

  2.  `epochs` and `save_path` are the two things you typically override
      at resume time — the former to extend training past the original
      target, the latter to overwrite the same encoder file.

- seed:

  RNG seed.

- verbose:

  Whether to print epoch progress.

## Value

List with `encoder`, `autoencoder`, `presence_head` (`NULL` unless
`presence_head = TRUE`), `loss_history` (reconstruction MSE),
`val_loss_history`, `real_val_loss_history`, `presence_loss_history`,
`val_presence_auc_history` (macro per-compound AUC on the synthetic val
set), `batch_loss_stats`, and `total_samples`.

## Details

Includes scale-up safety features: a fixed validation set held constant
across epochs for stable progress tracking, periodic encoder
checkpointing for crash recovery, per-epoch manifest save with all
hyperparameters for reproducibility, and per-epoch ETA estimate.

Total samples seen = `steps_per_epoch * epochs * batch_size`.
