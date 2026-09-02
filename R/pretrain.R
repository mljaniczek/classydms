# classydms: self-supervised pretraining of the encoder on synthetic data

#' Pre-train the encoder via batch denoising autoencoder
#'
#' Generates `n_synthetic` synthetic GC-DMS images, trains a denoising
#' autoencoder on them (noisy -> clean), and returns the trained encoder.
#' Memory-limited: the full synthetic dataset is held in RAM. For larger
#' totals, use [pretrain_denoising_online()].
#'
#' @param peak_params Output of [estimate_peak_params()].
#' @param H,W Synthetic image dimensions (must match downstream input).
#' @param n_synthetic Number of synthetic images.
#' @param epochs,batch_size,lr,weight_decay Optimization parameters.
#' @param add_noise If TRUE, denoising AE; if FALSE, reconstruction AE.
#' @param stem_stride_rt,stem_stride_cv Stem stride for the encoder.
#' @param device "cpu" or "cuda".
#' @param save_path If non-NULL, writes encoder and full autoencoder to disk.
#' @param seed RNG seed.
#' @return List with `encoder`, `autoencoder`, and `loss_history`.
#' @seealso [pretrain_denoising_online()]
#' @export
pretrain_autoencoder <- function(peak_params, H, W,
                                  n_synthetic = 5000L,
                                  epochs = 50L,
                                  batch_size = 32L,
                                  lr = 1e-3,
                                  weight_decay = 1e-4,
                                  add_noise = TRUE,
                                  stem_stride_rt = 4L,
                                  stem_stride_cv = 1L,
                                  device = if (torch::cuda_is_available()) "cuda" else "cpu",
                                  save_path = NULL,
                                  seed = 42L) {
  set.seed(seed); torch::torch_manual_seed(seed)

  message("Generating ", n_synthetic, " synthetic samples (", H, "x", W, ")...")
  message("Mode: ", if (add_noise) "denoising (noisy -> clean)"
          else "reconstruction (clean -> clean)")
  synth <- generate_synthetic_dataset(peak_params, N = n_synthetic, H = H,
                                       W = W, add_noise = add_noise,
                                       normalize = TRUE)
  message("Done. Building autoencoder...")

  ds_class <- torch::dataset(
    name = "denoising_dataset",
    initialize = function() {},
    .getitem = function(i) list(input = synth$noisy[i, ..],
                                target = synth$clean[i, ..]),
    .length = function() synth$noisy$size()[1]
  )
  ds <- ds_class()
  dl <- torch::dataloader(ds, batch_size = batch_size, shuffle = TRUE)

  model <- dms_denoising_autoencoder(target_H = H, target_W = W,
                                      stem_stride_rt = stem_stride_rt,
                                      stem_stride_cv = stem_stride_cv)
  model$to(device = device)
  opt <- torch::optim_adam(model$parameters, lr = lr, weight_decay = weight_decay)
  loss_fn <- torch::nn_mse_loss()

  loss_history <- numeric(epochs)
  for (ep in seq_len(epochs)) {
    model$train(); epoch_loss <- 0; n_batches <- 0
    coro::loop(for (batch in dl) {
      input <- batch$input$to(device = device)
      target <- batch$target$to(device = device)
      opt$zero_grad()
      recon <- model(input)
      loss <- loss_fn(recon, target)
      loss$backward(); opt$step()
      epoch_loss <- epoch_loss + as.numeric(loss$item())
      n_batches <- n_batches + 1
    })
    loss_history[ep] <- epoch_loss / n_batches
    if (ep %% 10 == 0 || ep == 1) {
      message("Epoch ", ep, "/", epochs, " | MSE: ",
              round(loss_history[ep], 6))
    }
  }

  if (!is.null(save_path)) {
    torch::torch_save(model$encoder, save_path)
    message("Encoder saved to: ", save_path)
    ae_path <- sub("\\.pt$", "_autoencoder.pt", save_path)
    if (ae_path == save_path) ae_path <- paste0(save_path, "_autoencoder.pt")
    torch::torch_save(model, ae_path)
    message("Autoencoder saved to: ", ae_path)
  }
  list(encoder = model$encoder, autoencoder = model, loss_history = loss_history)
}

#' @rdname pretrain_autoencoder
#' @export
pretrain_denoising <- pretrain_autoencoder

#' Pre-train the encoder from a peak catalog (smoke test only)
#'
#' Same idea as [pretrain_autoencoder()] but generates the synthetic
#' pretraining set via [generate_one_synthetic_from_catalog()]. Use
#' this when you have a catalog (perhaps loaded from disk) but not the
#' original `peak_params` object.
#'
#' **For real pretraining runs prefer
#' [pretrain_denoising_online()] with `catalog = your_catalog`** — that
#' generates fresh batches per training step, so the model can see
#' effectively unlimited unique synthetic samples with constant memory.
#' This pre-alloc version holds the full `n_synthetic` × H × W tensor
#' in RAM and reuses the same set every epoch, which caps the effective
#' training data. It stays available for smoke tests
#' (`n_synthetic <= a few thousand`).
#'
#' @param catalog A `peak_catalog` object.
#' @param H,W Synthetic image dimensions.
#' @param n_synthetic Number of synthetic samples to generate.
#' @param epochs,batch_size,lr,weight_decay Optimization parameters.
#' @param add_noise If TRUE, denoising autoencoder; if FALSE,
#'   reconstruction.
#' @param stem_stride_rt,stem_stride_cv Encoder stem strides.
#' @param device "cpu" or "cuda".
#' @param seed RNG seed.
#' @return List with `encoder`, `autoencoder`, `loss_history`.
#' @export
pretrain_autoencoder_from_catalog <- function(catalog, H, W,
                                                n_synthetic = 2000L,
                                                epochs = 30L,
                                                batch_size = 32L,
                                                lr = 1e-3,
                                                weight_decay = 1e-4,
                                                add_noise = TRUE,
                                                stem_stride_rt = 4L,
                                                stem_stride_cv = 1L,
                                                device = if (torch::cuda_is_available()) "cuda" else "cpu",
                                                seed = 42L) {
  set.seed(seed); torch::torch_manual_seed(seed)
  stopifnot(inherits(catalog, "peak_catalog"))
  message("Generating ", n_synthetic, " catalog-based synthetic samples ",
          "(", H, "x", W, ")...")

  # Generate stack of clean + noisy
  clean_arr <- array(0, dim = c(n_synthetic, 1L, H, W))
  noisy_arr <- array(0, dim = c(n_synthetic, 1L, H, W))
  for (i in seq_len(n_synthetic)) {
    res <- generate_one_synthetic_from_catalog(catalog, H = H, W = W,
                                                 add_noise = add_noise)
    clean_arr[i, 1L, , ] <- res$clean
    noisy_arr[i, 1L, , ] <- if (add_noise) res$noisy else res$clean
    if (i %% 200L == 0L) message("  generated ", i, "/", n_synthetic)
  }
  # Log-quantile normalize using the same convention as the peak_params path.
  # pmax(0, x) collapses dims when one arg is scalar, so clamp in-place instead.
  q_clean <- as.numeric(stats::quantile(clean_arr[clean_arr > 0], 0.95,
                                          na.rm = TRUE))
  if (!is.finite(q_clean) || q_clean == 0) q_clean <- 1
  clean_arr[clean_arr < 0] <- 0
  noisy_arr[noisy_arr < 0] <- 0
  clean_arr <- log1p(clean_arr) / log1p(q_clean)
  noisy_arr <- log1p(noisy_arr) / log1p(q_clean)

  clean_t <- torch::torch_tensor(clean_arr, dtype = torch::torch_float())
  noisy_t <- torch::torch_tensor(noisy_arr, dtype = torch::torch_float())

  ds_class <- torch::dataset(
    name = "catalog_pretrain_ds",
    initialize = function() {},
    .getitem = function(i) list(input = noisy_t[i, ..],
                                 target = clean_t[i, ..]),
    .length = function() clean_t$size()[1]
  )
  dl <- torch::dataloader(ds_class(), batch_size = batch_size,
                            shuffle = TRUE)

  model <- dms_denoising_autoencoder(target_H = H, target_W = W,
                                       stem_stride_rt = stem_stride_rt,
                                       stem_stride_cv = stem_stride_cv)
  model$to(device = device)
  opt <- torch::optim_adam(model$parameters, lr = lr,
                            weight_decay = weight_decay)
  loss_fn <- torch::nn_mse_loss()
  loss_history <- numeric(epochs)

  for (ep in seq_len(epochs)) {
    model$train(); epoch_loss <- 0; n_batches <- 0
    coro::loop(for (batch in dl) {
      input <- batch$input$to(device = device)
      target <- batch$target$to(device = device)
      opt$zero_grad()
      recon <- model(input)
      loss <- loss_fn(recon, target)
      loss$backward(); opt$step()
      epoch_loss <- epoch_loss + as.numeric(loss$item())
      n_batches <- n_batches + 1
    })
    loss_history[ep] <- epoch_loss / max(1, n_batches)
    if (ep %% 10 == 0 || ep == 1) {
      message("Epoch ", ep, "/", epochs, " | MSE: ",
              round(loss_history[ep], 6))
    }
  }
  list(encoder = model$encoder, autoencoder = model,
       loss_history = loss_history)
}

#' Pre-train the encoder via online denoising autoencoder
#'
#' Generates fresh synthetic batches on-the-fly during training, giving
#' constant memory usage regardless of total samples seen. Recommended
#' for scaling to hundreds of thousands or millions of total samples.
#'
#' Includes scale-up safety features: a fixed validation set held constant
#' across epochs for stable progress tracking, periodic encoder
#' checkpointing for crash recovery, per-epoch manifest save with all
#' hyperparameters for reproducibility, and per-epoch ETA estimate.
#'
#' Total samples seen = `steps_per_epoch * epochs * batch_size`.
#'
#' @param peak_params Output of [estimate_peak_params()].
#' @param H,W Synthetic image dimensions.
#' @param steps_per_epoch Number of training batches per epoch.
#' @param epochs Number of epochs.
#' @param batch_size Batch size.
#' @param lr,weight_decay Optimizer parameters.
#' @param add_noise If TRUE, denoising AE; if FALSE, reconstruction AE.
#' @param size_jitter Per-peak random log-normal size scale factor SD.
#' @param grad_clip Maximum L2 norm for gradient clipping (0 disables).
#' @param norm_clamp Maximum normalized pixel value (safety clamp).
#' @param dust_threshold Zero out pixels below this value in the CLEAN
#'   target BEFORE normalization. Matches the dust thresholding
#'   applied to real data via [baseline_basement()] so synthetic and
#'   real samples have the same sparsity profile entering the encoder.
#'   **Default `0.01`** — was `0` previously, which produced dense
#'   synthetic Z matrices that mismatched real (baseline_basement'ed)
#'   validation samples and pushed the encoder toward feature collapse.
#'   Set to `0` to disable. If your real preprocessing uses a
#'   different `basement_thr`, match it here.
#' @param noisy_dust_threshold Zero-out threshold for the NOISY input.
#'   `NULL` (default) chooses a mode-appropriate value: in split-layer
#'   mode (any `clean_layers` / `noise_layers` set) the noisy input
#'   is left UNTHRESHOLDED so contamination and sensor noise remain
#'   full-resolution corruption for the encoder to strip; in standard
#'   mode the same `dust_threshold` is applied to noisy as to clean
#'   (matches previous behavior). Pass an explicit scalar to override
#'   either default (e.g., `0` to force no thresholding, or a small
#'   positive value to keep some background floor in the noisy input).
#' @param location_mode Passed to [generate_one_synthetic()]. Default
#'   `"empirical"` samples peak locations from the real cohort's
#'   observed hotspots (`params$rt_loc_raw` / `cv_loc_raw`), producing
#'   synthetic samples whose peaks cluster at cohort-typical
#'   (RT, CV) coordinates rather than scattering across the image.
#'   `"marginal"` is the previous behavior kept for backward compatibility.
#' @param location_jitter_rt,location_jitter_cv Passed to
#'   [generate_one_synthetic()]. SD of per-peak location jitter.
#' @param attribute_mode Passed to [generate_one_synthetic()]. Default
#'   `"joint"` couples per-peak sigma_rt / sigma_cv / intensity by
#'   resampling the observed triple from `params$*_raw`, so wider
#'   peaks are drawn with correspondingly higher intensity (matching
#'   the physical coupling in real GC-DMS). `"marginal"` reverts to
#'   independent log-normal draws.
#' @param noise_scale Passed to [generate_one_synthetic()]. Multiplier
#'   on the noise SD used at synthesis time. `1.0` (default) uses
#'   real-cohort noise levels; larger values increase the corruption
#'   the denoising autoencoder must learn to remove, which typically
#'   improves encoder feature quality up to some optimum
#'   (Vincent et al. 2008). Sweep `c(1, 3, 5, 10)` and pick the
#'   largest value that doesn't degrade real-validation MSE. Accepts
#'   a **length-2 vector** `c(min, max)` to draw a fresh noise scale
#'   per sample uniformly from that range — teaches the encoder to
#'   handle varying real-instrument noise levels rather than assuming
#'   a fixed floor. Recommended range for the "wild" recipe:
#'   `c(1.0, 3.0)`.
#' @param anchor_ids Optional integer vector of catalog `compound_id`s.
#'   When provided (along with `catalog`), synthetic samples are
#'   generated via the three-layer noisy scheme (see
#'   [generate_noisy_pretrain_sample()]) rather than the default
#'   catalog-prevalence firing. The anchor set fires unconditionally
#'   in every sample; `variable_config` and `contamination_config`
#'   optionally add per-sample layers on top.
#' @param variable_config Optional list; see
#'   [generate_noisy_pretrain_sample()]. Requires `anchor_ids`. Adds
#'   a per-sample random subset of low-prevalence catalog compounds.
#' @param contamination_config Optional list; see
#'   [generate_noisy_pretrain_sample()]. Requires `anchor_ids`. Adds
#'   per-sample uniform-random peaks at positions outside catalog
#'   compound exclusion boxes.
#' @param clean_layers,noise_layers Optional character vectors naming
#'   which of the three catalog layers (`"anchors"`, `"variable"`,
#'   `"contamination"`) go into the CLEAN target vs. only into the
#'   NOISY input. Default `NULL / NULL` uses the returned
#'   `pair$clean` and `pair$noisy` directly, which is the current
#'   behavior where all three layers are in the clean target and only
#'   sensor noise (from `noise_scale`) is corruption. Setting either
#'   switches to split-layer denoising: `Z_clean = sum of clean_layers`
#'   and `Z_noisy = Z_clean + sum of noise_layers + sensor_noise`.
#'   Recommended recipe for skin GC-DMS pretraining:
#'   `clean_layers = c("anchors", "variable")`,
#'   `noise_layers = c("contamination")`, plus a stronger
#'   `contamination_config` (say 30–100 peaks, intensity 0.3–1.0) and
#'   `noise_scale` bumped ~2×. This teaches the encoder to strip
#'   apparatus contamination and unfamiliar peaks while preserving
#'   the catalog compound structure — the actual denoising objective
#'   the training loss is measuring. Layers can't appear in both
#'   vectors; naming a layer whose config is `NULL` produces a warning.
#' @param contamination_schedule Optional named list of regime weights
#'   giving stochastic per-sample variation in contamination level.
#'   Currently uses fixed scale multipliers: `clean = 0x`,
#'   `light = 0.3x`, `heavy = 1x` — where the multiplier scales BOTH
#'   `contamination_config$n_per_sample` and `intensity_range`. Weights
#'   normalize automatically. Example:
#'   `contamination_schedule = list(clean = 0.3, light = 0.4, heavy = 0.3)`
#'   makes 30% of samples clean, 40% lightly contaminated, 30% heavily
#'   contaminated. Prevents the encoder from overfitting to
#'   "always messy" — real skin GC-DMS ranges from clean to messy per
#'   subject, and the encoder should generalize to both extremes.
#'   Requires `contamination_config` (the config to scale) and
#'   `anchor_ids` (the three-layer path). `NULL` (default) disables
#'   scheduling — every sample uses the full `contamination_config`.
#' @param mask_config Optional list enabling MAE-style masking
#'   augmentation. When provided, zeros out random rectangles from
#'   the NOISY input only (target unchanged), forcing the encoder to
#'   reconstruct the masked regions from surrounding context. Fields:
#'   `n_rects = c(min, max)` (default `c(3L, 5L)`) — how many
#'   rectangles per sample, drawn uniformly; `size_frac = c(min, max)`
#'   (default `c(0.05, 0.15)`) — each rectangle's area as a fraction
#'   of image size, drawn uniformly per rectangle. Aspect ratio is
#'   randomized per rectangle in `[0.5, 2.0]`. `NULL` (default)
#'   disables masking. Regularizes against synthetic-pixel-pattern
#'   memorization; effective on top of split-layer denoising.
#' @param affine_shift_rt,affine_shift_cv Integer max magnitudes for
#'   coherent per-sample affine shift augmentation. When > 0, each
#'   sample's clean and noisy matrices are shifted by
#'   `dr ~ U(-affine_shift_rt, affine_shift_rt)` rows and
#'   `dc ~ U(-affine_shift_cv, affine_shift_cv)` cols, with zero-fill
#'   at exposed edges. Applied to BOTH clean and noisy identically —
#'   this is data augmentation, not drift correction. Teaches the
#'   encoder that a compound at row r is the same as one at row r+k.
#'   Defaults `0L / 0L` (no shift, backward-compatible). Reasonable
#'   values: `affine_shift_rt = 5L, affine_shift_cv = 1L` on
#'   1400 x 90 images.
#' @param affine_shift_prob Bernoulli probability that any given
#'   sample gets an affine shift. Default `1.0` — every sample when
#'   `affine_shift_rt`/`cv` > 0. Set to e.g. `0.6` to have 40% of
#'   samples receive NO shift (better mimics real acquisition where
#'   most runs are well-aligned and only some drift).
#' @param anchor_dropout Per-sample Bernoulli drop probability applied
#'   to the anchor set before rendering. Each sample, each anchor
#'   compound is independently dropped with this probability;
#'   surviving anchors fire as usual. Default `0.0` — all anchors
#'   always fire (backward-compatible). Recommended `0.1`–`0.2` to
#'   simulate real subject-to-subject compound absence.
#' @param anchor_jitter_variability Optional length-2 numeric range
#'   for a per-sample scalar multiplier applied to
#'   `location_jitter_rt` and `location_jitter_cv`. When provided,
#'   each sample draws `jm ~ Uniform(low, high)` and effective
#'   jitter for that sample is `location_jitter_rt * jm` on RT and
#'   `location_jitter_cv * jm` on CV. Whole-sample scalar
#'   (per-anchor variation deferred as a follow-up refactor).
#'   Default `NULL` disables (fixed jitter per sample).
#' @param stripe_noise_config Optional list controlling row- or
#'   column-correlated stripe noise added to the NOISY input on a
#'   fraction of samples. Simulates thermal / electrical acquisition
#'   artifacts. Fields:
#'   - `prob` (default `0.2`) — Bernoulli fraction of samples that
#'     get stripe noise this batch.
#'   - `direction` (`"rt"` or `"cv"`, default `"rt"`) — direction
#'     along which the noise is 1-D correlated. `"rt"` gives
#'     horizontal bands (one noise value per row, tiled across
#'     cols); `"cv"` gives vertical bands.
#'   - `intensity` (default `c(0.05, 0.15)`) — length-2 SD range
#'     for the 1-D noise. Drawn uniformly per sample.
#'   Applied to noisy only (before dust threshold).
#' @param rt_warp_config Optional list controlling smooth per-sample
#'   RT warp applied identically to clean and noisy. Simulates
#'   column-condition drift where later peaks shift more than early
#'   peaks. Fields:
#'   - `strength` (default `0.02`) — fractional shift at the LAST
#'     row (`s * H` pixels). Linear ramp from 0 at row 1. Sign is
#'     randomized per sample. `|strength|` must be `< 0.5`.
#'   - `prob` (default `1.0`) — Bernoulli fraction of samples that
#'     get warped.
#'   Nearest-neighbor row indexing (O(H) per sample).
#' @param val_real Optional list of real preprocessed, padded Z
#'   matrices (each of size `H x W`) held out as a real-data
#'   validation set. When provided, an additional per-epoch metric
#'   `real_val_MSE` is computed as the reconstruction MSE of the
#'   autoencoder on these real samples (target = input, no noise
#'   added). Directly measures whether encoder features learned from
#'   synthetic transfer to real. If `NULL` (default), only the
#'   synthetic validation metric is computed.
#' @param val_n Number of fixed validation samples (0 disables validation).
#' @param checkpoint_every Save encoder checkpoint every N epochs.
#' @param loss_diagnostics If TRUE, also print median, p95, max batch loss.
#' @param stem_stride_rt,stem_stride_cv Encoder stem strides.
#' @param encoder_channels Length-5 integer vector controlling encoder
#'   capacity: `c(stem_ch, layer1_ch, layer2_ch, layer3_ch, latent_ch)`.
#'   The final entry is the GAP-pooled feature dimension (input width
#'   for the downstream elastic-net classifier). Default
#'   `c(16L, 16L, 32L, 64L, 64L)` matches previous behavior exactly
#'   (no projection layer added). Wider examples:
#'   `c(16L, 32L, 64L, 128L, 128L)` for 128-dim latent;
#'   `c(16L, 32L, 64L, 128L, 256L)` for 256-dim latent (adds a 1x1
#'   projection on top). Same 32x/4x RT/CV spatial downsample
#'   regardless — only channel counts differ. Downstream classifier
#'   and feature extractor auto-detect the latent dim from the
#'   trained encoder, so no other args need updating.
#' @param encoder_dropout Channel-dropout probability for the encoder
#'   (see [dms_encoder()]). Applied via `nn_dropout2d` between
#'   residual block groups — whole feature-map channels are zeroed
#'   independently with probability `encoder_dropout`. Standard
#'   regularization to reduce overfit to the synthetic distribution
#'   and improve transferability of encoder features to real data.
#'   Default `0.0` disables dropout (backward-compatible).
#'   Recommended try-values `0.1` and `0.2` when real-validation MSE
#'   plateaus or begins rising past a certain epoch. Decoder is left
#'   clean — only encoder-side matters for downstream transfer.
#'   Dropout is training-only; encoder$eval() disables it for
#'   feature extraction and inference automatically.
#' @param device "cpu" or "cuda".
#' @param num_threads If non-NULL, sets the number of CPU threads
#'   used by torch for forward/backward passes via
#'   [torch::torch_set_num_threads()]. Apple Silicon machines benefit
#'   most from setting this to the number of performance cores (e.g. 4
#'   on an M-series Pro chip with 4P+10E layout), since the efficiency
#'   cores contribute less per thread. Default `NULL` leaves torch's
#'   automatic choice in place.
#' @param num_workers Number of parallel workers for R-side synthetic
#'   data generation (default `1L` = serial). When > 1, uses
#'   `future.apply::future_lapply` with a `future::multisession` plan
#'   — fresh R processes (not fork), so it side-steps the macOS
#'   libtorch fork crash where `mclapply` children inherit poisoned
#'   threadpool mutexes. Cross-platform (Linux, macOS, Windows) and
#'   handles global-variable discovery automatically. Most useful on
#'   GPU / MPS backends where the GPU is otherwise waiting on R
#'   between batches. A reasonable choice on a chip with N performance
#'   cores is `num_workers = N - 2` (leaving a couple cores for the
#'   GPU's communication threads, BLAS, and the OS). Each worker task
#'   gets an independent L'Ecuyer-CMRG substream via
#'   `future.seed = TRUE` for reproducibility. Requires the `future`
#'   and `future.apply` packages when > 1; falls back to serial with
#'   a message if either is missing.
#' @param save_path If non-NULL, writes encoder, autoencoder, and manifest
#'   to disk after every epoch (so a mid-run crash leaves recoverable
#'   artifacts).
#' @param resume_from If non-NULL, path to a previously-saved encoder
#'   (`*.pt`). The function looks for the corresponding `_autoencoder.pt`
#'   and `_manifest.Rdata` next to it and resumes training from
#'   `last_epoch_completed + 1`, inheriting prior loss history. Useful
#'   for a "trained 32 epochs, want 64" continuation and for picking up
#'   after a crash. All training-affecting args (`H`, `W`,
#'   `batch_size`, `lr`, `weight_decay`, `noise_scale`,
#'   `size_jitter`, `dust_threshold`, `location_mode`,
#'   `location_jitter_rt/cv`, `attribute_mode`, `stem_stride_rt/cv`,
#'   `anchor_ids`, `variable_config`, `contamination_config`,
#'   `clean_layers`, `noise_layers`, `mask_config`, `affine_shift_rt/cv`)
#'   should match the original run. Values that differ trigger a
#'   warning on load listing the mismatches — training continues with
#'   the CURRENT-call values, but the model was trained so far under
#'   the manifest values. To resume identically, read the manifest and
#'   copy its `hyperparams` list:
#'
#'   ```
#'   e <- new.env(); load("path/to/enc_manifest.Rdata", envir = e)
#'   hp <- e$training_manifest$hyperparams
#'   hp$epochs <- 64L                     # extend past original target
#'   hp$noise_scale_range <- NULL         # internal-only field, drop it
#'   do.call(pretrain_denoising_online,
#'     c(list(catalog = catalog,
#'             resume_from = "path/to/enc.pt",
#'             save_path   = "path/to/enc.pt"), hp))
#'   ```
#'
#'   Two housekeeping notes on the replay pattern:
#'   1. `noise_scale_range` is an internal-derived field on the manifest;
#'      drop it before splicing (the original `noise_scale` still lives
#'      in the same hp list, so scalar/range detection re-runs cleanly).
#'   2. `epochs` and `save_path` are the two things you typically override
#'      at resume time — the former to extend training past the original
#'      target, the latter to overwrite the same encoder file.
#' @param seed RNG seed.
#' @param verbose Whether to print epoch progress.
#' @return List with `encoder`, `autoencoder`, `loss_history`,
#'   `val_loss_history`, `batch_loss_stats`, and `total_samples`.
#' @export
pretrain_denoising_online <- function(peak_params = NULL, H, W,
                                       steps_per_epoch = 500L,
                                       epochs = 30L,
                                       batch_size = 32L,
                                       lr = 1e-3,
                                       weight_decay = 1e-4,
                                       add_noise = TRUE,
                                       size_jitter = 0.15,
                                       grad_clip = 1.0,
                                       norm_clamp = 10.0,
                                       dust_threshold = 0.01,
                                       noisy_dust_threshold = NULL,
                                       location_mode = c("empirical",
                                                          "marginal"),
                                       location_jitter_rt = 2,
                                       location_jitter_cv = 1,
                                       attribute_mode = c("joint",
                                                          "marginal"),
                                       noise_scale = 1.0,
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
                                       val_real = NULL,
                                       val_n = 200L,
                                       checkpoint_every = 10L,
                                       loss_diagnostics = TRUE,
                                       stem_stride_rt = 4L,
                                       stem_stride_cv = 1L,
                                       encoder_channels = c(16L, 16L, 32L, 64L, 64L),
                                       encoder_dropout = 0.0,
                                       anchor_dropout = 0.0,
                                       anchor_jitter_variability = NULL,
                                       affine_shift_prob = 1.0,
                                       stripe_noise_config = NULL,
                                       rt_warp_config = NULL,
                                       device = if (torch::cuda_is_available()) "cuda" else "cpu",
                                       num_threads = NULL,
                                       num_workers = 1L,
                                       save_path = NULL,
                                       resume_from = NULL,
                                       seed = 42L,
                                       verbose = TRUE) {
  if (is.null(catalog) && is.null(peak_params)) {
    stop("pretrain_denoising_online: pass either peak_params or catalog.")
  }
  if (!is.null(catalog)) {
    stopifnot(inherits(catalog, "peak_catalog"))
  }
  # Noisy three-layer mode requires catalog + anchor_ids. variable_config
  # and contamination_config are then optional and default to NULL each.
  use_noisy <- !is.null(anchor_ids)
  if (use_noisy) {
    if (is.null(catalog))
      stop("pretrain_denoising_online: anchor_ids requires catalog.")
    if (!all(anchor_ids %in% catalog$compounds$compound_id)) {
      bad <- setdiff(anchor_ids, catalog$compounds$compound_id)
      stop("anchor_ids references compound_id(s) not in the catalog: ",
           paste(head(bad, 5), collapse = ", "),
           if (length(bad) > 5) "..." else "")
    }
  }
  if ((!is.null(variable_config) || !is.null(contamination_config)) &&
      !use_noisy) {
    stop("variable_config / contamination_config require anchor_ids.")
  }
  # Split-layer denoising: when clean_layers or noise_layers is set,
  # compose the clean target and the noisy input from named layers of
  # the returned sample rather than using the default pair$clean /
  # pair$noisy. Enables true denoising (contamination as corruption
  # the encoder must strip, not signal to preserve).
  valid_layer_names <- c("anchors", "variable", "contamination")
  split_layers <- !is.null(clean_layers) || !is.null(noise_layers)
  if (split_layers) {
    if (!use_noisy)
      stop("clean_layers / noise_layers require anchor_ids ",
           "(the three-layer noisy path).")
    if (is.null(clean_layers)) clean_layers <- valid_layer_names
    if (is.null(noise_layers)) noise_layers <- character(0)
    bad_cl <- setdiff(clean_layers, valid_layer_names)
    bad_nl <- setdiff(noise_layers, valid_layer_names)
    if (length(bad_cl) || length(bad_nl))
      stop("Unknown layer name(s): ",
           paste(unique(c(bad_cl, bad_nl)), collapse = ", "),
           ". Valid: ", paste(valid_layer_names, collapse = ", "), ".")
    overlap <- intersect(clean_layers, noise_layers)
    if (length(overlap))
      stop("Layer(s) in BOTH clean_layers and noise_layers ",
           "(contradictory objective): ",
           paste(overlap, collapse = ", "), ".")
    # Warn if user names a layer whose config wasn't passed — the layer
    # will always render as zeros, which usually isn't what they meant.
    layer_configured <- c(
      anchors       = TRUE,
      variable      = !is.null(variable_config),
      contamination = !is.null(contamination_config)
    )
    named <- unique(c(clean_layers, noise_layers))
    missing_cfg <- named[!layer_configured[named]]
    if (length(missing_cfg))
      warning("clean_layers / noise_layers name(s) ",
              paste(missing_cfg, collapse = ", "),
              " but the corresponding config is NULL; those layers will ",
              "contribute zeros.")
  }
  # Contamination-regime schedule: per-sample stochastic selection of
  # {clean, light, heavy} contamination. Scale multipliers on the
  # config's n_per_sample and intensity_range are fixed (clean=0,
  # light=0.3, heavy=1.0). Named-list weights determine the sampling
  # probability; weights normalize automatically.
  if (!is.null(contamination_schedule)) {
    if (is.null(contamination_config))
      stop("contamination_schedule requires contamination_config.")
    if (!use_noisy)
      stop("contamination_schedule requires anchor_ids (three-layer path).")
    valid_regimes <- c("clean", "light", "heavy")
    bad_reg <- setdiff(names(contamination_schedule), valid_regimes)
    if (length(bad_reg))
      stop("contamination_schedule names must be from {clean, light, ",
           "heavy}; got: ", paste(bad_reg, collapse = ", "))
    w <- unlist(contamination_schedule)
    if (any(w < 0) || sum(w) <= 0)
      stop("contamination_schedule weights must be non-negative with sum > 0.")
    contamination_schedule <- as.list(w / sum(w))
  }
  # noise_scale can be a scalar (fixed sensor-noise strength) or a
  # length-2 vector (per-sample uniform draw between the two values).
  # Draw-per-sample makes the encoder robust to varying real-instrument
  # noise levels rather than assuming a fixed noise floor.
  noise_scale_range <- NULL
  if (length(noise_scale) == 2L) {
    noise_scale_range <- as.numeric(noise_scale)
    if (any(noise_scale_range < 0) ||
        noise_scale_range[1] > noise_scale_range[2])
      stop("noise_scale as a range must be c(min, max) with 0 <= min <= max.")
    noise_scale <- mean(noise_scale_range)   # for the log line
  } else if (length(noise_scale) != 1L) {
    stop("noise_scale must be scalar or length-2 numeric.")
  }
  # Masking augmentation. mask_config = NULL disables. Setting it
  # zeros out random rectangles from the NOISY input only (target
  # left unchanged) — MAE-style corruption. Encoder must reconstruct
  # the masked regions from surrounding context, which regularizes
  # against memorizing synthetic pixel patterns.
  if (!is.null(mask_config)) {
    mask_config$n_rects   <- mask_config$n_rects   %||% c(3L, 5L)
    mask_config$size_frac <- mask_config$size_frac %||% c(0.05, 0.15)
    if (length(mask_config$n_rects) == 1L)
      mask_config$n_rects <- c(mask_config$n_rects, mask_config$n_rects)
    if (any(mask_config$size_frac < 0) ||
        any(mask_config$size_frac > 1))
      stop("mask_config$size_frac must be within [0, 1].")
  }
  # Affine shift augmentation. Draws per-sample dRT ~ U(-affine_shift_rt,
  # affine_shift_rt) and dCV ~ U(-affine_shift_cv, affine_shift_cv) and
  # coherently shifts BOTH the clean target and noisy input by the same
  # amount (zero-fill at newly-exposed edges). Data augmentation, not
  # drift correction — the encoder learns that a compound at row r is
  # the same compound as one at row r+k, without penalizing itself.
  affine_shift_rt <- as.integer(affine_shift_rt)
  affine_shift_cv <- as.integer(affine_shift_cv)
  if (affine_shift_rt < 0 || affine_shift_cv < 0)
    stop("affine_shift_rt / affine_shift_cv must be >= 0.")
  use_affine <- affine_shift_rt > 0L || affine_shift_cv > 0L
  # Bernoulli gate on the affine shift — probability that a given
  # sample gets a shift at all. When 1.0 (default) every sample gets
  # a fresh dr/dc; when < 1, that fraction of samples gets NO shift
  # (dr = dc = 0). Simulates real acquisitions where alignment is
  # usually good and only some runs drift.
  if (affine_shift_prob < 0 || affine_shift_prob > 1)
    stop("affine_shift_prob must be in [0, 1].")
  # anchor_dropout validation (passthrough — sample-time enforcement is
  # inside generate_noisy_pretrain_sample).
  if (anchor_dropout < 0 || anchor_dropout >= 1)
    stop("anchor_dropout must be in [0, 1).")
  # stripe_noise_config validation.
  if (!is.null(stripe_noise_config)) {
    stripe_noise_config$prob            <- stripe_noise_config$prob            %||% 0.2
    stripe_noise_config$direction       <- stripe_noise_config$direction       %||% "rt"
    stripe_noise_config$intensity       <- stripe_noise_config$intensity       %||% c(0.05, 0.15)
    if (!stripe_noise_config$direction %in% c("rt", "cv"))
      stop("stripe_noise_config$direction must be 'rt' or 'cv'.")
    if (any(stripe_noise_config$intensity < 0))
      stop("stripe_noise_config$intensity must be non-negative.")
  }
  # rt_warp_config validation. strength expresses the fractional shift
  # applied at the LAST row (row H); the shift ramps linearly from 0
  # at row 1 up to strength * H. Applied identically to clean and noisy.
  if (!is.null(rt_warp_config)) {
    rt_warp_config$strength <- rt_warp_config$strength %||% 0.02
    rt_warp_config$prob     <- rt_warp_config$prob     %||% 1.0
    if (abs(rt_warp_config$strength) >= 0.5)
      stop("rt_warp_config$strength should be a small fraction (|strength| < 0.5).")
    if (rt_warp_config$prob < 0 || rt_warp_config$prob > 1)
      stop("rt_warp_config$prob must be in [0, 1].")
  }
  # Resolve noisy_dust_threshold once at entry:
  #   NULL default:
  #     split-layer mode -> 0    (keep contamination/noise in noisy)
  #     standard mode    -> same as dust_threshold (backward compat)
  #   explicit scalar   -> honored as-is
  if (is.null(noisy_dust_threshold)) {
    noisy_dust_threshold <- if (split_layers) 0 else dust_threshold
  } else {
    noisy_dust_threshold <- as.numeric(noisy_dust_threshold)
    if (length(noisy_dust_threshold) != 1L || noisy_dust_threshold < 0)
      stop("noisy_dust_threshold must be a non-negative scalar or NULL.")
  }
  if (!is.null(num_threads)) {
    torch::torch_set_num_threads(as.integer(num_threads))
  }
  # Parallel data-generation workers.
  #
  # Fork-based `parallel::mclapply` is unsafe here on macOS: libtorch's
  # internal threadpool mutexes don't survive fork(), and the child
  # processes hit `mutex lock failed: Invalid argument` on the first
  # batch. Use `future.apply::future_lapply` with `future::multisession`
  # instead — fresh R processes (not fork), globals resolved
  # automatically, no shared libtorch state to poison. Slightly higher
  # per-worker startup than fork on Linux but negligible over a
  # multi-epoch pretraining run.
  num_workers <- as.integer(num_workers)
  use_parallel <- num_workers > 1L
  if (use_parallel) {
    if (!requireNamespace("future", quietly = TRUE) ||
        !requireNamespace("future.apply", quietly = TRUE)) {
      message("num_workers > 1 requires the `future` and ",
              "`future.apply` packages. Install with ",
              "`install.packages(c('future', 'future.apply'))`. ",
              "Falling back to serial.")
      use_parallel <- FALSE
      num_workers <- 1L
    } else {
      # Multisession = fresh R processes on any platform. Save the
      # caller's plan so we don't clobber their outer future config.
      old_plan <- future::plan(future::multisession, workers = num_workers)
      on.exit(future::plan(old_plan), add = TRUE)
    }
  }
  location_mode  <- match.arg(location_mode)
  attribute_mode <- match.arg(attribute_mode)
  set.seed(seed); torch::torch_manual_seed(seed)

  total_samples <- as.numeric(steps_per_epoch) * epochs * batch_size
  message("Online denoising pretraining:")
  message("  ", steps_per_epoch, " steps/epoch x ", epochs, " epochs x ",
          batch_size, " batch = ",
          format(total_samples, big.mark = ",", scientific = FALSE),
          " total samples seen")
  message("  Memory footprint: ~1 batch at a time (",
          round(batch_size * H * W * 4 * 2 / 1e6, 1),
          " MB per batch, clean+noisy)")
  message("  Noise: ", if (add_noise) "enabled (denoising AE)"
          else "disabled (reconstruction AE)",
          if (add_noise) paste0(" (noise_scale = ", noise_scale, ")")
          else "")
  message("  CPU threads (torch): ", torch::torch_get_num_threads())
  message("  Encoder channels: [",
          paste(encoder_channels, collapse = ", "),
          "] (latent dim = ", tail(encoder_channels, 1L), ")")
  message("  Encoder dropout: ", encoder_dropout,
          if (encoder_dropout > 0)
            " (nn_dropout2d between residual block groups)"
          else " (off)")
  message("  R-side data workers: ", num_workers,
          if (use_parallel) " (parallel via future::multisession)" else " (serial)")
  message("  Source: ",
          if (use_noisy)
            paste0("peak_catalog three-layer (", nrow(catalog$compounds),
                   " compounds, ", length(anchor_ids), " anchors",
                   if (!is.null(variable_config))     ", +variable"     else "",
                   if (!is.null(contamination_config)) ", +contamination" else "",
                   ")")
          else if (!is.null(catalog))
            paste0("peak_catalog (", nrow(catalog$compounds), " compounds)")
          else "peak_params")
  if (split_layers) {
    message("  Denoising split — clean_layers: {",
            paste(clean_layers, collapse = ", "),
            "}, noise_layers: {",
            if (length(noise_layers)) paste(noise_layers, collapse = ", ")
            else "(none)",
            "}, sensor_noise: noisy-only")
  }
  if (!is.null(contamination_schedule)) {
    parts <- paste(names(contamination_schedule),
                    sprintf("=%.2f", unlist(contamination_schedule)),
                    sep = "", collapse = ", ")
    message("  Contamination schedule (per-sample): ", parts,
            "  (scale multipliers: clean=0x, light=0.3x, heavy=1x)")
  }
  if (!is.null(noise_scale_range)) {
    message("  noise_scale (per-sample uniform): [",
            noise_scale_range[1], ", ", noise_scale_range[2], "]")
  }
  if (!is.null(mask_config)) {
    message("  Mask augmentation: ",
            mask_config$n_rects[1], "-", mask_config$n_rects[2],
            " rects/sample, size_frac [",
            mask_config$size_frac[1], ", ", mask_config$size_frac[2],
            "] (zeros noisy input, target unchanged)")
  }
  if (use_affine) {
    message("  Affine shift augmentation: dRT +-", affine_shift_rt,
            " px, dCV +-", affine_shift_cv,
            " px (prob = ", affine_shift_prob,
            ", coherent shift of clean + noisy)")
  }
  if (anchor_dropout > 0) {
    message("  Anchor dropout: ", anchor_dropout,
            " (per-sample Bernoulli drop of anchor compounds)")
  }
  if (!is.null(anchor_jitter_variability)) {
    message("  Anchor jitter variability: [",
            anchor_jitter_variability[1], ", ",
            anchor_jitter_variability[2],
            "]x (per-sample scalar multiplier on location_jitter_rt/cv)")
  }
  if (!is.null(stripe_noise_config)) {
    message("  Stripe noise: prob = ", stripe_noise_config$prob,
            ", direction = ", stripe_noise_config$direction,
            ", intensity = [", stripe_noise_config$intensity[1], ", ",
            stripe_noise_config$intensity[2], "] (noisy-only)")
  }
  if (!is.null(rt_warp_config)) {
    message("  RT warp: strength = ", rt_warp_config$strength,
            ", prob = ", rt_warp_config$prob,
            " (linear-ramp, applied to clean + noisy)")
  }
  message("  size_jitter: ", size_jitter)
  message("  dust_threshold (clean target): ", dust_threshold,
          if (dust_threshold > 0)
            " (synthetic samples will match real-data sparsity)" else
            " (no dust thresholding; synthetic will be dense)")
  message("  noisy_dust_threshold (noisy input): ", noisy_dust_threshold,
          if (noisy_dust_threshold == 0)
            " (noisy input unclipped — contamination + sensor noise pass through)"
          else "")
  message("  location_mode: ", location_mode,
          if (location_mode == "empirical")
            paste0(" (jitter RT=", location_jitter_rt,
                    " px, CV=", location_jitter_cv, " px)")
          else "")
  message("  attribute_mode: ", attribute_mode,
          if (attribute_mode == "joint")
            " (per-peak sigma_rt / sigma_cv / intensity resampled as observed triple)"
          else " (per-peak morphology drawn independently)")
  if (val_n > 0L) message("  Validation set: ", val_n, " fixed synthetic samples")
  if (!is.null(val_real))
    message("  Real validation set: ", length(val_real), " held-out real samples")
  if (!is.null(save_path) && checkpoint_every > 0)
    message("  Checkpointing every ", checkpoint_every, " epochs to ", save_path)

  model <- dms_denoising_autoencoder(target_H = H, target_W = W,
                                      stem_stride_rt = stem_stride_rt,
                                      stem_stride_cv = stem_stride_cv,
                                      encoder_channels = encoder_channels,
                                      dropout_p = encoder_dropout)
  model$to(device = device)

  # ---- Resume from saved state if requested ----
  start_epoch <- 1L
  prior_loss_history <- NULL
  prior_val_loss_history <- NULL
  prior_batch_loss_stats <- NULL

  if (!is.null(resume_from)) {
    ae_path <- sub("\\.pt$", "_autoencoder.pt", resume_from)
    manifest_path <- sub("\\.pt$", "_manifest.Rdata", resume_from)
    if (!file.exists(ae_path)) {
      stop("resume_from: autoencoder file not found at ", ae_path,
           " (need *_autoencoder.pt next to the encoder *.pt to resume)")
    }
    message("Resuming from: ", ae_path)
    model <- torch::torch_load(ae_path)
    model$to(device = device)

    if (file.exists(manifest_path)) {
      e <- new.env()
      load(manifest_path, envir = e)
      tm <- e$training_manifest
      if (!is.null(tm$last_epoch_completed)) {
        start_epoch <- as.integer(tm$last_epoch_completed) + 1L
        prior_loss_history     <- tm$loss_history
        prior_val_loss_history <- tm$val_loss_history
        prior_batch_loss_stats <- tm$batch_loss_stats
        message("  Manifest shows last_epoch_completed = ",
                tm$last_epoch_completed,
                "; continuing from epoch ", start_epoch)
        # Verify that training-affecting args on THIS call match the ones
        # recorded in the manifest. Mismatch means continuing with a
        # DIFFERENT training distribution than what the model was trained
        # on so far — usually a mistake. Warn per mismatch so the user
        # can decide whether to abort and re-run.
        current_vals <- list(
          H = H, W = W, batch_size = batch_size, lr = lr,
          weight_decay = weight_decay, add_noise = add_noise,
          size_jitter = size_jitter, dust_threshold = dust_threshold,
          noisy_dust_threshold = noisy_dust_threshold,
          location_mode = location_mode,
          location_jitter_rt = location_jitter_rt,
          location_jitter_cv = location_jitter_cv,
          attribute_mode = attribute_mode,
          noise_scale = noise_scale,
          noise_scale_range = noise_scale_range,
          anchor_ids = anchor_ids,
          variable_config = variable_config,
          contamination_config = contamination_config,
          contamination_schedule = contamination_schedule,
          clean_layers = clean_layers,
          noise_layers = noise_layers,
          mask_config = mask_config,
          affine_shift_rt = affine_shift_rt,
          affine_shift_cv = affine_shift_cv,
          affine_shift_prob = affine_shift_prob,
          anchor_dropout = anchor_dropout,
          anchor_jitter_variability = anchor_jitter_variability,
          stripe_noise_config = stripe_noise_config,
          rt_warp_config = rt_warp_config,
          stem_stride_rt = stem_stride_rt,
          stem_stride_cv = stem_stride_cv,
          encoder_channels = encoder_channels,
          encoder_dropout = encoder_dropout
        )
        hp <- tm$hyperparams
        mismatches <- character(0)
        for (nm in names(current_vals)) {
          if (nm %in% names(hp) && !identical(current_vals[[nm]], hp[[nm]])) {
            mismatches <- c(mismatches, nm)
          }
        }
        if (length(mismatches)) {
          warning("resume_from: the following args differ between the ",
                  "current call and the saved manifest: ",
                  paste(mismatches, collapse = ", "),
                  ". Continuing with the CURRENT-call values (the model ",
                  "was trained so far with the manifest values). If you ",
                  "intended to resume identically, re-invoke with the ",
                  "manifest's args (they're in tm$hyperparams).")
        }
        if (start_epoch > epochs) {
          message("  Already at or past requested epochs; nothing to do.")
          return(list(encoder = model$encoder, autoencoder = model,
                       loss_history = tm$loss_history,
                       val_loss_history = tm$val_loss_history,
                       real_val_loss_history = tm$real_val_loss_history,
                       batch_loss_stats = tm$batch_loss_stats,
                       total_samples = tm$total_samples))
        }
      }
    } else {
      message("  No manifest found at ", manifest_path,
              "; restarting epoch counter at 1 with loaded weights.")
    }
  }

  # Build optimizer AFTER any resume so it binds to the (possibly loaded)
  # model's parameters. Note: Adam moment estimates are not restored
  # from disk; resumed training effectively re-warms Adam from scratch.
  opt <- torch::optim_adam(model$parameters, lr = lr, weight_decay = weight_decay)
  loss_fn <- torch::nn_mse_loss()

  # generate_batch: builds clean + noisy synthetic batches. When
  # use_parallel = TRUE, the per-sample loop is dispatched to
  # `num_workers` fresh R processes via future.apply::future_lapply
  # backed by future::multisession. Each task gets an independent
  # L'Ecuyer-CMRG substream (via future.seed = TRUE) so worker samples
  # are independent and reproducible.
  generate_one_pair <- function(i) {
    # Per-sample noise_scale draw (if range specified).
    ns <- if (!is.null(noise_scale_range))
            stats::runif(1, noise_scale_range[1], noise_scale_range[2])
          else noise_scale
    # Per-sample contamination-regime draw. NULL schedule -> config
    # passed through unchanged. Otherwise resolve to a scaled config:
    #   clean -> NULL (no contamination layer this sample)
    #   light -> 0.3x n_per_sample and intensity_range
    #   heavy -> full config as-is
    active_cc <- contamination_config
    if (!is.null(contamination_schedule)) {
      regime <- sample(names(contamination_schedule), size = 1L,
                        prob = unlist(contamination_schedule))
      scale <- switch(regime, clean = 0, light = 0.3, heavy = 1.0)
      if (scale == 0) {
        active_cc <- NULL
      } else if (scale < 1) {
        active_cc$n_per_sample <-
          pmax(0L, as.integer(active_cc$n_per_sample * scale))
        active_cc$intensity_range <- active_cc$intensity_range * scale
      }
    }
    pair <- if (use_noisy) {
      generate_noisy_pretrain_sample(catalog, H, W,
                                       anchor_ids = anchor_ids,
                                       variable_config = variable_config,
                                       contamination_config = active_cc,
                                       add_noise = add_noise,
                                       size_jitter = size_jitter,
                                       location_jitter_rt = location_jitter_rt,
                                       location_jitter_cv = location_jitter_cv,
                                       noise_scale = ns,
                                       anchor_dropout = anchor_dropout,
                                       anchor_jitter_variability = anchor_jitter_variability)
    } else if (!is.null(catalog)) {
      generate_one_synthetic_from_catalog(catalog, H, W,
                                    add_noise = add_noise,
                                    size_jitter = size_jitter,
                                    location_jitter_rt = location_jitter_rt,
                                    location_jitter_cv = location_jitter_cv,
                                    noise_scale = ns)
    } else {
      generate_one_synthetic(peak_params, H, W,
                                    add_noise = add_noise,
                                    size_jitter = size_jitter,
                                    location_mode = location_mode,
                                    location_jitter_rt = location_jitter_rt,
                                    location_jitter_cv = location_jitter_cv,
                                    attribute_mode = attribute_mode,
                                    noise_scale = ns)
    }
    # Split-layer denoising: compose Z_clean and Z_noisy from the
    # per-layer matrices returned by generate_noisy_pretrain_sample.
    # Sensor noise is always in Z_noisy only (i.e. always corruption)
    # so the encoder can never learn to preserve it.
    if (split_layers) {
      zero <- matrix(0, nrow = H, ncol = W)
      add_layers <- function(names_vec) {
        Reduce(`+`, pair$layers[names_vec], init = zero)
      }
      Z_clean       <- add_layers(clean_layers)
      Z_noise_add   <- add_layers(noise_layers)
      Z_sensor      <- pair$layers$sensor_noise
      pair$clean <- Z_clean
      pair$noisy <- Z_clean + Z_noise_add + Z_sensor
    }
    # Coherent per-sample affine shift (data augmentation). Applied
    # to BOTH clean and noisy so the encoder is not asked to correct
    # the drift — it just sees peaks at different absolute row/col
    # positions across samples. Bernoulli-gated by affine_shift_prob:
    # at 1.0 (default) every sample gets a shift; at 0.6 only 60% do.
    if (use_affine && stats::runif(1) < affine_shift_prob) {
      dr <- if (affine_shift_rt > 0L)
              sample(-affine_shift_rt:affine_shift_rt, 1L) else 0L
      dc <- if (affine_shift_cv > 0L)
              sample(-affine_shift_cv:affine_shift_cv, 1L) else 0L
      shift_mat <- function(M, dr, dc) {
        out <- matrix(0, nrow = nrow(M), ncol = ncol(M))
        src_r <- max(1L, 1L - dr):min(nrow(M), nrow(M) - dr)
        src_c <- max(1L, 1L - dc):min(ncol(M), ncol(M) - dc)
        if (length(src_r) == 0L || length(src_c) == 0L) return(out)
        out[src_r + dr, src_c + dc] <- M[src_r, src_c]
        out
      }
      pair$clean <- shift_mat(pair$clean, dr, dc)
      pair$noisy <- shift_mat(pair$noisy, dr, dc)
    }
    # Non-linear RT warp: linear-ramp warp along RT applied to BOTH
    # clean and noisy identically. Simulates smooth column-condition
    # drift (later peaks shift more than early peaks). Nearest-neighbor
    # row indexing so it's O(H) per sample.
    if (!is.null(rt_warp_config) &&
        stats::runif(1) < rt_warp_config$prob) {
      s <- rt_warp_config$strength
      # Randomize sign per sample so warp can compress OR stretch
      s <- s * sample(c(-1, 1), 1L)
      r_in <- seq_len(H)
      # Sample source row for each output row: linear ramp from 0
      # (top) to s*H (bottom). Nearest-neighbor with clipping.
      src <- as.integer(round(r_in + s * (r_in - 1)))
      src <- pmax(1L, pmin(H, src))
      pair$clean <- pair$clean[src, , drop = FALSE]
      pair$noisy <- pair$noisy[src, , drop = FALSE]
    }
    # Mask augmentation: zero out random rectangles from the NOISY
    # input only. Encoder reconstructs from surrounding context.
    if (!is.null(mask_config)) {
      n_rects <- if (mask_config$n_rects[1] == mask_config$n_rects[2])
                   mask_config$n_rects[1]
                 else
                   as.integer(stats::runif(1,
                     mask_config$n_rects[1],
                     mask_config$n_rects[2] + 1))
      for (k in seq_len(n_rects)) {
        area_frac <- stats::runif(1, mask_config$size_frac[1],
                                       mask_config$size_frac[2])
        # Aspect ratio ~ uniform 0.5..2.0 so rectangles vary in shape
        aspect <- exp(stats::runif(1, log(0.5), log(2.0)))
        rect_h <- max(1L, as.integer(sqrt(area_frac * H * W * aspect)))
        rect_w <- max(1L, as.integer(sqrt(area_frac * H * W / aspect)))
        rect_h <- min(rect_h, H); rect_w <- min(rect_w, W)
        r0 <- sample.int(H - rect_h + 1L, 1L)
        c0 <- sample.int(W - rect_w + 1L, 1L)
        pair$noisy[r0:(r0 + rect_h - 1L),
                    c0:(c0 + rect_w - 1L)] <- 0
      }
    }
    # Stripe noise: row- or column-correlated band-limited noise added
    # to the NOISY input only (fraction of samples set by
    # stripe_noise_config$prob). Simulates thermal / electrical
    # acquisition artifacts. RT direction gives horizontal bands
    # (constant along CV within each row); CV gives vertical bands.
    if (!is.null(stripe_noise_config) &&
        stats::runif(1) < stripe_noise_config$prob) {
      sn_amp <- stats::runif(1, stripe_noise_config$intensity[1],
                                    stripe_noise_config$intensity[2])
      if (stripe_noise_config$direction == "rt") {
        # One noise value per row, tiled across all cols
        stripe1d <- stats::rnorm(H, mean = 0, sd = sn_amp)
        stripe <- matrix(stripe1d, nrow = H, ncol = W)
      } else {
        stripe1d <- stats::rnorm(W, mean = 0, sd = sn_amp)
        stripe <- matrix(stripe1d, nrow = H, ncol = W, byrow = TRUE)
      }
      pair$noisy <- pair$noisy + stripe
    }
    # Dust thresholding matches real-data preprocessing (see
    # baseline_basement) so synthetic samples have the same sparsity
    # profile as real samples entering the encoder.
    if (dust_threshold > 0) {
      pair$clean[pair$clean < dust_threshold] <- 0
    }
    if (noisy_dust_threshold > 0) {
      pair$noisy[pair$noisy < noisy_dust_threshold] <- 0
    }
    list(
      clean = pmin(normalize_sample(pair$clean), norm_clamp),
      noisy = pmin(normalize_sample(pair$noisy), norm_clamp)
    )
  }

  generate_batch <- function(n = batch_size) {
    pairs <- if (use_parallel) {
      # future.seed = TRUE gives each parallel task an independent
      # L'Ecuyer-CMRG substream, so worker samples don't collide.
      future.apply::future_lapply(seq_len(n), generate_one_pair,
                                    future.seed = TRUE)
    } else {
      lapply(seq_len(n), generate_one_pair)
    }

    clean_arr <- array(0, dim = c(n, 1L, H, W))
    noisy_arr <- array(0, dim = c(n, 1L, H, W))
    for (i in seq_len(n)) {
      clean_arr[i, 1, , ] <- pairs[[i]]$clean
      noisy_arr[i, 1, , ] <- pairs[[i]]$noisy
    }
    list(
      clean = torch::torch_tensor(clean_arr, dtype = torch::torch_float()),
      noisy = torch::torch_tensor(noisy_arr, dtype = torch::torch_float())
    )
  }

  val_batches <- NULL
  if (val_n > 0L) {
    val_seed_state <- .Random.seed
    set.seed(seed + 1L)
    n_val_batches <- ceiling(val_n / batch_size)
    val_batches <- vector("list", n_val_batches)
    for (b in seq_len(n_val_batches)) val_batches[[b]] <- generate_batch(batch_size)
    .Random.seed <<- val_seed_state
  }

  eval_val <- function() {
    if (is.null(val_batches)) return(NA_real_)
    model$eval()
    val_loss_total <- 0
    torch::with_no_grad({
      for (vb in val_batches) {
        n <- vb$noisy$to(device = device)
        c <- vb$clean$to(device = device)
        val_loss_total <- val_loss_total + as.numeric(loss_fn(model(n), c)$item())
      }
    })
    model$train()
    val_loss_total / length(val_batches)
  }

  # Real-data validation: build torch batches from val_real once. Per-epoch
  # eval_val_real() computes MSE of autoencoder reconstruction on real
  # samples (target = input, no added noise). Non-training test that the
  # encoder-decoder actually generalizes off synthetic.
  val_real_batches <- NULL
  if (!is.null(val_real) && length(val_real) > 0L) {
    n_val_real <- length(val_real)
    n_batches_real <- ceiling(n_val_real / batch_size)
    val_real_batches <- vector("list", n_batches_real)
    for (b in seq_len(n_batches_real)) {
      s_lo <- (b - 1L) * batch_size + 1L
      s_hi <- min(b * batch_size, n_val_real)
      chunk <- val_real[s_lo:s_hi]
      n_chunk <- length(chunk)
      arr <- array(0, dim = c(n_chunk, 1L, H, W))
      for (i in seq_len(n_chunk)) {
        Z_i <- as.matrix(chunk[[i]])
        if (nrow(Z_i) != H || ncol(Z_i) != W) {
          stop("val_real entry ", s_lo + i - 1L, " has dim ",
                nrow(Z_i), "x", ncol(Z_i),
                "; expected ", H, "x", W,
                " (pad your val_real matrices to match pad_dims).")
        }
        arr[i, 1, , ] <- Z_i
      }
      val_real_batches[[b]] <- torch::torch_tensor(arr,
                                    dtype = torch::torch_float())
    }
  }

  eval_val_real <- function() {
    if (is.null(val_real_batches)) return(NA_real_)
    model$eval()
    total <- 0; n <- 0L
    torch::with_no_grad({
      for (vb in val_real_batches) {
        x <- vb$to(device = device)
        # Target = input (reconstruction check, no added noise).
        total <- total + as.numeric(loss_fn(model(x), x)$item())
        n <- n + 1L
      }
    })
    model$train()
    total / max(n, 1L)
  }

  save_manifest <- function(loss_h, val_h, real_val_h, stats, last_ep,
                             final = FALSE) {
    if (is.null(save_path)) return(invisible(NULL))
    manifest_path <- sub("\\.pt$", "_manifest.Rdata", save_path)
    training_manifest <- list(
      peak_params = peak_params,
      hyperparams = list(
        H = H, W = W, steps_per_epoch = steps_per_epoch, epochs = epochs,
        batch_size = batch_size, lr = lr, weight_decay = weight_decay,
        add_noise = add_noise, size_jitter = size_jitter,
        grad_clip = grad_clip, norm_clamp = norm_clamp, val_n = val_n,
        dust_threshold = dust_threshold,
        noisy_dust_threshold = noisy_dust_threshold,
        location_mode = location_mode,
        location_jitter_rt = location_jitter_rt,
        location_jitter_cv = location_jitter_cv,
        attribute_mode = attribute_mode,
        noise_scale = noise_scale,
        noise_scale_range = noise_scale_range,
        # Three-layer noisy / split denoising / augmentation config
        anchor_ids = anchor_ids,
        variable_config = variable_config,
        contamination_config = contamination_config,
        contamination_schedule = contamination_schedule,
        clean_layers = clean_layers,
        noise_layers = noise_layers,
        mask_config = mask_config,
        affine_shift_rt = affine_shift_rt,
        affine_shift_cv = affine_shift_cv,
        val_real_n = if (is.null(val_real)) 0L else length(val_real),
        checkpoint_every = checkpoint_every,
        stem_stride_rt = stem_stride_rt, stem_stride_cv = stem_stride_cv,
        encoder_channels = encoder_channels,
        encoder_dropout = encoder_dropout,
        anchor_dropout = anchor_dropout,
        anchor_jitter_variability = anchor_jitter_variability,
        affine_shift_prob = affine_shift_prob,
        stripe_noise_config = stripe_noise_config,
        rt_warp_config = rt_warp_config,
        device = device, seed = seed),
      loss_history = loss_h, val_loss_history = val_h,
      real_val_loss_history = real_val_h,
      batch_loss_stats = stats, total_samples = total_samples,
      last_epoch_completed = last_ep, complete = final,
      timestamp = Sys.time(), r_version = R.version.string)
    save(training_manifest, file = manifest_path)
  }

  # Initialize loss-history vectors; if resuming, splice in prior values.
  loss_history <- numeric(epochs)
  val_loss_history <- rep(NA_real_, epochs)
  real_val_loss_history <- rep(NA_real_, epochs)
  batch_loss_stats <- vector("list", epochs)
  if (!is.null(prior_loss_history)) {
    n_prior <- min(length(prior_loss_history), epochs)
    loss_history[seq_len(n_prior)]     <- prior_loss_history[seq_len(n_prior)]
    val_loss_history[seq_len(n_prior)] <- prior_val_loss_history[seq_len(n_prior)]
    batch_loss_stats[seq_len(n_prior)] <- prior_batch_loss_stats[seq_len(n_prior)]
  }
  training_start <- Sys.time()

  for (ep in seq(start_epoch, epochs)) {
    model$train()
    ep_start <- Sys.time()
    batch_losses <- numeric(steps_per_epoch)
    for (step in seq_len(steps_per_epoch)) {
      batch <- generate_batch()
      noisy <- batch$noisy$to(device = device)
      clean <- batch$clean$to(device = device)
      opt$zero_grad()
      recon <- model(noisy)
      loss <- loss_fn(recon, clean)
      loss$backward()
      if (grad_clip > 0) {
        torch::nn_utils_clip_grad_norm_(model$parameters, max_norm = grad_clip)
      }
      opt$step()
      batch_losses[step] <- as.numeric(loss$item())
    }
    batch_loss_stats[[ep]] <- list(
      mean = mean(batch_losses), median = stats::median(batch_losses),
      max = max(batch_losses), min = min(batch_losses),
      p95 = as.numeric(stats::quantile(batch_losses, 0.95)))
    loss_history[ep] <- batch_loss_stats[[ep]]$mean
    val_loss_history[ep] <- eval_val()
    real_val_loss_history[ep] <- eval_val_real()
    ep_time <- as.numeric(difftime(Sys.time(), ep_start, units = "mins"))
    elapsed_total <- as.numeric(difftime(Sys.time(), training_start, units = "mins"))
    eta_min <- elapsed_total / ep * (epochs - ep)
    if (verbose) {
      val_str <- if (val_n > 0L) paste0(" | val MSE: ",
                                          round(val_loss_history[ep], 4)) else ""
      real_val_str <- if (!is.null(val_real_batches))
        paste0(" | real val MSE: ",
                round(real_val_loss_history[ep], 4)) else ""
      diag_str <- if (loss_diagnostics) paste0(
        " | median: ", round(batch_loss_stats[[ep]]$median, 4),
        " | p95: ",    round(batch_loss_stats[[ep]]$p95,    4),
        " | max: ",    round(batch_loss_stats[[ep]]$max,    4)) else ""
      message("Epoch ", ep, "/", epochs,
              " | train MSE: ", round(batch_loss_stats[[ep]]$mean, 4),
              val_str, real_val_str, diag_str,
              " | Time: ", round(ep_time, 2), " min",
              " | ETA: ", round(eta_min, 0), " min")
    }
    save_manifest(loss_history, val_loss_history, real_val_loss_history,
                   batch_loss_stats,
                   last_ep = ep, final = (ep == epochs))

    # After EVERY epoch: overwrite encoder.pt and autoencoder.pt with the
    # current state so a mid-run crash always leaves recoverable
    # artifacts. The files are small (~ a few MB) so the per-epoch I/O
    # cost is negligible (milliseconds) even for million-sample runs.
    if (!is.null(save_path)) {
      torch::torch_save(model$encoder, save_path)
      ae_path <- sub("\\.pt$", "_autoencoder.pt", save_path)
      if (ae_path == save_path) ae_path <- paste0(save_path, "_autoencoder.pt")
      torch::torch_save(model, ae_path)
    }

    # Additionally: save a numbered checkpoint every `checkpoint_every`
    # epochs so the user can track training progression and roll back
    # to a specific earlier epoch if needed.
    if (!is.null(save_path) && checkpoint_every > 0L &&
        ep %% checkpoint_every == 0L && ep < epochs) {
      ckpt_path <- sub("\\.pt$", paste0("_epoch", ep, ".pt"), save_path)
      torch::torch_save(model$encoder, ckpt_path)
      if (verbose) message("  -> Checkpoint saved: ", ckpt_path)
    }

    # Memory hygiene: explicitly free per-epoch tensors and nudge R to
    # collect. Helps prevent slow memory growth across long runs and
    # reduces the chance of OS-level memory-pressure kills, especially
    # on machines with limited RAM (e.g., laptops).
    rm(batch, noisy, clean, recon, loss, batch_losses)
    gc(verbose = FALSE, full = FALSE)
  }

  if (!is.null(save_path)) {
    message("Encoder saved to: ", save_path)
    ae_path <- sub("\\.pt$", "_autoencoder.pt", save_path)
    if (ae_path == save_path) ae_path <- paste0(save_path, "_autoencoder.pt")
    message("Autoencoder saved to: ", ae_path)
    manifest_path <- sub("\\.pt$", "_manifest.Rdata", save_path)
    message("Manifest saved to: ", manifest_path)
  }

  list(encoder = model$encoder, autoencoder = model,
       loss_history = loss_history, val_loss_history = val_loss_history,
       real_val_loss_history = real_val_loss_history,
       batch_loss_stats = batch_loss_stats, total_samples = total_samples)
}

#' Pre-train the encoder from a peak catalog using three-layer noisy samples (smoke test only)
#'
#' Same training loop as [pretrain_autoencoder_from_catalog()] but each
#' synthetic sample is generated via
#' [generate_noisy_pretrain_sample()] — an anchor backbone plus
#' optional variable-catalog and contamination layers. The extra
#' variability is intended to teach the encoder to expect real-world
#' skin-VOC characteristics: consistent backbone, per-sample dietary /
#' host-genetic peaks, and contamination-like unfamiliar peaks.
#'
#' **For real pretraining runs prefer
#' [pretrain_denoising_online()] with `catalog = your_catalog`,
#' `anchor_ids = ids`, and the same `variable_config` /
#' `contamination_config`** — that path is online (fresh batches per
#' step) so training data isn't capped at `n_synthetic`. This pre-alloc
#' version stays available for smoke tests.
#'
#' To reproduce the current `pretrain_autoencoder_from_catalog()`
#' behavior, pass `variable_config = NULL` and `contamination_config = NULL`
#' — the samples will contain only the anchor set (fired at 100 percent).
#'
#' @param catalog A `peak_catalog` object.
#' @param H,W Synthetic image dimensions.
#' @param anchor_ids Integer vector of compound_ids to force-fire in
#'   every sample (the invariant backbone). Typically the ~100 anchors
#'   at prevalence >= 0.9 from the universal catalog.
#' @param variable_config,contamination_config Passed through to
#'   [generate_noisy_pretrain_sample()]. Pass `NULL` for either to
#'   skip that layer. See that function's docs for field details.
#' @param n_synthetic Number of synthetic samples generated.
#' @param epochs,batch_size,lr,weight_decay Optimization parameters.
#' @param add_noise If TRUE, denoising autoencoder; if FALSE,
#'   reconstruction.
#' @param stem_stride_rt,stem_stride_cv Encoder stem strides.
#' @param device "cpu" or "cuda".
#' @param seed RNG seed.
#' @return List with `encoder`, `autoencoder`, `loss_history`, and the
#'   `variable_config` / `contamination_config` actually used
#'   (echoed for provenance).
#' @export
pretrain_autoencoder_from_catalog_noisy <- function(
  catalog, H, W,
  anchor_ids,
  variable_config = NULL,
  contamination_config = NULL,
  n_synthetic = 2000L,
  epochs = 30L,
  batch_size = 32L,
  lr = 1e-3,
  weight_decay = 1e-4,
  add_noise = TRUE,
  stem_stride_rt = 4L,
  stem_stride_cv = 1L,
  device = if (torch::cuda_is_available()) "cuda" else "cpu",
  seed = 42L
) {
  set.seed(seed); torch::torch_manual_seed(seed)
  stopifnot(inherits(catalog, "peak_catalog"))
  message("Generating ", n_synthetic, " noisy synthetic samples ",
          "(", H, "x", W, ", ", length(anchor_ids), " anchors",
          if (!is.null(variable_config))       ", +variable"    else "",
          if (!is.null(contamination_config)) ", +contamination" else "",
          ")...")

  clean_arr <- array(0, dim = c(n_synthetic, 1L, H, W))
  noisy_arr <- array(0, dim = c(n_synthetic, 1L, H, W))
  for (i in seq_len(n_synthetic)) {
    res <- generate_noisy_pretrain_sample(catalog, H = H, W = W,
      anchor_ids = anchor_ids,
      variable_config = variable_config,
      contamination_config = contamination_config,
      add_noise = add_noise)
    clean_arr[i, 1L, , ] <- res$clean
    noisy_arr[i, 1L, , ] <- if (add_noise) res$noisy else res$clean
    if (i %% 200L == 0L) message("  generated ", i, "/", n_synthetic)
  }

  q_clean <- as.numeric(stats::quantile(clean_arr[clean_arr > 0], 0.95,
                                          na.rm = TRUE))
  if (!is.finite(q_clean) || q_clean == 0) q_clean <- 1
  clean_arr[clean_arr < 0] <- 0
  noisy_arr[noisy_arr < 0] <- 0
  clean_arr <- log1p(clean_arr) / log1p(q_clean)
  noisy_arr <- log1p(noisy_arr) / log1p(q_clean)

  clean_t <- torch::torch_tensor(clean_arr, dtype = torch::torch_float())
  noisy_t <- torch::torch_tensor(noisy_arr, dtype = torch::torch_float())

  ds_class <- torch::dataset(
    name = "noisy_catalog_pretrain_ds",
    initialize = function() {},
    .getitem = function(i) list(input = noisy_t[i, ..],
                                 target = clean_t[i, ..]),
    .length = function() clean_t$size()[1]
  )
  dl <- torch::dataloader(ds_class(), batch_size = batch_size,
                            shuffle = TRUE)

  model <- dms_denoising_autoencoder(target_H = H, target_W = W,
                                       stem_stride_rt = stem_stride_rt,
                                       stem_stride_cv = stem_stride_cv)
  model$to(device = device)
  opt <- torch::optim_adam(model$parameters, lr = lr,
                            weight_decay = weight_decay)
  loss_fn <- torch::nn_mse_loss()
  loss_history <- numeric(epochs)

  for (ep in seq_len(epochs)) {
    model$train(); epoch_loss <- 0; n_batches <- 0
    coro::loop(for (batch in dl) {
      input <- batch$input$to(device = device)
      target <- batch$target$to(device = device)
      opt$zero_grad()
      recon <- model(input)
      loss <- loss_fn(recon, target)
      loss$backward(); opt$step()
      epoch_loss <- epoch_loss + as.numeric(loss$item())
      n_batches <- n_batches + 1
    })
    loss_history[ep] <- epoch_loss / max(1, n_batches)
    if (ep %% 10 == 0 || ep == 1) {
      message("Epoch ", ep, "/", epochs, " | MSE: ",
              round(loss_history[ep], 6))
    }
  }
  list(encoder = model$encoder, autoencoder = model,
       loss_history = loss_history,
       variable_config = variable_config,
       contamination_config = contamination_config,
       anchor_ids = anchor_ids)
}
