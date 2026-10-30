# Parallel synthetic-data workers for pretrain_denoising_online().
#
# Each batch used to ship closures defined inside pretrain_denoising_online()
# to the workers via future_lapply(). A closure serialises its whole
# enclosing frame, so every future carried the catalog, val_real, the model
# and optimiser (and, under do.call(), a second copy of each argument held
# in its promise): ~0.2-0.5 GiB per future per batch. Generation now runs
# in package-level functions that take their inputs explicitly:
#   - the static inputs (`cfg`) go to each worker once per run;
#   - each batch sends only the per-sample augment specs and RNG seeds.
# Dispatch uses a persistent PSOCK cluster and parallel::clusterApply()
# directly; the caller's future plan is never touched.

# Inputs .pretrain_pair() reads from pretrain_denoising_online(); the
# snapshot is taken once, after every one of them is final.
.pretrain_pair_fields <- c(
  "H", "W", "catalog", "peak_params", "use_noisy", "anchor_ids",
  "variable_config", "contamination_config", "contamination_schedule",
  "variable_source_peak_params", "add_noise", "noise_scale",
  "noise_scale_range", "size_jitter", "location_mode",
  "location_jitter_rt", "location_jitter_cv", "attribute_mode",
  "anchor_dropout", "anchor_draw", "anchor_jitter_variability",
  "peak_shape", "emg_tau_range", "emg_mix_prob", "split_layers",
  "clean_layers", "noise_layers", "rt_warp_config", "mask_config",
  "dust_threshold", "noisy_dust_threshold", "norm_clamp", "presence_head",
  "presence_layers", "presence_ids")

# One (clean, noisy, presence) training pair. Same steps and RNG draws, in
# the same order, as the closure it replaces, so serial output is unchanged.
.pretrain_pair <- function(augment_spec, cfg) {
  H                           <- cfg[["H"]]
  W                           <- cfg[["W"]]
  catalog                     <- cfg[["catalog"]]
  peak_params                 <- cfg[["peak_params"]]
  use_noisy                   <- cfg[["use_noisy"]]
  anchor_ids                  <- cfg[["anchor_ids"]]
  variable_config             <- cfg[["variable_config"]]
  contamination_config        <- cfg[["contamination_config"]]
  contamination_schedule      <- cfg[["contamination_schedule"]]
  variable_source_peak_params <- cfg[["variable_source_peak_params"]]
  add_noise                   <- cfg[["add_noise"]]
  noise_scale                 <- cfg[["noise_scale"]]
  noise_scale_range           <- cfg[["noise_scale_range"]]
  size_jitter                 <- cfg[["size_jitter"]]
  location_mode               <- cfg[["location_mode"]]
  location_jitter_rt          <- cfg[["location_jitter_rt"]]
  location_jitter_cv          <- cfg[["location_jitter_cv"]]
  attribute_mode              <- cfg[["attribute_mode"]]
  anchor_dropout              <- cfg[["anchor_dropout"]]
  anchor_draw                 <- cfg[["anchor_draw"]]
  anchor_jitter_variability   <- cfg[["anchor_jitter_variability"]]
  peak_shape                  <- cfg[["peak_shape"]]
  emg_tau_range               <- cfg[["emg_tau_range"]]
  emg_mix_prob                <- cfg[["emg_mix_prob"]]
  split_layers                <- cfg[["split_layers"]]
  clean_layers                <- cfg[["clean_layers"]]
  noise_layers                <- cfg[["noise_layers"]]
  rt_warp_config              <- cfg[["rt_warp_config"]]
  mask_config                 <- cfg[["mask_config"]]
  dust_threshold              <- cfg[["dust_threshold"]]
  noisy_dust_threshold        <- cfg[["noisy_dust_threshold"]]
  norm_clamp                  <- cfg[["norm_clamp"]]
  presence_head               <- cfg[["presence_head"]]
  presence_layers             <- cfg[["presence_layers"]]
  presence_ids                <- cfg[["presence_ids"]]
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
      # Fill generate_noisy_pretrain_sample()'s defaults before scaling:
      # NULL * scale is numeric(0), which would slip past its `%||%`.
      active_cc$n_per_sample <-
        pmax(0L, as.integer((active_cc$n_per_sample %||% c(5, 30)) * scale))
      active_cc$intensity_range <-
        (active_cc$intensity_range %||% c(0.1, 0.5)) * scale
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
                                     anchor_draw = anchor_draw,
                                     anchor_jitter_variability = anchor_jitter_variability,
                                     variable_source_peak_params = variable_source_peak_params,
                                     peak_shape = peak_shape,
                                     emg_tau_range = emg_tau_range,
                                     emg_mix_prob = emg_mix_prob)
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
  # Coherent affine shift (data augmentation). Draw is precomputed
  # in generate_batch so the same (dr, dc) can be shared across a
  # shift group (shift_group_size > 1). Applied to BOTH clean and
  # noisy so the encoder is not asked to correct the drift — it just
  # sees peaks at different absolute row/col positions across samples.
  if (isTRUE(augment_spec$shift$active) &&
      (augment_spec$shift$dr != 0L || augment_spec$shift$dc != 0L)) {
    dr <- augment_spec$shift$dr; dc <- augment_spec$shift$dc
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
  # to the NOISY input only. The stripe matrix is precomputed in
  # generate_batch: per-sample by default, or once per shift group
  # when stripe_noise_config$grouped = TRUE. Simulates thermal /
  # electrical acquisition artifacts.
  if (!is.null(augment_spec$stripe)) {
    pair$noisy <- pair$noisy + augment_spec$stripe
  }
  # Baseline drift: smooth low-frequency baseline added to NOISY
  # input only (clean target untouched). The encoder learns to strip
  # residual baseline that AsLS didn't fully remove.
  if (!is.null(augment_spec$baseline)) {
    pair$noisy <- pair$noisy + augment_spec$baseline
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
  presence_vec <- NULL
  if (presence_head) {
    present_ids <- unlist(pair$rendered_ids[presence_layers],
                          use.names = FALSE)
    presence_vec <- as.numeric(presence_ids %in% present_ids)
  }
  list(
    clean = pmin(normalize_sample(pair$clean), norm_clamp),
    noisy = pmin(normalize_sample(pair$noisy), norm_clamp),
    presence = presence_vec
  )
}

# Static inputs for the current run on a worker process. Lives in the
# namespace (not the worker's global env) so it persists across futures.
.pretrain_worker_state <- new.env(parent = emptyenv())

# Fingerprint of this namespace's function bodies. The main session sends
# its own; a worker whose installed classydms differs (e.g. the main
# session uses devtools::load_all()) would generate data with other code.
.pretrain_code_fingerprint <- function() {
  ns <- asNamespace("classydms")
  nms <- sort(ls(ns, all.names = TRUE))
  txt <- unlist(lapply(nms, function(n) {
    f <- get(n, envir = ns)
    if (is.function(f)) c(n, deparse(f, control = NULL)) else NULL
  }))
  path <- tempfile(fileext = ".txt")
  on.exit(unlink(path))
  writeLines(txt, path)
  unname(tools::md5sum(path))
}

# Runs on each worker once per run (via parallel::clusterCall).
.pretrain_worker_init <- function(cfg, run_id, fingerprint) {
  ns <- asNamespace("classydms")
  if (!exists(".pretrain_worker_state", envir = ns, inherits = FALSE))
    stop("the worker's installed classydms (",
         as.character(utils::packageVersion("classydms")),
         ") predates the parallel pretraining workers; reinstall ",
         "classydms so the workers run the same code as this session.")
  st <- get(".pretrain_worker_state", envir = ns)
  st$cfg <- cfg
  st$run_id <- run_id
  list(pid = Sys.getpid(),
       same_code = identical(get(".pretrain_code_fingerprint", envir = ns)(),
                             fingerprint))
}

# Per-sample L'Ecuyer-CMRG seeds, the same ones
# future.apply::future_lapply(X, FUN, future.seed = TRUE) assigns to n
# elements, and the same effect on the main RNG stream (advanced by two
# draws), so parallel output matches what the future-based dispatch gave.
.pretrain_rng_seeds <- function(n) {
  genv <- globalenv()
  oseed <- get(".Random.seed", envir = genv)
  if (oseed[1L] %% 10000L == 407L) {
    lec <- oseed
  } else {
    # Switching kind seeds L'Ecuyer-CMRG from the current stream; take that
    # seed and put the caller's kind and state back.
    okind <- RNGkind("L'Ecuyer-CMRG")[1L]
    lec <- get(".Random.seed", envir = genv)
    RNGkind(okind)
    assign(".Random.seed", oseed, envir = genv)
  }
  sample.int(1L, 1L)
  sample.int(1L, 1L)
  seeds <- vector("list", n)
  for (i in seq_len(n)) {
    seeds[[i]] <- parallel::nextRNGSubStream(lec)
    lec <- parallel::nextRNGStream(lec)
  }
  seeds
}

# Runs on a worker: one chunk of samples, each from its own RNG stream.
# `job` holds only this chunk's augment specs and seeds.
.pretrain_worker_chunk <- function(job, run_id) {
  st <- .pretrain_worker_state
  if (!identical(st$run_id, run_id))
    stop("pretrain worker ", Sys.getpid(), " holds no inputs for this run.")
  genv <- globalenv()
  out <- vector("list", length(job$specs))
  for (i in seq_along(job$specs)) {
    assign(".Random.seed", job$seeds[[i]], envir = genv)
    out[[i]] <- .pretrain_pair(job$specs[[i]], st$cfg)
  }
  out
}

# Generate one batch on the workers: contiguous chunks, one per worker.
.pretrain_parallel_pairs <- function(workers, augment_specs, run_id) {
  n <- length(augment_specs)
  seeds <- .pretrain_rng_seeds(n)
  jobs <- lapply(parallel::splitIndices(n, min(n, length(workers$cl))),
                 function(ii) list(specs = augment_specs[ii], seeds = seeds[ii]))
  res <- parallel::clusterApply(workers$cl, jobs, .pretrain_worker_chunk,
                                run_id = run_id)
  unlist(res, recursive = FALSE)
}

# Start `n` fresh R processes (PSOCK, never fork: libtorch's threadpool
# mutexes don't survive fork() on macOS). Workers get this session's
# library paths, as future::multisession() does.
.pretrain_start_workers <- function(n) {
  list(cl = future::makeClusterPSOCK(n, rscript_libs = .libPaths()))
}

# Send the static inputs to every worker, once per run.
.pretrain_init_workers <- function(workers, cfg, run_id) {
  res <- parallel::clusterCall(workers$cl, .pretrain_worker_init, cfg, run_id,
                               .pretrain_code_fingerprint())
  if (!all(vapply(res, `[[`, logical(1), "same_code")))
    warning("pretrain_denoising_online: the parallel workers' installed ",
            "classydms differs from the code loaded in this session (e.g. ",
            "devtools::load_all()); workers generate data with the ",
            "installed version. Reinstall classydms to make them match.",
            call. = FALSE)
  invisible(vapply(res, `[[`, integer(1), "pid"))
}

# Shut the workers down. Safe to call from on.exit() after an error or
# interrupt.
.pretrain_stop_workers <- function(workers) {
  try(parallel::stopCluster(workers$cl), silent = TRUE)
  invisible(NULL)
}
