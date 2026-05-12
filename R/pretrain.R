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
#' @param val_n Number of fixed validation samples (0 disables validation).
#' @param checkpoint_every Save encoder checkpoint every N epochs.
#' @param loss_diagnostics If TRUE, also print median, p95, max batch loss.
#' @param stem_stride_rt,stem_stride_cv Encoder stem strides.
#' @param device "cpu" or "cuda".
#' @param num_threads If non-NULL, sets the number of CPU threads
#'   used by torch for forward/backward passes via
#'   [torch::torch_set_num_threads()]. Apple Silicon machines benefit
#'   most from setting this to the number of performance cores (e.g. 4
#'   on an M-series Pro chip with 4P+10E layout), since the efficiency
#'   cores contribute less per thread. Default `NULL` leaves torch's
#'   automatic choice in place.
#' @param save_path If non-NULL, writes encoder, autoencoder, and manifest
#'   to disk after every epoch (so a mid-run crash leaves recoverable
#'   artifacts).
#' @param resume_from If non-NULL, path to a previously-saved encoder
#'   (`*.pt`). The function looks for the corresponding `_autoencoder.pt`
#'   and `_manifest.Rdata` next to it and resumes training from
#'   `last_epoch_completed + 1`, inheriting prior loss history. Useful
#'   for picking up after a crash without losing the work done so far.
#'   The other hyperparameters (`H`, `W`, `epochs`, etc.) must match
#'   the original run.
#' @param seed RNG seed.
#' @param verbose Whether to print epoch progress.
#' @return List with `encoder`, `autoencoder`, `loss_history`,
#'   `val_loss_history`, `batch_loss_stats`, and `total_samples`.
#' @export
pretrain_denoising_online <- function(peak_params, H, W,
                                       steps_per_epoch = 500L,
                                       epochs = 30L,
                                       batch_size = 32L,
                                       lr = 1e-3,
                                       weight_decay = 1e-4,
                                       add_noise = TRUE,
                                       size_jitter = 0.6,
                                       grad_clip = 1.0,
                                       norm_clamp = 10.0,
                                       val_n = 200L,
                                       checkpoint_every = 10L,
                                       loss_diagnostics = TRUE,
                                       stem_stride_rt = 4L,
                                       stem_stride_cv = 1L,
                                       device = if (torch::cuda_is_available()) "cuda" else "cpu",
                                       num_threads = NULL,
                                       save_path = NULL,
                                       resume_from = NULL,
                                       seed = 42L,
                                       verbose = TRUE) {
  if (!is.null(num_threads)) {
    torch::torch_set_num_threads(as.integer(num_threads))
  }
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
          else "disabled (reconstruction AE)")
  message("  CPU threads (torch): ", torch::torch_get_num_threads())
  message("  size_jitter: ", size_jitter)
  if (val_n > 0L) message("  Validation set: ", val_n, " fixed synthetic samples")
  if (!is.null(save_path) && checkpoint_every > 0)
    message("  Checkpointing every ", checkpoint_every, " epochs to ", save_path)

  model <- dms_denoising_autoencoder(target_H = H, target_W = W,
                                      stem_stride_rt = stem_stride_rt,
                                      stem_stride_cv = stem_stride_cv)
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
        if (start_epoch > epochs) {
          message("  Already at or past requested epochs; nothing to do.")
          return(list(encoder = model$encoder, autoencoder = model,
                       loss_history = tm$loss_history,
                       val_loss_history = tm$val_loss_history,
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

  generate_batch <- function(n = batch_size) {
    clean_arr <- array(0, dim = c(n, 1L, H, W))
    noisy_arr <- array(0, dim = c(n, 1L, H, W))
    for (i in seq_len(n)) {
      pair <- generate_one_synthetic(peak_params, H, W,
                                      add_noise = add_noise,
                                      size_jitter = size_jitter)
      clean_arr[i, 1, , ] <- pmin(normalize_sample(pair$clean), norm_clamp)
      noisy_arr[i, 1, , ] <- pmin(normalize_sample(pair$noisy), norm_clamp)
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

  save_manifest <- function(loss_h, val_h, stats, last_ep, final = FALSE) {
    if (is.null(save_path)) return(invisible(NULL))
    manifest_path <- sub("\\.pt$", "_manifest.Rdata", save_path)
    training_manifest <- list(
      peak_params = peak_params,
      hyperparams = list(
        H = H, W = W, steps_per_epoch = steps_per_epoch, epochs = epochs,
        batch_size = batch_size, lr = lr, weight_decay = weight_decay,
        add_noise = add_noise, size_jitter = size_jitter,
        grad_clip = grad_clip, norm_clamp = norm_clamp, val_n = val_n,
        checkpoint_every = checkpoint_every,
        stem_stride_rt = stem_stride_rt, stem_stride_cv = stem_stride_cv,
        device = device, seed = seed),
      loss_history = loss_h, val_loss_history = val_h,
      batch_loss_stats = stats, total_samples = total_samples,
      last_epoch_completed = last_ep, complete = final,
      timestamp = Sys.time(), r_version = R.version.string)
    save(training_manifest, file = manifest_path)
  }

  # Initialize loss-history vectors; if resuming, splice in prior values.
  loss_history <- numeric(epochs)
  val_loss_history <- rep(NA_real_, epochs)
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
    ep_time <- as.numeric(difftime(Sys.time(), ep_start, units = "mins"))
    elapsed_total <- as.numeric(difftime(Sys.time(), training_start, units = "mins"))
    eta_min <- elapsed_total / ep * (epochs - ep)
    if (verbose) {
      val_str <- if (val_n > 0L) paste0(" | val MSE: ",
                                          round(val_loss_history[ep], 4)) else ""
      diag_str <- if (loss_diagnostics) paste0(
        " | median: ", round(batch_loss_stats[[ep]]$median, 4),
        " | p95: ",    round(batch_loss_stats[[ep]]$p95,    4),
        " | max: ",    round(batch_loss_stats[[ep]]$max,    4)) else ""
      message("Epoch ", ep, "/", epochs,
              " | train MSE: ", round(batch_loss_stats[[ep]]$mean, 4),
              val_str, diag_str,
              " | Time: ", round(ep_time, 2), " min",
              " | ETA: ", round(eta_min, 0), " min")
    }
    save_manifest(loss_history, val_loss_history, batch_loss_stats,
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
       batch_loss_stats = batch_loss_stats, total_samples = total_samples)
}
