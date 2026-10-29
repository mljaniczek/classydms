# Short pretrain_denoising_online() run on the toy catalog with every
# per-sample augmentation on. Returns each batch exactly as generated: the
# arrays passed to torch::torch_tensor() (clean / noisy 4-D batches and
# presence-label matrices), in order. Pure R values, so comparable
# bit-for-bit across machines.
capture_pretrain_batches <- function(num_workers = 1L, seed = 11L,
                                     val_real = NULL, ...) {
  got <- new.env()
  suppressMessages(trace(
    "torch_tensor", where = asNamespace("torch"), print = FALSE,
    tracer = bquote(if (is.array(data) && length(dim(data)) %in% c(2L, 4L))
      assign(sprintf("%04d", length(ls(.(got))) + 1L), data, envir = .(got)))))
  on.exit(suppressMessages(untrace("torch_tensor",
                                   where = asNamespace("torch"))), add = TRUE)
  suppressMessages(pretrain_denoising_online(
    H = 128L, W = 32L, catalog = make_toy_catalog(), anchor_ids = c(1L, 2L),
    variable_config = list(draw = "prevalence", prevalence_range = c(0, 1)),
    contamination_config = list(n_per_sample = c(2L, 6L),
                                intensity_range = c(0.05, 0.2)),
    contamination_schedule = list(clean = 0.3, light = 0.4, heavy = 0.3),
    noise_scale = c(1, 3), anchor_draw = "prevalence",
    mask_config = list(n_rects = c(1L, 2L), size_frac = c(0.02, 0.05)),
    affine_shift_rt = c(1L, 4L), affine_shift_cv = c(0L, 2L),
    affine_shift_prob = 0.6, shift_group_size = c(2L, 3L),
    stripe_noise_config = list(prob = 0.5, direction = "rt",
                               intensity = c(0.05, 0.1), grouped = TRUE),
    baseline_drift_config = list(prob = 0.5, direction = c("rt", "cv"),
                                 shape = "polynomial", order = c(2L, 3L),
                                 amplitude = c(0.02, 0.05)),
    rt_warp_config = list(strength = 0.02, prob = 0.5),
    peak_shape = "mix", presence_head = TRUE,
    encoder_channels = c(4L, 4L, 8L, 8L, 8L),
    epochs = 1L, steps_per_epoch = 2L, batch_size = 5L, val_n = 5L,
    val_real = val_real, num_workers = num_workers, device = "cpu",
    num_threads = 1L, seed = seed, verbose = FALSE, ...))
  mget(sort(ls(got)), envir = got)
}

# Parallel workers load the *installed* classydms. Skip when it isn't the
# code under test (e.g. devtools::test() against an older install).
skip_if_workers_differ <- function() {
  skip_on_cran()
  cl <- future::makeClusterPSOCK(1L, rscript_libs = .libPaths())
  on.exit(parallel::stopCluster(cl))
  probe <- function()
    get(".pretrain_code_fingerprint", envir = asNamespace("classydms"))()
  environment(probe) <- globalenv()   # don't ship the test env
  fp <- tryCatch(parallel::clusterCall(cl, probe)[[1]],
                 error = function(e) NA_character_)
  skip_if_not(identical(fp, .pretrain_code_fingerprint()),
              "parallel workers load an installed classydms that differs from the code under test")
}

# Bit-level digest of each captured batch: dims plus the md5 of its values
# as little-endian doubles (independent of R's serialize() header).
batch_digests <- function(batches) {
  vapply(batches, function(x) {
    f <- tempfile()
    on.exit(unlink(f))
    writeBin(as.double(x), f, endian = "little")
    paste(paste(dim(x), collapse = "x"), unname(tools::md5sum(f)))
  }, character(1))
}
