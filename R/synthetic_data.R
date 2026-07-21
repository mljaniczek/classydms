# classydms: parametric synthetic GC-DMS data generation
# Estimates peak distributions from real data (label-agnostic) and generates
# synthetic images statistically matched to the real cohort.

#' Estimate peak parameter distributions from real GC-DMS samples
#'
#' Pools all samples (cases and controls together; no labels used) and
#' extracts empirical distributions of: number of peaks per sample, peak
#' locations (as fractions of image dimensions), peak widths (sigma_RT,
#' sigma_CV) measured at `frac_height` height and converted to Gaussian
#' sigmas, peak intensities, and background noise mean / SD. Widths are
#' fitted to log-normal distributions.
#'
#' If `Z_pretrim_list` is provided (pre-dust-threshold matrices), peak
#' width measurement uses the pre-thresholded data as a fallback to
#' preserve peak flanks that the dust thresholding zeroed out.
#'
#' @param Z_list List of trimmed, normalized intensity matrices (real data).
#' @param Z_pretrim_list Optional list of pre-dust-threshold matrices.
#' @param top_k Maximum peaks to detect per sample (default 150).
#' @param frac_height Fraction of peak height defining the width (default 0.25;
#'   0.5 corresponds to FWHM).
#' @param min_sep_rt,min_sep_cv Minimum pixel separation between peaks.
#' @param verbose Whether to print diagnostic messages (default TRUE).
#' @return A list of fitted distributions and raw observation vectors,
#'   suitable for passing to [generate_one_synthetic()] or
#'   [generate_synthetic_dataset()].
#' @export
estimate_peak_params <- function(Z_list,
                                  Z_pretrim_list = NULL,
                                  top_k = 150,
                                  frac_height = 0.25,
                                  min_sep_rt = 8,
                                  min_sep_cv = 2,
                                  verbose = TRUE) {
  stopifnot(length(Z_list) > 0)
  use_pretrim <- !is.null(Z_pretrim_list) &&
    length(Z_pretrim_list) == length(Z_list)

  all_rt_loc <- c(); all_cv_loc <- c()
  all_sigma_rt <- c(); all_sigma_cv <- c()
  all_intensity <- c(); all_n_peaks <- c()
  all_noise_mean <- c(); all_noise_sd <- c()
  fwhm_to_sigma <- 1 / (2 * sqrt(2 * log(2)))
  n_width_fail <- 0L

  for (i in seq_along(Z_list)) {
    Z <- as.matrix(Z_list[[i]]); storage.mode(Z) <- "numeric"
    R <- nrow(Z); C <- ncol(Z)
    eps_i <- robust_eps(Z)
    if (use_pretrim) {
      Z_width <- as.matrix(Z_pretrim_list[[i]])
      storage.mode(Z_width) <- "numeric"
      Rw <- min(nrow(Z_width), R); Cw <- min(ncol(Z_width), C)
      Z_width <- Z_width[seq_len(Rw), seq_len(Cw), drop = FALSE]
    } else {
      Z_width <- Z
    }
    if (verbose && i == 1) {
      zpos <- Z[Z > 0]
      message("  [estimate_peak_params] sample 1: dim=", R, "x", C,
              ", eps=", signif(eps_i, 3),
              ", n_positive=", length(zpos),
              ", Q95=", signif(stats::quantile(zpos, 0.95), 3),
              ", max=", signif(max(zpos), 3),
              ", using_pretrim=", use_pretrim)
    }
    centers <- pick_peak_centers(Z, top_k = top_k, eps = eps_i,
                                  min_sep_rt = min_sep_rt,
                                  min_sep_cv = min_sep_cv)
    all_n_peaks <- c(all_n_peaks, nrow(centers))
    if (nrow(centers) > 0) {
      max_radius_rt <- 4L; max_radius_cv <- 2L
      for (j in seq_len(nrow(centers))) {
        r0 <- centers[j, 1]; c0 <- centers[j, 2]
        h0 <- Z[r0, c0]
        w <- peak_width_at(Z, r0, c0, frac_height = frac_height, eps = 0,
                           max_radius_rt = max_radius_rt,
                           max_radius_cv = max_radius_cv)
        w_rt <- if (is.finite(w["rt"]) && w["rt"] > 1) w["rt"] else NA_real_
        w_cv <- if (is.finite(w["cv"]) && w["cv"] > 1) w["cv"] else NA_real_

        if (is.na(w_rt) && r0 <= nrow(Z_width) && c0 <= ncol(Z_width)) {
          h_center <- Z_width[r0, c0]
          lo <- r0; prev_val <- h_center
          for (step in seq_len(max_radius_rt)) {
            rr <- r0 - step; if (rr < 1) break
            cur_val <- Z_width[rr, c0]
            if (cur_val <= 0 || cur_val > prev_val) break
            lo <- rr; prev_val <- cur_val
          }
          hi <- r0; prev_val <- h_center
          for (step in seq_len(max_radius_rt)) {
            rr <- r0 + step; if (rr > nrow(Z_width)) break
            cur_val <- Z_width[rr, c0]
            if (cur_val <= 0 || cur_val > prev_val) break
            hi <- rr; prev_val <- cur_val
          }
          w_rt <- hi - lo + 1
        }
        if (is.na(w_cv) && r0 <= nrow(Z_width) && c0 <= ncol(Z_width)) {
          h_center <- Z_width[r0, c0]
          lo <- c0; prev_val <- h_center
          for (step in seq_len(max_radius_cv)) {
            cc <- c0 - step; if (cc < 1) break
            cur_val <- Z_width[r0, cc]
            if (cur_val <= 0 || cur_val > prev_val) break
            lo <- cc; prev_val <- cur_val
          }
          hi <- c0; prev_val <- h_center
          for (step in seq_len(max_radius_cv)) {
            cc <- c0 + step; if (cc > ncol(Z_width)) break
            cur_val <- Z_width[r0, cc]
            if (cur_val <= 0 || cur_val > prev_val) break
            hi <- cc; prev_val <- cur_val
          }
          w_cv <- hi - lo + 1
        }
        if (is.na(w_rt)) { w_rt <- 1; n_width_fail <- n_width_fail + 1L }
        if (is.na(w_cv)) { w_cv <- 1; n_width_fail <- n_width_fail + 1L }
        if (is.finite(h0) && h0 > 0) {
          all_rt_loc <- c(all_rt_loc, r0 / R)
          all_cv_loc <- c(all_cv_loc, c0 / C)
          all_sigma_rt <- c(all_sigma_rt, w_rt * fwhm_to_sigma)
          all_sigma_cv <- c(all_sigma_cv, w_cv * fwhm_to_sigma)
          all_intensity <- c(all_intensity, h0)
        }
      }
    }
    # Noise estimation: mask out detected peak regions and use the
    # remaining pixels as the noise distribution. This is more
    # principled than the previous intensity-threshold approach, which
    # didn't distinguish between "dim peak shoulders" and "actual
    # background noise". With sufficient top_k (e.g. >= 500), the
    # detected peaks cover most of the real peak signal and the
    # unmasked region really is non-peak background.
    #
    # Mask half-extents: max_radius_* + a small buffer to capture peak
    # tails the bounded width-measurement walk may have missed.
    Z_noise <- Z_width
    if (nrow(centers) > 0) {
      mask_half_rt <- max_radius_rt + 2L
      mask_half_cv <- max_radius_cv + 1L
      peak_mask <- matrix(FALSE, nrow = nrow(Z_noise), ncol = ncol(Z_noise))
      for (j in seq_len(nrow(centers))) {
        r0 <- centers[j, 1]; c0 <- centers[j, 2]
        r_lo <- max(1L, r0 - mask_half_rt)
        r_hi <- min(nrow(Z_noise), r0 + mask_half_rt)
        c_lo <- max(1L, c0 - mask_half_cv)
        c_hi <- min(ncol(Z_noise), c0 + mask_half_cv)
        peak_mask[r_lo:r_hi, c_lo:c_hi] <- TRUE
      }
      non_peak <- Z_noise[!peak_mask]
      # Restrict to positive pixels (zeros are background regions of the
      # canvas that don't reflect instrument noise).
      bg <- non_peak[non_peak > 0]
    } else {
      # Fallback if no peaks were detected for this sample
      bg <- Z_noise[Z_noise > 0 & Z_noise <= eps_i]
    }
    if (length(bg) > 50) {
      all_noise_mean <- c(all_noise_mean, mean(bg))
      all_noise_sd <- c(all_noise_sd, stats::sd(bg))
    }
  }

  if (verbose) {
    message("  [estimate_peak_params] width measurement failed for ",
            n_width_fail, " / ", length(all_intensity),
            " peak-center checks (used fallback width=1)")
    if (length(all_sigma_rt) > 0) {
      message("  [estimate_peak_params] sigma_rt range: ",
              signif(min(all_sigma_rt), 3), " - ", signif(max(all_sigma_rt), 3),
              " (median ", signif(stats::median(all_sigma_rt), 3), ")")
      message("  [estimate_peak_params] sigma_cv range: ",
              signif(min(all_sigma_cv), 3), " - ", signif(max(all_sigma_cv), 3),
              " (median ", signif(stats::median(all_sigma_cv), 3), ")")
    }
  }

  fit_lognormal <- function(x) {
    x <- x[is.finite(x) & x > 0]
    if (length(x) < 10) return(list(meanlog = 0, sdlog = 1))
    list(meanlog = mean(log(x)), sdlog = stats::sd(log(x)))
  }
  safe_mean <- function(x) if (length(x) == 0) NA_real_ else mean(x)
  safe_sd   <- function(x) if (length(x) < 2) 0 else stats::sd(x)

  # Robust default for noise stats when no background pixels were
  # collected (e.g., aggressive dust thresholding zeros out the
  # sub-eps range). Use 0.1% of the median peak intensity as a small
  # but non-trivial fallback so synthetic generation does not produce
  # NaNs downstream.
  noise_mean_est <- safe_mean(all_noise_mean)
  noise_sd_est   <- safe_mean(all_noise_sd)
  if (!is.finite(noise_mean_est) || !is.finite(noise_sd_est)) {
    fallback <- if (length(all_intensity) > 0)
      stats::median(all_intensity) * 0.001 else 1e-4
    if (!is.finite(noise_mean_est)) noise_mean_est <- fallback
    if (!is.finite(noise_sd_est))   noise_sd_est   <- fallback / 2
    if (verbose) message("  [estimate_peak_params] no background pixels ",
                          "found; using fallback noise mean=",
                          signif(noise_mean_est, 3), ", sd=",
                          signif(noise_sd_est, 3))
  }

  result <- list(
    n_peaks = list(values = all_n_peaks, mean = safe_mean(all_n_peaks),
                   sd = safe_sd(all_n_peaks)),
    rt_loc = list(values = all_rt_loc, mean = safe_mean(all_rt_loc),
                  sd = safe_sd(all_rt_loc)),
    cv_loc = list(values = all_cv_loc, mean = safe_mean(all_cv_loc),
                  sd = safe_sd(all_cv_loc)),
    sigma_rt = fit_lognormal(all_sigma_rt),
    sigma_cv = fit_lognormal(all_sigma_cv),
    intensity = fit_lognormal(all_intensity),
    noise = list(mean = noise_mean_est, sd = noise_sd_est),
    n_samples = length(Z_list),
    n_peaks_detected = length(all_intensity),
    sigma_rt_raw = all_sigma_rt,
    sigma_cv_raw = all_sigma_cv,
    intensity_raw = all_intensity,
    rt_loc_raw = all_rt_loc,
    cv_loc_raw = all_cv_loc
  )
  if (verbose && result$n_peaks_detected == 0) {
    warning("estimate_peak_params: 0 peaks detected.")
  }
  result
}

#' Generate one synthetic GC-DMS image
#'
#' Generates a single synthetic image as a sum of 2-D Gaussian peaks
#' (sampled from the empirical peak distributions in `params`) plus
#' optional background noise.
#'
#' @param params Output of [estimate_peak_params()].
#' @param H,W Output image dimensions.
#' @param add_noise If TRUE, add spatially-varying background noise.
#' @param size_jitter Per-peak random log-normal scale factor SD for size
#'   diversity (default 0.6).
#' @param location_mode How to sample peak (RT, CV) locations. `"empirical"`
#'   (default, recommended) draws each peak's location by picking a random
#'   observed peak location from the real cohort (from `params$rt_loc_raw`
#'   / `params$cv_loc_raw`) and adding small per-peak jitter. This produces
#'   synthetic images whose peaks cluster at cohort-typical hotspots
#'   rather than scattering uniformly. `"marginal"` uses independent
#'   Normal draws on each axis with mean/SD from `params$rt_loc$mean`/`sd`
#'   — the previous behavior, kept for backward compatibility.
#' @param location_jitter_rt,location_jitter_cv SD of per-peak location
#'   jitter (in pixels) added on top of the empirical draw. Represents
#'   instrument reproducibility. Default 2 pixels RT, 1 pixel CV.
#' @return List with `clean` and `noisy` matrices (`clean == noisy` if
#'   `add_noise = FALSE`).
#' @export
generate_one_synthetic <- function(params, H, W, add_noise = TRUE,
                                    size_jitter = 0.6,
                                    location_mode = c("empirical", "marginal"),
                                    location_jitter_rt = 2,
                                    location_jitter_cv = 1) {
  location_mode <- match.arg(location_mode)

  if (params$n_peaks_detected == 0) {
    stop("generate_one_synthetic: params contain 0 detected peaks.")
  }

  # Empirical location mode requires stored raw observations.
  if (location_mode == "empirical") {
    if (is.null(params$rt_loc_raw) || is.null(params$cv_loc_raw) ||
        length(params$rt_loc_raw) == 0L) {
      warning("location_mode = 'empirical' requested but params lacks ",
              "rt_loc_raw / cv_loc_raw; falling back to 'marginal'.")
      location_mode <- "marginal"
    }
  }

  n_peaks <- max(1L, round(stats::rnorm(1,
    mean = params$n_peaks$mean,
    sd = max(params$n_peaks$sd, params$n_peaks$mean * 0.3))))
  clean <- matrix(0, nrow = H, ncol = W)
  rows <- seq_len(H); cols <- seq_len(W)

  # Pre-sample peak-location indices for empirical mode (one draw per
  # peak, with replacement). Adding per-peak jitter reflects instrument
  # reproducibility around each hotspot.
  if (location_mode == "empirical") {
    n_obs <- length(params$rt_loc_raw)
    loc_idx <- sample.int(n_obs, size = n_peaks, replace = TRUE)
  }

  for (k in seq_len(n_peaks)) {
    if (location_mode == "empirical") {
      j <- loc_idx[k]
      mu_rt <- params$rt_loc_raw[j] * H +
                 stats::rnorm(1, mean = 0, sd = location_jitter_rt)
      mu_cv <- params$cv_loc_raw[j] * W +
                 stats::rnorm(1, mean = 0, sd = location_jitter_cv)
    } else {
      mu_rt <- stats::rnorm(1, mean = params$rt_loc$mean * H,
                                sd = params$rt_loc$sd * H)
      mu_cv <- stats::rnorm(1, mean = params$cv_loc$mean * W,
                                sd = params$cv_loc$sd * W)
    }
    mu_rt <- max(1, min(H, mu_rt))
    mu_cv <- max(1, min(W, mu_cv))
    sig_rt <- stats::rlnorm(1, meanlog = params$sigma_rt$meanlog,
                              sdlog = params$sigma_rt$sdlog)
    sig_cv <- stats::rlnorm(1, meanlog = params$sigma_cv$meanlog,
                              sdlog = params$sigma_cv$sdlog)
    scale_factor <- exp(stats::rnorm(1, mean = 0, sd = size_jitter))
    sig_rt <- max(sig_rt * scale_factor, 0.8)
    sig_cv <- max(sig_cv * scale_factor, 0.8)
    amp <- stats::rlnorm(1, meanlog = params$intensity$meanlog,
                              sdlog = params$intensity$sdlog) *
      (0.5 + 0.5 * scale_factor)
    g_rt <- exp(-0.5 * ((rows - mu_rt) / sig_rt)^2)
    g_cv <- exp(-0.5 * ((cols - mu_cv) / sig_cv)^2)
    clean <- clean + amp * outer(g_rt, g_cv)
  }
  if (!add_noise) return(list(clean = clean, noisy = clean))
  # Belt-and-suspenders: if noise stats are somehow still non-finite,
  # skip noise rather than produce NaN samples that would crash training.
  if (!is.finite(params$noise$mean) || !is.finite(params$noise$sd)) {
    return(list(clean = clean, noisy = clean))
  }
  noise_sd <- max(1e-8, params$noise$sd)
  noise <- matrix(abs(stats::rnorm(H * W, mean = params$noise$mean,
                                    sd = noise_sd)),
                  nrow = H, ncol = W)
  list(clean = clean, noisy = clean + noise)
}

#' Generate N synthetic pairs as torch tensors
#'
#' Calls [generate_one_synthetic()] N times, applies per-sample
#' log-quantile normalization, and stacks the result as torch tensors
#' of shape `(N, 1, H, W)`.
#'
#' @inheritParams generate_one_synthetic
#' @param N Number of synthetic samples to generate.
#' @param normalize Whether to apply [normalize_sample()] (default TRUE).
#' @param dust_threshold If > 0, zero out pixels below this value in both
#'   the clean and noisy synthetic images BEFORE normalization. This
#'   matches the dust-thresholding step applied to real data via
#'   [baseline_basement()] and produces synthetic samples whose sparsity
#'   profile matches real samples. Default 0 (no dust thresholding,
#'   backward-compatible). To match real-data preprocessing, set this
#'   to the same `basement_thr` you passed to [baseline_basement()].
#' @return List with `clean` and `noisy` 4-D torch tensors.
#' @export
generate_synthetic_dataset <- function(params, N, H, W, add_noise = TRUE,
                                        size_jitter = 0.6,
                                        normalize = TRUE,
                                        dust_threshold = 0,
                                        location_mode = c("empirical",
                                                          "marginal"),
                                        location_jitter_rt = 2,
                                        location_jitter_cv = 1) {
  location_mode <- match.arg(location_mode)
  clean_arr <- array(0, dim = c(N, 1L, H, W))
  noisy_arr <- array(0, dim = c(N, 1L, H, W))
  for (i in seq_len(N)) {
    pair <- generate_one_synthetic(params, H, W,
                                    add_noise = add_noise,
                                    size_jitter = size_jitter,
                                    location_mode = location_mode,
                                    location_jitter_rt = location_jitter_rt,
                                    location_jitter_cv = location_jitter_cv)
    # Match real-data preprocessing: zero pixels below dust threshold
    # before normalization. Real data goes through baseline_basement at
    # this same point in the pipeline.
    if (dust_threshold > 0) {
      pair$clean[pair$clean < dust_threshold] <- 0
      pair$noisy[pair$noisy < dust_threshold] <- 0
    }
    if (normalize) {
      pair$clean <- normalize_sample(pair$clean)
      pair$noisy <- normalize_sample(pair$noisy)
    }
    clean_arr[i, 1, , ] <- pair$clean
    noisy_arr[i, 1, , ] <- pair$noisy
  }
  list(
    clean = torch::torch_tensor(clean_arr, dtype = torch::torch_float()),
    noisy = torch::torch_tensor(noisy_arr, dtype = torch::torch_float())
  )
}

#' Quantitative comparison of synthetic vs real GC-DMS samples
#'
#' Generates a batch of synthetic samples from `peak_params` and
#' compares their distributional summary statistics to those of real
#' samples. Use this to verify the synthetic generator's output
#' resembles real spectra before committing to a long pretraining run.
#'
#' Five comparisons are reported, each with a Kolmogorov-Smirnov D
#' statistic and p-value (small p-value indicates the real and
#' synthetic distributions differ):
#'
#' \enumerate{
#'   \item \strong{pixel_intensity_pooled}: distribution of all positive
#'     pixel values, pooled across samples. Tests whether the intensity
#'     dynamic range matches.
#'   \item \strong{total_signal}: per-sample sum of all pixels. Tests
#'     overall signal magnitude.
#'   \item \strong{mean_positive}: per-sample mean of positive pixels.
#'     Tests typical peak amplitude scale.
#'   \item \strong{sparsity}: per-sample fraction of zero pixels. Tests
#'     whether the proportion of background-vs-signal region matches.
#'   \item \strong{n_peaks}: number of detected peaks per sample (via
#'     [pick_peak_centers()]). Tests whether the peak density matches.
#' }
#'
#' @param Z_real_list List of real preprocessed Z matrices. Typically
#'   the same `Z_raw` you would pass to [estimate_peak_params()].
#' @param peak_params Output of [estimate_peak_params()].
#' @param H,W Dimensions for synthetic generation. Default uses the
#'   dimensions of the first real sample.
#' @param n_synthetic Number of synthetic samples to generate.
#' @param size_jitter Same as in [generate_one_synthetic()].
#' @param normalize If TRUE, apply [normalize_sample()] to both real
#'   and synthetic matrices before comparing. Default FALSE compares
#'   on the raw (post-baseline-correction, pre-normalization) scale.
#' @param top_k Maximum peaks counted per sample for the n_peaks
#'   comparison via [pick_peak_centers()]. Set high enough to NOT cap
#'   detected peak counts (e.g., 10000) so the comparison reflects
#'   actual peak density rather than the cap. Default 10000.
#' @param dust_threshold If > 0, apply the same dust threshold to
#'   synthetic samples that was applied to real data. Pass the same
#'   value used in [baseline_basement()] so synthetic and real are
#'   compared on the same preprocessed sparsity profile.
#' @param location_mode Passed through to [generate_one_synthetic()].
#'   Default `"empirical"` samples peak locations from the real
#'   cohort's observed hotspots.
#' @param location_jitter_rt,location_jitter_cv Passed through to
#'   [generate_one_synthetic()].
#' @param seed Optional RNG seed for reproducibility.
#' @return A list with:
#'   \item{summary}{data.frame of comparison metrics, one row per
#'     statistic.}
#'   \item{intensity_real}{pooled positive pixel intensities from real
#'     samples (numeric vector).}
#'   \item{intensity_synth}{same for synthetic samples.}
#'   \item{per_sample}{data.frame of per-sample stats for both groups,
#'     useful for plotting.}
#' @seealso [plot_synthetic_quality()] for a visual companion.
#' @export
synthetic_quality_check <- function(Z_real_list, peak_params,
                                      H = NULL, W = NULL,
                                      n_synthetic = 50L,
                                      size_jitter = 0.6,
                                      normalize = FALSE,
                                      top_k = 10000L,
                                      dust_threshold = 0,
                                      location_mode = c("empirical",
                                                        "marginal"),
                                      location_jitter_rt = 2,
                                      location_jitter_cv = 1,
                                      seed = NULL) {
  location_mode <- match.arg(location_mode)
  stopifnot(length(Z_real_list) > 0)
  if (is.null(H)) H <- nrow(as.matrix(Z_real_list[[1]]))
  if (is.null(W)) W <- ncol(as.matrix(Z_real_list[[1]]))
  if (!is.null(seed)) set.seed(seed)

  # Generate synthetic samples (noisy version, matching what pretraining sees)
  Z_synth_list <- vector("list", n_synthetic)
  for (i in seq_len(n_synthetic)) {
    pair <- generate_one_synthetic(peak_params, H, W,
                                    add_noise = TRUE,
                                    size_jitter = size_jitter,
                                    location_mode = location_mode,
                                    location_jitter_rt = location_jitter_rt,
                                    location_jitter_cv = location_jitter_cv)
    z <- pair$noisy
    if (dust_threshold > 0) z[z < dust_threshold] <- 0
    Z_synth_list[[i]] <- z
  }

  # Optional normalization to compare on the encoder-input scale
  if (normalize) {
    Z_real_use  <- lapply(Z_real_list,  normalize_sample)
    Z_synth_use <- lapply(Z_synth_list, normalize_sample)
  } else {
    Z_real_use  <- lapply(Z_real_list,  as.matrix)
    Z_synth_use <- Z_synth_list
  }

  # Per-sample summary stats
  summarize_sample <- function(Z) {
    pos <- Z[Z > 0]
    data.frame(
      total_signal  = sum(Z),
      mean_positive = if (length(pos) > 0) mean(pos) else 0,
      sparsity      = mean(Z == 0),
      n_peaks       = nrow(pick_peak_centers(Z, top_k = top_k))
    )
  }
  real_stats  <- do.call(rbind, lapply(Z_real_use,  summarize_sample))
  synth_stats <- do.call(rbind, lapply(Z_synth_use, summarize_sample))
  real_stats$group  <- "real"
  synth_stats$group <- "synthetic"
  per_sample <- rbind(real_stats, synth_stats)

  # Pooled positive-pixel intensities
  intensity_real  <- unlist(lapply(Z_real_use,  function(Z) Z[Z > 0]))
  intensity_synth <- unlist(lapply(Z_synth_use, function(Z) Z[Z > 0]))

  # KS test helper (returns NA's for too-small samples)
  ks_for <- function(x, y) {
    if (length(x) < 5 || length(y) < 5) return(list(D = NA_real_, p = NA_real_))
    res <- suppressWarnings(stats::ks.test(x, y))
    list(D = as.numeric(res$statistic), p = res$p.value)
  }

  # Build summary table
  stat_cols <- c("total_signal", "mean_positive", "sparsity", "n_peaks")
  rows <- lapply(stat_cols, function(s) {
    r <- real_stats[[s]]; sy <- synth_stats[[s]]
    ks <- ks_for(r, sy)
    data.frame(
      statistic   = s,
      real_mean   = mean(r),   real_sd  = stats::sd(r),
      synth_mean  = mean(sy),  synth_sd = stats::sd(sy),
      ks_D        = ks$D,      ks_p     = ks$p,
      stringsAsFactors = FALSE
    )
  })
  # Pooled-pixel-intensity comparison: subsample large vectors to keep
  # ks.test fast and avoid numerical issues with ties.
  ip_r <- if (length(intensity_real)  > 1e5) sample(intensity_real,  1e5) else intensity_real
  ip_s <- if (length(intensity_synth) > 1e5) sample(intensity_synth, 1e5) else intensity_synth
  ks_int <- ks_for(ip_r, ip_s)
  rows[[length(rows) + 1L]] <- data.frame(
    statistic = "pixel_intensity_pooled",
    real_mean = mean(intensity_real), real_sd = stats::sd(intensity_real),
    synth_mean = mean(intensity_synth), synth_sd = stats::sd(intensity_synth),
    ks_D = ks_int$D, ks_p = ks_int$p,
    stringsAsFactors = FALSE
  )

  list(
    summary         = do.call(rbind, rows),
    intensity_real  = intensity_real,
    intensity_synth = intensity_synth,
    per_sample      = per_sample
  )
}

#' Plot real vs synthetic per-sample summary statistics
#'
#' Faceted ggplot showing the distribution of each per-sample summary
#' statistic (total signal, mean positive intensity, sparsity, peak
#' count) for real and synthetic samples side-by-side. Use as a
#' visual companion to [synthetic_quality_check()].
#'
#' @param qc_result Output of [synthetic_quality_check()].
#' @return A ggplot object.
#' @export
plot_synthetic_quality <- function(qc_result) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting.")
  }
  if (!requireNamespace("tidyr", quietly = TRUE)) {
    stop("Package 'tidyr' is required (for pivoting per-sample stats).")
  }
  ps <- qc_result$per_sample
  long <- tidyr::pivot_longer(ps,
    cols = c("total_signal", "mean_positive", "sparsity", "n_peaks"),
    names_to = "statistic", values_to = "value")
  ggplot2::ggplot(long, ggplot2::aes(x = .data$value,
                                       fill = .data$group)) +
    ggplot2::geom_density(alpha = 0.45) +
    ggplot2::facet_wrap(~ statistic, scales = "free", ncol = 2) +
    ggplot2::scale_fill_manual(values = c("real"      = "steelblue",
                                            "synthetic" = "coral")) +
    ggplot2::labs(
      title = "Per-sample summary statistics: real vs synthetic",
      x = NULL, y = "density", fill = NULL
    ) +
    ggplot2::theme_minimal()
}
