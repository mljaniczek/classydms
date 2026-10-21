# classydms: mapping peak positions onto the generation / model canvas.
#
# estimate_peak_params() stores positions as fractions of each sample's
# UNPADDED matrix: rt_frac = r0 / R, cv_frac = c0 / C, with r0, c0 the
# 1-based row / column and R x C that sample's data dims. Real images are
# then centre-padded by pad_to_target() to H x W, which puts data row r0 at
# canvas row floor((H - R) / 2) + r0. Placing a peak at frac * H instead
# stretches the data block over the padding bands (about 1.31x in CV for
# 73 data columns padded to 96), so synthetic and real images disagree.
#
# Every generator and the presence head therefore place positions through
# .pad_frac_to_px(), using the per-sample data dims recorded on
# peak_params$data_dims and carried onto catalogs as
# catalog$parameters$data_dims plus the per-observation sample_idx_obs
# list-column. When the target canvas is smaller than the data grid along
# an axis (no padding possible), positions are scaled instead, which was
# the only behaviour before data dims were recorded.

.classydms_state <- new.env(parent = emptyenv())

# Map fraction positions to 1-based canvas pixels. `n_data` is the data
# extent along that axis (recycled to length(frac)); NA means unknown and
# falls back to scaling.
.pad_frac_to_px <- function(frac, n_data, n_target) {
  n_data <- rep_len(as.numeric(n_data), length(frac))
  out <- frac * n_target
  pad <- !is.na(n_data) & n_target >= n_data
  out[pad] <- floor((n_target - n_data[pad]) / 2) + frac[pad] * n_data[pad]
  out
}

.warn_scaled_placement <- function(fn) {
  key <- paste0("warned_", fn)
  if (isTRUE(.classydms_state[[key]])) return(invisible())
  assign(key, TRUE, envir = .classydms_state)
  warning(fn, ": no data dims recorded (catalog$parameters$data_dims / ",
          "peak_params$data_dims), so positions are placed at fraction * ",
          "H / W. That misregisters peaks against centre-padded real images ",
          "(pad_to_target). Rebuild the catalog with this classydms ",
          "version. Shown once per session.", call. = FALSE)
}

# Data dims of a peak_params object: the recorded `data_dims`, else the
# axis lengths in time_ranges / cv_ranges (estimate_peak_params requires
# axes at least as long as the matrices; trim / resample keep them equal).
.pp_data_dims <- function(peak_params, fn = "build_peak_catalog") {
  if (!is.null(peak_params$data_dims)) return(peak_params$data_dims)
  tr <- peak_params$time_ranges; cr <- peak_params$cv_ranges
  if (is.null(tr) || is.null(cr) || !"n" %in% colnames(tr) ||
      !"n" %in% colnames(cr)) return(NULL)
  if (!isTRUE(.classydms_state$inferred_dims_msg)) {
    assign("inferred_dims_msg", TRUE, envir = .classydms_state)
    message(fn, ": peak_params has no data_dims; using the axis lengths in ",
            "time_ranges / cv_ranges as the unpadded data dims (exact when ",
            "the axes passed to estimate_peak_params matched the ",
            "matrices). Shown once per session.")
  }
  cbind(rows = unname(tr[, "n"]), cols = unname(cr[, "n"]))
}

# Fraction columns holding a catalog compound's centroid position.
.catalog_frac_cols <- function(catalog) {
  cmp <- catalog$compounds
  has <- c("rt_frac", "cv_frac") %in% names(cmp)
  if ((catalog$parameters$coord_mode %||% "fraction") == "physical" ||
      any(has)) {
    if (!all(has))
      stop("physical-mode catalog must carry rt_frac and cv_frac columns. ",
           "Rebuild the catalog with a recent build_peak_catalog().",
           call. = FALSE)
    list(rt = cmp$rt_frac, cv = cmp$cv_frac)
  } else {
    list(rt = cmp$rt_loc, cv = cmp$cv_loc)
  }
}

# One-time placement setup for a catalog on an H x W canvas: centroid
# pixels (placed with the median data dims) and what per-observation
# placement needs. Generators call this once per sample, then .place_obs()
# per firing, which is O(1).
.catalog_placement <- function(catalog, compounds, H, W, fn) {
  dims <- catalog$parameters$data_dims
  if (is.null(dims)) .warn_scaled_placement(fn)
  has_obs <- all(c("rt_frac_obs", "cv_frac_obs") %in% names(compounds))
  # Per-observation placement needs each observation's sample; without
  # sample_idx_obs, observations fall back to scaling (centroids still
  # use the median data dims).
  obs_dims <- if ("sample_idx_obs" %in% names(compounds)) dims else NULL
  fr <- .catalog_frac_cols(list(compounds = compounds,
                                parameters = catalog$parameters))
  ref_r <- if (is.null(dims)) NA_real_ else stats::median(dims[, 1])
  ref_c <- if (is.null(dims)) NA_real_ else stats::median(dims[, 2])
  list(H = H, W = W, dims = obs_dims, has_obs = has_obs,
       ref_rows = ref_r, ref_cols = ref_c,
       rt_cent = .pad_frac_to_px(fr$rt, ref_r, H),
       cv_cent = .pad_frac_to_px(fr$cv, ref_c, W))
}

# Canvas-fraction range covering `frac_range` of the (median) data block,
# e.g. c(0.05, 0.95) of the data rather than of the padded canvas.
.data_band_range <- function(pl, frac_range, axis = c("rt", "cv")) {
  axis <- match.arg(axis)
  n_t <- if (axis == "rt") pl$H else pl$W
  n_d <- if (axis == "rt") pl$ref_rows else pl$ref_cols
  .pad_frac_to_px(frac_range, n_d, n_t) / n_t
}

# Canvas pixel (unjittered) for observation j of compound row i, falling
# back to the centroid when per-observation positions are unavailable.
.place_obs <- function(pl, compounds, i, j) {
  if (pl$has_obs && length(compounds$rt_frac_obs[[i]]) >= j) {
    rf <- compounds$rt_frac_obs[[i]][j]
    cf <- compounds$cv_frac_obs[[i]][j]
    if (is.finite(rf) && is.finite(cf)) {
      if (is.null(pl$dims))
        return(c(rf * pl$H, cf * pl$W))
      s <- compounds$sample_idx_obs[[i]][j]
      return(c(.pad_frac_to_px(rf, pl$dims[s, 1], pl$H),
               .pad_frac_to_px(cf, pl$dims[s, 2], pl$W)))
    }
  }
  c(pl$rt_cent[i], pl$cv_cent[i])
}

#' Catalog compound positions on a padded canvas
#'
#' Maps every compound's observations onto an `H x W` canvas laid out the
#' way [pad_to_target()] lays out real samples: each observation's
#' fraction position is converted with its own sample's unpadded data
#' dims (`catalog$parameters$data_dims`) and centre-padding offset. This
#' is the frame the synthetic generators render in and the frame the
#' encoder sees for real samples.
#'
#' Catalogs built before data dims were recorded fall back to
#' `fraction * H` (with a once-per-session warning).
#'
#' @param catalog A `peak_catalog`.
#' @param H,W Canvas dimensions in pixels.
#' @return A tibble with one row per compound: `compound_id`, `rt_px`,
#'   `cv_px` (median 1-based canvas pixel over the compound's
#'   observations, or the centroid when it has none), and `rt_px_sd`,
#'   `cv_px_sd` (SD of the observations' canvas pixels; 0 with fewer
#'   than two).
#' @export
catalog_pixel_positions <- function(catalog, H, W) {
  stopifnot(inherits(catalog, "peak_catalog"))
  cmp <- catalog$compounds
  pl <- .catalog_placement(catalog, cmp, H, W, "catalog_pixel_positions")
  n <- nrow(cmp)
  rt_px <- pl$rt_cent; cv_px <- pl$cv_cent
  rt_sd <- cv_sd <- numeric(n)
  if (pl$has_obs && n > 0L) {
    len <- lengths(cmp$rt_frac_obs)
    row <- rep(seq_len(n), len)
    rf <- unlist(cmp$rt_frac_obs, use.names = FALSE)
    cf <- unlist(cmp$cv_frac_obs, use.names = FALSE)
    if (is.null(pl$dims)) {
      rp <- rf * H; cp <- cf * W
    } else {
      s  <- unlist(cmp$sample_idx_obs, use.names = FALSE)
      rp <- .pad_frac_to_px(rf, pl$dims[s, 1], H)
      cp <- .pad_frac_to_px(cf, pl$dims[s, 2], W)
    }
    ok <- is.finite(rp) & is.finite(cp)
    rp_l <- split(rp[ok], factor(row[ok], levels = seq_len(n)))
    cp_l <- split(cp[ok], factor(row[ok], levels = seq_len(n)))
    has <- lengths(rp_l) > 0L
    med <- function(v) if (length(v)) stats::median(v) else NA_real_
    sdv <- function(v) if (length(v) > 1L) stats::sd(v) else 0
    rt_px[has] <- vapply(rp_l[has], med, numeric(1))
    cv_px[has] <- vapply(cp_l[has], med, numeric(1))
    rt_sd <- vapply(rp_l, sdv, numeric(1))
    cv_sd <- vapply(cp_l, sdv, numeric(1))
  }
  tibble::tibble(compound_id = cmp$compound_id,
                 rt_px = unname(rt_px), cv_px = unname(cv_px),
                 rt_px_sd = unname(rt_sd), cv_px_sd = unname(cv_sd))
}

# Canvas pixel positions for peak_params peaks `j` (vectorised).
.place_pp_peaks <- function(params, j, H, W) {
  dims <- params$data_dims
  if (is.null(dims) || is.null(params$sample_idx_raw)) {
    .warn_scaled_placement("generate_one_synthetic")
    return(cbind(params$rt_loc_raw[j] * H, params$cv_loc_raw[j] * W))
  }
  s <- params$sample_idx_raw[j]
  cbind(.pad_frac_to_px(params$rt_loc_raw[j], dims[s, 1], H),
        .pad_frac_to_px(params$cv_loc_raw[j], dims[s, 2], W))
}
