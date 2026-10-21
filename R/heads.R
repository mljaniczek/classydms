# classydms: spatially-aware readout heads that replace global average
# pooling (GAP). GAP averages the encoder feature map over all of RT x CV,
# which discards WHERE activation occurred — and in GC-DMS, location is
# compound identity. These heads keep it:
#   - dms_presence_head: one learned query per catalog compound, attending
#     around that compound's catalog location. Trained during pretraining
#     to predict which compounds are present; its per-compound logits are
#     then an interpretable feature vector for classification.
#   - attention_pool_readout: a few learned queries with 2D positional
#     encoding, trained end-to-end in the downstream classifier.

#' Catalog geometry for a compound-anchored head
#'
#' Collects each compound's position on the `H x W` canvas, expressed as
#' a fraction in the same convention as the head's token centres (pixel
#' `p` sits at `(p - 0.5) / H`), and a positional spread used as the
#' width of the head's location prior. Positions come from
#' [catalog_pixel_positions()], which maps every observation onto the
#' canvas the way [pad_to_target()] centre-pads real samples, so the
#' prior is centred where the peak sits in both synthetic and real
#' images. The spread is the SD of the observations' canvas positions
#' combined with any extra displacement the pretraining augmentation can
#' introduce (affine shift, RT warp), so the prior stays wide enough to
#' find a shifted peak. (The head additionally floors the width at one
#' feature-map cell.)
#'
#' @param catalog A `peak_catalog`.
#' @param H,W Image dimensions in pixels.
#' @param compound_ids Optional subset of `compound_id`s. Default: every
#'   compound with a finite position.
#' @param extra_rt_px,extra_cv_px Additional displacement (pixels) to fold
#'   into the prior width, e.g. the maximum affine shift.
#' @return List with `compound_id`, `rt_frac`, `cv_frac`, `sd_rt_frac`,
#'   `sd_cv_frac`, all as fractions of `H` / `W`.
#' @export
compound_head_geometry <- function(catalog, H, W, compound_ids = NULL,
                                    extra_rt_px = 0, extra_cv_px = 0) {
  stopifnot(inherits(catalog, "peak_catalog"))
  cmp <- catalog$compounds
  if (!is.null(compound_ids)) {
    missing <- setdiff(compound_ids, cmp$compound_id)
    if (length(missing))
      stop("compound_head_geometry: compound_id(s) not in catalog: ",
           paste(utils::head(missing, 5), collapse = ", "))
  }
  pos <- catalog_pixel_positions(catalog, H, W)
  rt <- (pos$rt_px - 0.5) / H
  cv <- (pos$cv_px - 0.5) / W
  keep <- is.finite(rt) & is.finite(cv)
  if (!is.null(compound_ids)) keep <- keep & cmp$compound_id %in% compound_ids
  if (!any(keep)) stop("compound_head_geometry: no compounds with finite positions.")
  sd_rt <- sqrt((pos$rt_px_sd / H)^2 + (extra_rt_px / H)^2)
  sd_cv <- sqrt((pos$cv_px_sd / W)^2 + (extra_cv_px / W)^2)
  list(compound_id = cmp$compound_id[keep],
       rt_frac = rt[keep], cv_frac = cv[keep],
       sd_rt_frac = sd_rt[keep], sd_cv_frac = sd_cv[keep])
}

#' Rank-based (Mann-Whitney) AUC; NA when only one class is present.
#' @keywords internal
rank_auc <- function(y, s) {
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  if (n1 == 0L || n0 == 0L) return(NA_real_)
  r <- rank(s)
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

#' 2D sinusoidal positional encoding
#'
#' Returns a `(d, h, w)` tensor: the first half of the channels encode
#' the row (RT) position, the second half the column (CV) position, using
#' the standard transformer sin/cos frequency ladder.
#' @keywords internal
positional_encoding_2d <- function(h, w, d, device = "cpu") {
  stopifnot(d %% 4L == 0L)
  half <- d %/% 2L
  enc_1d <- function(n) {
    pos <- seq_len(n) - 1
    freqs <- 1 / (10000^(seq(0, half - 2, by = 2) / half))
    ang <- outer(pos, freqs)
    m <- matrix(0, n, half)
    m[, seq(1, half, by = 2)] <- sin(ang)
    m[, seq(2, half, by = 2)] <- cos(ang)
    m
  }
  pe_r <- enc_1d(h); pe_c <- enc_1d(w)
  arr <- array(0, dim = c(d, h, w))
  for (k in seq_len(half)) {
    arr[k, , ] <- matrix(pe_r[, k], h, w)
    arr[half + k, , ] <- matrix(pe_c[, k], h, w, byrow = TRUE)
  }
  torch::torch_tensor(arr, dtype = torch::torch_float(), device = device)
}

#' Compound-anchored attention pooling
#'
#' One learned query per catalog compound. Attention logits are
#' `q_k . key_t / sqrt(d) + log_prior_k(t)`, where the log-prior is a
#' Gaussian centred on compound k's catalog location with width taken
#' from [compound_head_geometry()] (floored at one feature-map cell, and
#' scaled by a single learnable factor). The prior tells each query WHERE
#' to look; the content term picks the peak within that window. Values
#' are the raw encoder features, so each compound's output token lives in
#' the encoder's channel space.
#'
#' @param in_channels Encoder output channels.
#' @param geometry Output of [compound_head_geometry()].
#' @param d_attn Query/key dimension.
#' @return Module whose forward returns `list(tokens = (B, K, C),
#'   attn = (B, K, h*w))`.
#' @export
compound_attention_pool <- torch::nn_module(
  "compound_attention_pool",
  initialize = function(in_channels, geometry, d_attn = 32L) {
    self$compound_ids <- geometry$compound_id
    self$rt_frac <- as.numeric(geometry$rt_frac)
    self$cv_frac <- as.numeric(geometry$cv_frac)
    self$sd_rt_frac <- as.numeric(geometry$sd_rt_frac)
    self$sd_cv_frac <- as.numeric(geometry$sd_cv_frac)
    self$d_attn <- as.integer(d_attn)
    K <- length(self$compound_ids)
    self$key <- torch::nn_conv2d(in_channels, d_attn, kernel_size = 1L)
    # Shared query + per-compound offset (zero-init). "Find the peak in
    # my window" is the same job for every compound — the location prior
    # supplies identity — so every sample trains the shared query, and
    # offsets only specialise where the data supports it. Independent
    # per-compound queries would each see only their own rare positives.
    self$query_shared <- torch::nn_parameter(torch::torch_randn(1, d_attn) * 0.02)
    self$query_offset <- torch::nn_parameter(torch::torch_zeros(K, d_attn))
    self$log_width <- torch::nn_parameter(torch::torch_zeros(1))
  },
  log_prior = function(h, w, device) {
    tt <- function(x) torch::torch_tensor(x, dtype = torch::torch_float(),
                                            device = device)
    r_tok <- tt((seq_len(h) - 0.5) / h)
    c_tok <- tt((seq_len(w) - 0.5) / w)
    scale <- torch::torch_exp(self$log_width)
    sd_rt <- torch::torch_clamp(tt(self$sd_rt_frac), min = 1 / h) * scale
    sd_cv <- torch::torch_clamp(tt(self$sd_cv_frac), min = 1 / w) * scale
    dr <- (r_tok$unsqueeze(1) - tt(self$rt_frac)$unsqueeze(2)) / sd_rt$unsqueeze(2)
    dc <- (c_tok$unsqueeze(1) - tt(self$cv_frac)$unsqueeze(2)) / sd_cv$unsqueeze(2)
    lp <- -0.5 * dr$pow(2)$unsqueeze(3) - 0.5 * dc$pow(2)$unsqueeze(2)
    lp$flatten(start_dim = 2)
  },
  forward = function(fmap) {
    h <- fmap$size(3); w <- fmap$size(4)
    keys <- self$key(fmap)$flatten(start_dim = 3)
    query <- self$query_shared + self$query_offset
    logits <- torch::torch_matmul(query, keys) / sqrt(self$d_attn)
    logits <- logits + self$log_prior(h, w, fmap$device)$unsqueeze(1)
    attn <- torch::nnf_softmax(logits, dim = 3)
    vals <- fmap$flatten(start_dim = 3)$transpose(2, 3)
    list(tokens = torch::torch_matmul(attn, vals), attn = attn)
  }
)

#' Compound-presence head
#'
#' Compound-anchored attention pooling ([compound_attention_pool()])
#' followed by a per-compound linear readout, producing one presence
#' logit per catalog compound. Trained as an auxiliary pretraining task in
#' [pretrain_denoising_online()] (`presence_head = TRUE`) against the
#' generator's ground-truth list of compounds it rendered. Downstream, the
#' logits form an `(N, K)` feature matrix whose columns ARE catalog
#' compounds — see [extract_encoder_features()] with
#' `readout = "presence"`.
#'
#' @inheritParams compound_attention_pool
#' @export
dms_presence_head <- torch::nn_module(
  "dms_presence_head",
  initialize = function(in_channels, geometry, d_attn = 32L) {
    self$compound_ids <- geometry$compound_id
    self$in_channels <- as.integer(in_channels)
    K <- length(self$compound_ids)
    self$pool <- compound_attention_pool(in_channels, geometry, d_attn)
    # Shared peak detector + zero-init per-compound offsets (same
    # rationale as the shared query); per-compound bias absorbs each
    # compound's base rate.
    self$weight_shared <- torch::nn_parameter(
      torch::torch_randn(1, in_channels) / sqrt(in_channels))
    self$weight_offset <- torch::nn_parameter(torch::torch_zeros(K, in_channels))
    self$bias <- torch::nn_parameter(torch::torch_zeros(K))
    # Running firing counts from pretraining, used for rarity-weighted
    # BCE. Plain R fields so they persist through torch_save().
    self$fire_count <- numeric(K)
    self$n_seen <- 0
  },
  forward = function(fmap) {
    tokens <- self$pool(fmap)$tokens
    w <- self$weight_shared + self$weight_offset
    (tokens * w$unsqueeze(1))$sum(dim = 3) + self$bias
  }
)

#' Learned-query attention pooling readout
#'
#' Replacement for GAP in a downstream classifier: `n_queries` learned
#' queries attend over the encoder feature map (keys carry a fixed 2D
#' sinusoidal positional encoding, so queries can become location-
#' specific). Output is the concatenation of the `n_queries` pooled
#' feature vectors, `(B, n_queries * C)`. Trained end-to-end with the
#' classifier; keep `n_queries` small when labeled n is small.
#'
#' @param in_channels Encoder output channels.
#' @param n_queries Number of learned queries.
#' @param d_attn Query/key dimension (multiple of 4).
#' @export
attention_pool_readout <- torch::nn_module(
  "attention_pool_readout",
  initialize = function(in_channels, n_queries = 4L, d_attn = 32L) {
    self$n_queries <- as.integer(n_queries)
    self$d_attn <- as.integer(d_attn)
    self$out_dim <- as.integer(n_queries * in_channels)
    self$key <- torch::nn_conv2d(in_channels, d_attn, kernel_size = 1L)
    self$query <- torch::nn_parameter(torch::torch_randn(n_queries, d_attn) * 0.02)
  },
  forward = function(fmap) {
    h <- fmap$size(3); w <- fmap$size(4)
    pe <- positional_encoding_2d(h, w, self$d_attn, device = fmap$device)
    keys <- (self$key(fmap) + pe$unsqueeze(1))$flatten(start_dim = 3)
    attn <- torch::nnf_softmax(
      torch::torch_matmul(self$query, keys) / sqrt(self$d_attn), dim = 3)
    vals <- fmap$flatten(start_dim = 3)$transpose(2, 3)
    torch::torch_matmul(attn, vals)$flatten(start_dim = 2)
  }
)
