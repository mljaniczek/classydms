# classydms: interpretability and diagnostic functions

#' Reconstruct real images through a pre-trained autoencoder
#'
#' Feeds each requested sample through the autoencoder and returns the
#' original, the reconstruction, and the per-sample MSE. Useful as a
#' qualitative check that synthetic pretraining produced an autoencoder
#' that generalizes to real data.
#'
#' @param autoencoder A trained `dms_denoising_autoencoder`.
#' @param Z_list List of preprocessed, padded intensity matrices.
#' @param sample_indices Which samples to reconstruct (default first 6).
#' @param device "cpu" or "cuda".
#' @return List of `list(index, original, reconstructed, mse)`.
#' @export
reconstruct_real_samples <- function(autoencoder, Z_list,
                                      sample_indices = 1:6,
                                      device = "cpu") {
  autoencoder$eval(); autoencoder$to(device = device)
  results <- list()
  for (idx in sample_indices) {
    Z <- as.matrix(Z_list[[idx]])
    x <- torch::torch_tensor(Z, dtype = torch::torch_float())$
      unsqueeze(1)$unsqueeze(1)$to(device = device)
    torch::with_no_grad({ recon <- autoencoder(x) })
    recon_mat <- as.matrix(recon$squeeze()$to(device = "cpu"))
    mse <- mean((Z - recon_mat)^2)
    results[[length(results) + 1L]] <- list(
      index = idx, original = Z, reconstructed = recon_mat, mse = mse)
  }
  results
}

#' Compute a Class Activation Map (CAM) for one sample
#'
#' Weights each encoder channel's spatial activation by its elastic net
#' coefficient (`coefs`), summing across channels to produce a 2-D
#' saliency map. Because the elastic net classifier head is linear in
#' the global-average-pooled features, this construction is mathematically
#' equivalent to Class Activation Mapping using elastic net coefficients
#' as head weights (and similar in spirit to Grad-CAM, without requiring
#' backpropagation).
#'
#' Returns the saliency at the encoder output resolution (raw_cam) and
#' bilinearly upsampled to the input dimensions (upsampled_cam).
#'
#' @param encoder A pre-trained `dms_encoder`.
#' @param Z One sample's preprocessed, padded matrix.
#' @param coefs Numeric vector of length 64 (elastic net coefficients).
#' @param device "cpu" or "cuda".
#' @return List with `raw_cam` and `upsampled_cam`.
#' @export
compute_saliency_map <- function(encoder, Z, coefs, device = "cpu") {
  encoder$eval(); encoder$to(device = device)
  x <- torch::torch_tensor(as.matrix(Z), dtype = torch::torch_float())$
    unsqueeze(1)$unsqueeze(1)$to(device = device)
  torch::with_no_grad({ fmap <- encoder(x) })
  fmap_arr <- as.array(fmap$squeeze(1)$to(device = "cpu"))
  raw_cam <- matrix(0, nrow = dim(fmap_arr)[2], ncol = dim(fmap_arr)[3])
  for (c in seq_len(64L)) {
    if (coefs[c] != 0) raw_cam <- raw_cam + coefs[c] * fmap_arr[c, , ]
  }
  H_in <- nrow(Z); W_in <- ncol(Z)
  cam_t <- torch::torch_tensor(raw_cam, dtype = torch::torch_float())$
    unsqueeze(1)$unsqueeze(1)
  cam_up <- torch::nnf_interpolate(cam_t, size = c(H_in, W_in),
                                    mode = "bilinear",
                                    align_corners = FALSE)
  list(raw_cam = raw_cam, upsampled_cam = as.matrix(cam_up$squeeze()))
}

#' Aggregate saliency maps across samples, masking padded regions
#'
#' Computes a single saliency map averaged over `sample_indices`. For
#' each sample, only the original (non-padded) region contributes; the
#' divisor at each pixel is the number of samples that actually had data
#' there. Prevents data/padding boundary artifacts from dominating the
#' aggregate map.
#'
#' @param encoder Pre-trained `dms_encoder`.
#' @param Z_list List of padded matrices.
#' @param coefs Elastic net coefficients (length 64).
#' @param sample_indices Which samples to aggregate over.
#' @param orig_dims List of `c(H, W)` pairs giving each sample's
#'   pre-padding dimensions.
#' @param device "cpu" or "cuda".
#' @return Matrix of the same shape as the padded image.
#' @export
aggregate_saliency_masked <- function(encoder, Z_list, coefs,
                                        sample_indices, orig_dims,
                                        device = "cpu") {
  H <- nrow(Z_list[[1]]); W <- ncol(Z_list[[1]])
  sum_mat <- matrix(0, H, W); count_mat <- matrix(0, H, W)
  for (i in sample_indices) {
    cam <- compute_saliency_map(encoder, Z_list[[i]], coefs,
                                  device = device)$upsampled_cam
    h0 <- orig_dims[[i]]["H"]; w0 <- orig_dims[[i]]["W"]
    mask <- matrix(0, H, W)
    mask[seq_len(h0), seq_len(w0)] <- 1
    sum_mat <- sum_mat + cam * mask
    count_mat <- count_mat + mask
  }
  count_mat[count_mat == 0] <- 1
  sum_mat / count_mat
}
