# classydms: classification via pre-trained encoder
# Both neural-network-head (cv_pretrained_oof) and elastic-net-head
# (cv_encoder_elastic_net) classifiers, plus the feature extraction
# helper that bridges encoder output to the elastic net.

#' Attach a classifier head to a pre-trained encoder
#'
#' Deep-copies the encoder so each fold gets its own independent weights
#' (prevents cross-fold leakage when `unfreeze_last_block = TRUE`), wraps
#' it with adaptive average pooling and a dropout + linear classifier head.
#'
#' @param encoder A pre-trained `dms_encoder`.
#' @param num_classes Number of output classes (default 2).
#' @param dropout_p Dropout probability before the linear layer.
#' @param freeze_backbone If TRUE, encoder weights are not updated.
#' @param unfreeze_last_block If TRUE, the encoder's last residual block
#'   (layer3) is set trainable even if `freeze_backbone = TRUE`.
#' @return An instantiated torch module.
#' @export
build_pretrained_classifier <- function(encoder,
                                         num_classes = 2L,
                                         dropout_p = 0.3,
                                         freeze_backbone = TRUE,
                                         unfreeze_last_block = TRUE) {
  enc <- encoder$clone(deep = TRUE)
  # Auto-detect latent dim so encoders built with wider
  # encoder_channels work without needing an explicit latent_dim arg.
  # Peek at the encoder_channels attr if present (recent dms_encoder),
  # otherwise fall back to 64 (legacy encoders).
  latent_dim <- if (!is.null(enc$encoder_channels))
                  as.integer(tail(enc$encoder_channels, 1L)) else 64L
  classifier <- torch::nn_module(
    initialize = function() {
      self$encoder <- enc
      self$pool <- torch::nn_adaptive_avg_pool2d(output_size = c(1L, 1L))
      self$head <- torch::nn_sequential(
        torch::nn_flatten(),
        torch::nn_dropout(p = dropout_p),
        torch::nn_linear(latent_dim, num_classes)
      )
      if (freeze_backbone) {
        for (p in self$encoder$parameters) p$requires_grad_(FALSE)
      }
      if (unfreeze_last_block && !is.null(self$encoder$layer3)) {
        for (p in self$encoder$layer3$parameters) p$requires_grad_(TRUE)
      }
      for (p in self$head$parameters) p$requires_grad_(TRUE)
    },
    forward = function(x) {
      feat <- self$encoder(x)
      pooled <- self$pool(feat)
      self$head(pooled)
    }
  )
  classifier()
}

#' Train one classifier fold with a pre-trained encoder
#' @keywords internal
#' @export
train_one_fold_pretrained <- function(ds, train_idx, test_idx, encoder,
                                       epochs = 20L, batch_size = 4L,
                                       lr = 1e-3, weight_decay = 1e-4,
                                       dropout_p = 0.3,
                                       freeze_backbone = TRUE,
                                       unfreeze_last_block = TRUE,
                                       augment = TRUE,
                                       rt_shift = 25L, cv_shift = 5L,
                                       intensity_scale = c(0.9, 1.1),
                                       noise_sd = 0.0,
                                       device = if (torch::cuda_is_available()) "cuda" else "cpu",
                                       seed = 1L) {
  set.seed(seed); torch::torch_manual_seed(seed)
  ds_train <- subset_dataset_view(ds, train_idx)
  ds_test  <- subset_dataset_view(ds, test_idx)
  train_dl <- torch::dataloader(ds_train, batch_size = batch_size, shuffle = TRUE)
  test_dl  <- torch::dataloader(ds_test,  batch_size = batch_size, shuffle = FALSE)

  model <- build_pretrained_classifier(encoder, num_classes = 2L,
                                        dropout_p = dropout_p,
                                        freeze_backbone = freeze_backbone,
                                        unfreeze_last_block = unfreeze_last_block)
  model$to(device = device)
  params <- purrr::keep(model$parameters, ~.x$requires_grad)
  opt <- torch::optim_adam(params, lr = lr, weight_decay = weight_decay)
  loss_fn <- torch::nn_cross_entropy_loss()

  model$train()
  for (ep in seq_len(epochs)) {
    coro::loop(for (b in train_dl) {
      x <- b$x$to(device = device); y <- b$y$to(device = device)
      if (augment) {
        x <- augment_batch(x, rt_shift = rt_shift, cv_shift = cv_shift,
                            intensity_scale = intensity_scale,
                            noise_sd = noise_sd, device = device)
      }
      opt$zero_grad()
      logits <- model(x)
      loss <- loss_fn(logits, y)
      loss$backward(); opt$step()
    })
  }

  model$eval()
  probs <- numeric(length(test_idx))
  cursor <- 1L
  coro::loop(for (b in test_dl) {
    x <- b$x$to(device = device)
    logits <- model(x)
    p <- torch::nnf_softmax(logits, dim = 2)[, 2]
    p <- as.numeric(p$to(device = "cpu"))
    n <- length(p)
    probs[cursor:(cursor + n - 1L)] <- p
    cursor <- cursor + n
  })
  list(probs = probs)
}

#' 5-fold CV classification with a pre-trained encoder
#'
#' Out-of-fold predictions are aggregated into a single AUC. Each fold
#' independently deep-clones the encoder via
#' [build_pretrained_classifier()] to prevent cross-fold leakage.
#'
#' @param Z_list List of preprocessed, padded intensity matrices.
#' @param y Binary labels (integer or factor, two unique values).
#' @param encoder A pre-trained `dms_encoder`.
#' @param k Number of CV folds (default 5).
#' @param epochs Number of training epochs per fold.
#' @param batch_size Mini-batch size.
#' @param lr Adam learning rate.
#' @param weight_decay Adam L2 weight decay.
#' @param dropout_p Dropout probability before the linear classifier.
#' @param freeze_backbone If TRUE, encoder weights are not updated
#'   during fine-tuning.
#' @param unfreeze_last_block If TRUE, the last encoder residual block
#'   (`layer3`) is set trainable even when `freeze_backbone = TRUE`.
#' @param augment If TRUE, apply training-time augmentation
#'   ([augment_batch()]: random RT/CV shifts, intensity scaling, noise).
#' @param rt_shift,cv_shift Maximum augmentation shift magnitudes
#'   (pixels) along RT and CV.
#' @param intensity_scale Length-2 numeric vector giving the bounds of
#'   random intensity rescaling.
#' @param noise_sd SD of Gaussian noise added to augmented images
#'   (0 disables).
#' @param seed RNG seed for fold assignment and training.
#' @return List with `probs` (out-of-fold predictions) and `auc`.
#' @export
cv_pretrained_oof <- function(Z_list, y, encoder,
                                k = 5L, epochs = 20L, batch_size = 4L,
                                lr = 1e-3, weight_decay = 1e-4,
                                dropout_p = 0.3, freeze_backbone = TRUE,
                                unfreeze_last_block = TRUE,
                                augment = TRUE,
                                rt_shift = 25L, cv_shift = 5L,
                                intensity_scale = c(0.9, 1.1),
                                noise_sd = 0.0, seed = 1L) {
  ds <- build_tensor_dataset(Z_list, y)
  folds <- make_stratified_folds(y, k = k)
  oof <- rep(NA_real_, length(y))
  for (j in seq_len(k)) {
    message("Fold ", j, "/", k)
    tr <- folds[[j]]$train; te <- folds[[j]]$test
    fit <- train_one_fold_pretrained(
      ds = ds, train_idx = tr, test_idx = te, encoder = encoder,
      epochs = epochs, batch_size = batch_size,
      lr = lr, weight_decay = weight_decay,
      dropout_p = dropout_p,
      freeze_backbone = freeze_backbone,
      unfreeze_last_block = unfreeze_last_block,
      augment = augment, rt_shift = rt_shift, cv_shift = cv_shift,
      intensity_scale = intensity_scale, noise_sd = noise_sd,
      seed = 100L + j
    )
    oof[te] <- fit$probs
  }
  list(probs = oof, auc = auc_from_probs(y, oof))
}

#' Extract encoder features for all samples (GAP-pooled)
#'
#' Runs each sample through the (frozen) encoder, applies global average
#' pooling to its feature map, and stacks the result into an
#' `(N, latent_dim)` matrix suitable for use as input to a classical
#' classifier. Latent dim is auto-detected from the encoder (either via
#' its `encoder_channels` attribute or by peeking at its output shape
#' on the first sample), so encoders built with wider
#' `encoder_channels` work without additional args.
#'
#' @param encoder A pre-trained `dms_encoder`.
#' @param Z_list List of preprocessed, padded intensity matrices.
#' @param device "cpu" or "cuda".
#' @return Numeric matrix with one row per sample, `latent_dim` columns.
#' @export
extract_encoder_features <- function(encoder, Z_list, device = "cpu") {
  encoder$eval(); encoder$to(device = device)
  gap <- torch::nn_adaptive_avg_pool2d(c(1L, 1L))
  n <- length(Z_list)
  # Determine latent_dim: prefer the encoder_channels attribute (recent
  # encoders record it); otherwise probe the first sample's output.
  latent_dim <- if (!is.null(encoder$encoder_channels))
                  as.integer(tail(encoder$encoder_channels, 1L))
                else {
                  Z0 <- as.matrix(Z_list[[1]])
                  x0 <- torch::torch_tensor(Z0, dtype = torch::torch_float())$
                    unsqueeze(1)$unsqueeze(1)$to(device = device)
                  torch::with_no_grad({
                    dim <- as.integer(encoder(x0)$size()[2])
                  })
                  dim
                }
  feat_mat <- matrix(NA_real_, nrow = n, ncol = latent_dim)
  for (i in seq_len(n)) {
    Z <- as.matrix(Z_list[[i]])
    x <- torch::torch_tensor(Z, dtype = torch::torch_float())$
      unsqueeze(1)$unsqueeze(1)$to(device = device)
    torch::with_no_grad({
      fmap <- encoder(x)
      pooled <- gap(fmap)
      vec <- pooled$squeeze()
    })
    feat_mat[i, ] <- as.numeric(vec$to(device = "cpu"))
  }
  feat_mat
}

#' Elastic net classification on 64-dim encoder features
#'
#' Performs stratified k-fold CV: extract features once via the frozen
#' encoder, then fit `cv.glmnet` independently on each fold. Returns
#' out-of-fold predictions, AUC, the feature matrix, and elastic net
#' coefficients (which map directly to encoder channels and are the
#' starting point for spatial saliency interpretation).
#'
#' @param encoder A pre-trained `dms_encoder`.
#' @param Z_list List of preprocessed, padded intensity matrices.
#' @param y Binary labels.
#' @param k Number of CV folds.
#' @param alpha Elastic net mixing (1 = LASSO, 0 = ridge).
#' @param nlambda Number of lambda values for `cv.glmnet`.
#' @param seed RNG seed for fold assignment.
#' @param device "cpu" or "cuda".
#' @return List with `probs`, `auc`, `features`, and `coefs`.
#' @export
cv_encoder_elastic_net <- function(encoder, Z_list, y,
                                     k = 5L, alpha = 0.5,
                                     nlambda = 100L, seed = 1L,
                                     device = "cpu") {
  feat_mat <- extract_encoder_features(encoder, Z_list, device = device)
  folds <- make_stratified_folds(y, k = k, seed = seed)
  oof <- rep(NA_real_, length(y))
  coef_list <- list()
  for (j in seq_len(k)) {
    set.seed(seed + j)
    train_idx <- folds[[j]]$train; test_idx <- folds[[j]]$test
    X_train <- feat_mat[train_idx, , drop = FALSE]
    X_test  <- feat_mat[test_idx,  , drop = FALSE]
    y_train <- y[train_idx]
    fit <- glmnet::cv.glmnet(X_train, y_train, family = "binomial",
                              alpha = alpha, nlambda = nlambda)
    preds <- as.numeric(stats::predict(fit, X_test, s = "lambda.min",
                                         type = "response"))
    oof[test_idx] <- preds
    coef_list[[j]] <- as.numeric(stats::coef(fit, s = "lambda.min"))[-1]
  }
  avg_coefs <- Reduce("+", coef_list) / length(coef_list)
  names(avg_coefs) <- paste0("encoder_ch_", seq_len(64))
  list(probs = oof, auc = auc_from_probs(y, oof),
       features = feat_mat, coefs = avg_coefs)
}
