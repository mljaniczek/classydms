# classydms: denoising autoencoder architecture
# Encoder + decoder + composed denoising autoencoder. The encoder is the
# part that is transferred to downstream classification; the decoder is
# discarded after pretraining (except for reconstruction QC).

#' GC-DMS encoder (ResNet-style, asymmetric stride)
#'
#' Three residual block groups with an asymmetric stem stride that
#' aggressively downsamples RT while preserving CV resolution. Output
#' is a `(B, encoder_channels[5], h, w)` feature map suitable as input
#' to a classification head or for spatial saliency mapping.
#'
#' Channel schedule via `encoder_channels`, a length-5 integer vector
#' `c(stem_ch, layer1_ch, layer2_ch, layer3_ch, latent_ch)`:
#' - Stem: `1 -> stem_ch` (kernel 7, stride `(stem_stride_rt, stem_stride_cv)`)
#' - MaxPool: stride `(2, 1)`
#' - Layer1: `stem_ch -> layer1_ch`, 2 residual blocks at stride 1
#' - Layer2: `layer1_ch -> layer2_ch`, 2 blocks with first stride 2
#' - Layer3: `layer2_ch -> layer3_ch`, 2 blocks with first stride 2
#' - Optional 1x1 projection: `layer3_ch -> latent_ch` (skipped when equal)
#'
#' Total spatial downsample is fixed at 32x on RT and 4x on CV
#' regardless of channel counts (stem 4x1, maxpool 2x1, layer2 2x2,
#' layer3 2x2). Only channel counts change with different
#' `encoder_channels` values.
#'
#' @param in_channels Number of input channels (default 1).
#' @param stem_stride_rt,stem_stride_cv Stride in the stem conv layer
#'   (default 4 and 1, reflecting the typical 30:1 RT:CV pixel ratio).
#' @param encoder_channels Length-5 integer vector of per-stage channel
#'   counts. Default `c(16L, 16L, 32L, 64L, 64L)` matches the previous
#'   fixed architecture exactly (no projection layer added). Wider
#'   examples: `c(16L, 32L, 64L, 128L, 128L)` for 128-dim latent;
#'   `c(16L, 32L, 64L, 128L, 256L)` for 256-dim latent with a 1x1
#'   projection on top.
#' @param dropout_p Channel-dropout probability applied via
#'   `nn_dropout2d` between residual block groups (after layer1,
#'   layer2, layer3). Whole feature-map channels are zeroed independently
#'   with probability `dropout_p`. Standard regularization for
#'   convolutional encoders; more transferable than element-wise dropout
#'   for spatially-correlated feature maps. Default `0.0` disables
#'   dropout (identity behavior). Recommended try-values `0.1` and `0.2`
#'   when overfit past a certain epoch is observed on real-validation
#'   MSE. Applied only in `training` mode (encoder$eval() disables it
#'   automatically for downstream feature extraction and inference).
#' @export
dms_encoder <- torch::nn_module(
  initialize = function(in_channels = 1L,
                        stem_stride_rt = 4L, stem_stride_cv = 1L,
                        encoder_channels = c(16L, 16L, 32L, 64L, 64L),
                        dropout_p = 0.0) {
    encoder_channels <- as.integer(encoder_channels)
    stopifnot(length(encoder_channels) == 5L,
              all(encoder_channels > 0))
    if (dropout_p < 0 || dropout_p >= 1)
      stop("dropout_p must be in [0, 1).")
    self$encoder_channels <- encoder_channels
    self$dropout_p <- dropout_p
    c1 <- encoder_channels[1]; c2 <- encoder_channels[2]
    c3 <- encoder_channels[3]; c4 <- encoder_channels[4]
    c5 <- encoder_channels[5]
    self$stem <- torch::nn_sequential(
      torch::nn_conv2d(in_channels, c1, kernel_size = 7L,
                stride = c(stem_stride_rt, stem_stride_cv),
                padding = 3L, bias = FALSE),
      torch::nn_batch_norm2d(c1),
      torch::nn_relu(),
      torch::nn_max_pool2d(kernel_size = 3L, stride = c(2L, 1L), padding = 1L)
    )
    self$layer1 <- torch::nn_sequential(
      res_block(c1, c2, stride = 1L),
      res_block(c2, c2, stride = 1L)
    )
    self$layer2 <- torch::nn_sequential(
      res_block(c2, c3, stride = 2L),
      res_block(c3, c3, stride = 1L)
    )
    self$layer3 <- torch::nn_sequential(
      res_block(c3, c4, stride = 2L),
      res_block(c4, c4, stride = 1L)
    )
    # Optional 1x1 projection to latent_ch. Skipped when c5 == c4 so
    # default architecture is exactly preserved (no extra parameters).
    self$has_projection <- c5 != c4
    if (self$has_projection) {
      self$projection <- torch::nn_sequential(
        torch::nn_conv2d(c4, c5, kernel_size = 1L, bias = FALSE),
        torch::nn_batch_norm2d(c5),
        torch::nn_relu()
      )
    }
    # Single stateless dropout module reused across the three
    # inter-layer positions. Instantiated only when dropout is on so
    # the default-arg encoder has zero extra parameters/modules.
    self$use_dropout <- dropout_p > 0
    if (self$use_dropout) {
      self$drop <- torch::nn_dropout2d(p = dropout_p)
    }
  },
  forward = function(x) {
    x <- self$stem(x)
    x <- self$layer1(x)
    if (self$use_dropout) x <- self$drop(x)
    x <- self$layer2(x)
    if (self$use_dropout) x <- self$drop(x)
    x <- self$layer3(x)
    if (self$use_dropout) x <- self$drop(x)
    if (self$has_projection) x <- self$projection(x)
    x
  }
)

#' GC-DMS decoder (mirror of encoder, transposed convolutions)
#'
#' Used during pretraining to reconstruct synthetic inputs through the
#' encoder bottleneck. Crops output to the exact `target_H` x `target_W`
#' if transposed-conv strides produce slightly larger output.
#'
#' Channel schedule mirrors the encoder: input is `latent_ch`
#' (= `encoder_channels[5]`), optional 1x1 unprojection to `layer3_ch`
#' (skipped when equal), then three transposed convs
#' `layer3_ch -> layer2_ch -> layer1_ch -> stem_ch -> 1`.
#'
#' @param target_H,target_W Input image dimensions.
#' @param stem_stride_rt,stem_stride_cv Must match encoder configuration.
#' @param encoder_channels Length-5 integer vector, same as the encoder's.
#'   Default `c(16L, 16L, 32L, 64L, 64L)` matches previous behavior.
#' @export
dms_decoder <- torch::nn_module(
  initialize = function(target_H, target_W,
                        stem_stride_rt = 4L, stem_stride_cv = 1L,
                        encoder_channels = c(16L, 16L, 32L, 64L, 64L)) {
    encoder_channels <- as.integer(encoder_channels)
    stopifnot(length(encoder_channels) == 5L)
    self$target_H <- target_H
    self$target_W <- target_W
    self$encoder_channels <- encoder_channels
    c1 <- encoder_channels[1]; c2 <- encoder_channels[2]
    c3 <- encoder_channels[3]; c4 <- encoder_channels[4]
    c5 <- encoder_channels[5]
    # Optional 1x1 unprojection: c5 -> c4 (skipped when equal).
    self$has_unprojection <- c5 != c4
    if (self$has_unprojection) {
      self$unprojection <- torch::nn_sequential(
        torch::nn_conv2d(c5, c4, kernel_size = 1L, bias = FALSE),
        torch::nn_batch_norm2d(c4),
        torch::nn_relu()
      )
    }
    self$up3 <- torch::nn_sequential(
      torch::nn_conv_transpose2d(c4, c3, kernel_size = 3L, stride = 2L,
                          padding = 1L, output_padding = 1L),
      torch::nn_batch_norm2d(c3), torch::nn_relu()
    )
    self$up2 <- torch::nn_sequential(
      torch::nn_conv_transpose2d(c3, c2, kernel_size = 3L, stride = 2L,
                          padding = 1L, output_padding = 1L),
      torch::nn_batch_norm2d(c2), torch::nn_relu()
    )
    self$up1 <- torch::nn_sequential(
      torch::nn_conv_transpose2d(c2, c1, kernel_size = 3L,
                          stride = c(2L, 1L), padding = 1L,
                          output_padding = c(1L, 0L)),
      torch::nn_batch_norm2d(c1), torch::nn_relu()
    )
    self$up0 <- torch::nn_conv_transpose2d(c1, 1L, kernel_size = 7L,
                                     stride = c(stem_stride_rt, stem_stride_cv),
                                     padding = 3L,
                                     output_padding = c(stem_stride_rt - 1L, 0L))
  },
  forward = function(x) {
    if (self$has_unprojection) x <- self$unprojection(x)
    x <- self$up3(x); x <- self$up2(x); x <- self$up1(x); x <- self$up0(x)
    h <- x$size()[3]; w <- x$size()[4]
    if (h != self$target_H || w != self$target_W) {
      x <- x[, , 1:self$target_H, 1:self$target_W]
    }
    x
  }
)

#' Composed denoising autoencoder (encoder + decoder)
#'
#' Forward pass: `x -> encoder -> z -> decoder -> reconstruction`. Trained
#' on synthetic images (noisy input, clean target) via MSE loss in the
#' pretraining functions; the resulting encoder transfers to downstream
#' classification.
#'
#' @inheritParams dms_decoder
#' @param in_channels Input channels (default 1).
#' @param dropout_p Encoder-side channel dropout probability
#'   (see [dms_encoder()]). Decoder stays clean — dropout is only
#'   applied to the encoder because that's the part transferred to
#'   downstream classification. Default `0.0`.
#' @export
dms_denoising_autoencoder <- torch::nn_module(
  initialize = function(target_H, target_W, in_channels = 1L,
                        stem_stride_rt = 4L, stem_stride_cv = 1L,
                        encoder_channels = c(16L, 16L, 32L, 64L, 64L),
                        dropout_p = 0.0) {
    self$encoder_channels <- as.integer(encoder_channels)
    self$encoder <- dms_encoder(in_channels = in_channels,
                                 stem_stride_rt = stem_stride_rt,
                                 stem_stride_cv = stem_stride_cv,
                                 encoder_channels = encoder_channels,
                                 dropout_p = dropout_p)
    self$decoder <- dms_decoder(target_H = target_H, target_W = target_W,
                                 stem_stride_rt = stem_stride_rt,
                                 stem_stride_cv = stem_stride_cv,
                                 encoder_channels = encoder_channels)
  },
  forward = function(x) {
    z <- self$encoder(x)
    self$decoder(z)
  }
)
