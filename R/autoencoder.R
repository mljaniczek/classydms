# classydms: denoising autoencoder architecture
# Encoder + decoder + composed denoising autoencoder. The encoder is the
# part that is transferred to downstream classification; the decoder is
# discarded after pretraining (except for reconstruction QC).

#' GC-DMS encoder (ResNet-style, asymmetric stride)
#'
#' Three residual block groups (16 -> 32 -> 64 channels) with an asymmetric
#' stem stride that aggressively downsamples RT while preserving CV
#' resolution. Output is a `(B, 64, h, w)` feature map suitable as input to
#' a classification head or for spatial saliency mapping.
#'
#' @param in_channels Number of input channels (default 1).
#' @param stem_stride_rt,stem_stride_cv Stride in the stem conv layer
#'   (default 4 and 1, reflecting the typical 30:1 RT:CV pixel ratio).
#' @export
dms_encoder <- torch::nn_module(
  initialize = function(in_channels = 1L,
                        stem_stride_rt = 4L, stem_stride_cv = 1L) {
    self$stem <- torch::nn_sequential(
      torch::nn_conv2d(in_channels, 16L, kernel_size = 7L,
                stride = c(stem_stride_rt, stem_stride_cv),
                padding = 3L, bias = FALSE),
      torch::nn_batch_norm2d(16L),
      torch::nn_relu(),
      torch::nn_max_pool2d(kernel_size = 3L, stride = c(2L, 1L), padding = 1L)
    )
    self$layer1 <- torch::nn_sequential(
      res_block(16L, 16L, stride = 1L),
      res_block(16L, 16L, stride = 1L)
    )
    self$layer2 <- torch::nn_sequential(
      res_block(16L, 32L, stride = 2L),
      res_block(32L, 32L, stride = 1L)
    )
    self$layer3 <- torch::nn_sequential(
      res_block(32L, 64L, stride = 2L),
      res_block(64L, 64L, stride = 1L)
    )
  },
  forward = function(x) {
    x <- self$stem(x)
    x <- self$layer1(x)
    x <- self$layer2(x)
    x <- self$layer3(x)
    x
  }
)

#' GC-DMS decoder (mirror of encoder, transposed convolutions)
#'
#' Used during pretraining to reconstruct synthetic inputs through the
#' encoder bottleneck. Crops output to the exact `target_H` x `target_W`
#' if transposed-conv strides produce slightly larger output.
#'
#' @param target_H,target_W Input image dimensions.
#' @param stem_stride_rt,stem_stride_cv Must match encoder configuration.
#' @export
dms_decoder <- torch::nn_module(
  initialize = function(target_H, target_W,
                        stem_stride_rt = 4L, stem_stride_cv = 1L) {
    self$target_H <- target_H
    self$target_W <- target_W
    self$up3 <- torch::nn_sequential(
      torch::nn_conv_transpose2d(64L, 32L, kernel_size = 3L, stride = 2L,
                          padding = 1L, output_padding = 1L),
      torch::nn_batch_norm2d(32L), torch::nn_relu()
    )
    self$up2 <- torch::nn_sequential(
      torch::nn_conv_transpose2d(32L, 16L, kernel_size = 3L, stride = 2L,
                          padding = 1L, output_padding = 1L),
      torch::nn_batch_norm2d(16L), torch::nn_relu()
    )
    self$up1 <- torch::nn_sequential(
      torch::nn_conv_transpose2d(16L, 16L, kernel_size = 3L,
                          stride = c(2L, 1L), padding = 1L,
                          output_padding = c(1L, 0L)),
      torch::nn_batch_norm2d(16L), torch::nn_relu()
    )
    self$up0 <- torch::nn_conv_transpose2d(16L, 1L, kernel_size = 7L,
                                     stride = c(stem_stride_rt, stem_stride_cv),
                                     padding = 3L,
                                     output_padding = c(stem_stride_rt - 1L, 0L))
  },
  forward = function(x) {
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
#' @export
dms_denoising_autoencoder <- torch::nn_module(
  initialize = function(target_H, target_W, in_channels = 1L,
                        stem_stride_rt = 4L, stem_stride_cv = 1L) {
    self$encoder <- dms_encoder(in_channels = in_channels,
                                 stem_stride_rt = stem_stride_rt,
                                 stem_stride_cv = stem_stride_cv)
    self$decoder <- dms_decoder(target_H = target_H, target_W = target_W,
                                 stem_stride_rt = stem_stride_rt,
                                 stem_stride_cv = stem_stride_cv)
  },
  forward = function(x) {
    z <- self$encoder(x)
    self$decoder(z)
  }
)
