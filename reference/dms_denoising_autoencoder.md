# Composed denoising autoencoder (encoder + decoder)

Forward pass: `x -> encoder -> z -> decoder -> reconstruction`. Trained
on synthetic images (noisy input, clean target) via MSE loss in the
pretraining functions; the resulting encoder transfers to downstream
classification.

## Usage

``` r
dms_denoising_autoencoder(
  target_H,
  target_W,
  in_channels = 1L,
  stem_stride_rt = 4L,
  stem_stride_cv = 1L,
  encoder_channels = c(16L, 16L, 32L, 64L, 64L),
  dropout_p = 0
)
```

## Arguments

- target_H, target_W:

  Input image dimensions.

- in_channels:

  Input channels (default 1).

- stem_stride_rt, stem_stride_cv:

  Must match encoder configuration.

- encoder_channels:

  Length-5 integer vector, same as the encoder's. Default
  `c(16L, 16L, 32L, 64L, 64L)` matches previous behavior.

- dropout_p:

  Encoder-side channel dropout probability (see
  [`dms_encoder()`](https://mljaniczek.github.io/classydms/reference/dms_encoder.md)).
  Decoder stays clean — dropout is only applied to the encoder because
  that's the part transferred to downstream classification. Default
  `0.0`.
