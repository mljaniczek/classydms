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
  stem_stride_cv = 1L
)
```

## Arguments

- target_H, target_W:

  Input image dimensions.

- in_channels:

  Input channels (default 1).

- stem_stride_rt, stem_stride_cv:

  Must match encoder configuration.
