# GC-DMS decoder (mirror of encoder, transposed convolutions)

Used during pretraining to reconstruct synthetic inputs through the
encoder bottleneck. Crops output to the exact `target_H` x `target_W` if
transposed-conv strides produce slightly larger output.

## Usage

``` r
dms_decoder(
  target_H,
  target_W,
  stem_stride_rt = 4L,
  stem_stride_cv = 1L,
  encoder_channels = c(16L, 16L, 32L, 64L, 64L)
)
```

## Arguments

- target_H, target_W:

  Input image dimensions.

- stem_stride_rt, stem_stride_cv:

  Must match encoder configuration.

- encoder_channels:

  Length-5 integer vector, same as the encoder's. Default
  `c(16L, 16L, 32L, 64L, 64L)` matches previous behavior.

## Details

Channel schedule mirrors the encoder: input is `latent_ch` (=
`encoder_channels[5]`), optional 1x1 unprojection to `layer3_ch`
(skipped when equal), then three transposed convs
`layer3_ch -> layer2_ch -> layer1_ch -> stem_ch -> 1`.
