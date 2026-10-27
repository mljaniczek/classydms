# GC-DMS encoder (ResNet-style, asymmetric stride)

Three residual block groups with an asymmetric stem stride that
aggressively downsamples RT while preserving CV resolution. Output is a
`(B, encoder_channels[5], h, w)` feature map suitable as input to a
classification head or for spatial saliency mapping.

## Usage

``` r
dms_encoder(
  in_channels = 1L,
  stem_stride_rt = 4L,
  stem_stride_cv = 1L,
  encoder_channels = c(16L, 16L, 32L, 64L, 64L),
  dropout_p = 0
)
```

## Arguments

- in_channels:

  Number of input channels (default 1).

- stem_stride_rt, stem_stride_cv:

  Stride in the stem conv layer (default 4 and 1, reflecting the typical
  30:1 RT:CV pixel ratio).

- encoder_channels:

  Length-5 integer vector of per-stage channel counts. Default
  `c(16L, 16L, 32L, 64L, 64L)` matches the previous fixed architecture
  exactly (no projection layer added). Wider examples:
  `c(16L, 32L, 64L, 128L, 128L)` for 128-dim latent;
  `c(16L, 32L, 64L, 128L, 256L)` for 256-dim latent with a 1x1
  projection on top.

- dropout_p:

  Channel-dropout probability applied via `nn_dropout2d` between
  residual block groups (after layer1, layer2, layer3). Whole
  feature-map channels are zeroed independently with probability
  `dropout_p`. Standard regularization for convolutional encoders; more
  transferable than element-wise dropout for spatially-correlated
  feature maps. Default `0.0` disables dropout (identity behavior).
  Recommended try-values `0.1` and `0.2` when overfit past a certain
  epoch is observed on real-validation MSE. Applied only in `training`
  mode (encoder\$eval() disables it automatically for downstream feature
  extraction and inference).

## Details

Channel schedule via `encoder_channels`, a length-5 integer vector
`c(stem_ch, layer1_ch, layer2_ch, layer3_ch, latent_ch)`:

- Stem: `1 -> stem_ch` (kernel 7, stride
  `(stem_stride_rt, stem_stride_cv)`)

- MaxPool: stride `(2, 1)`

- Layer1: `stem_ch -> layer1_ch`, 2 residual blocks at stride 1

- Layer2: `layer1_ch -> layer2_ch`, 2 blocks with first stride 2

- Layer3: `layer2_ch -> layer3_ch`, 2 blocks with first stride 2

- Optional 1x1 projection: `layer3_ch -> latent_ch` (skipped when equal)

Total spatial downsample is fixed at 32x on RT and 4x on CV regardless
of channel counts (stem 4x1, maxpool 2x1, layer2 2x2, layer3 2x2). Only
channel counts change with different `encoder_channels` values.
