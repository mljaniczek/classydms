# GC-DMS encoder (ResNet-style, asymmetric stride)

Three residual block groups (16 -\> 32 -\> 64 channels) with an
asymmetric stem stride that aggressively downsamples RT while preserving
CV resolution. Output is a `(B, 64, h, w)` feature map suitable as input
to a classification head or for spatial saliency mapping.

## Usage

``` r
dms_encoder(in_channels = 1L, stem_stride_rt = 4L, stem_stride_cv = 1L)
```

## Arguments

- in_channels:

  Number of input channels (default 1).

- stem_stride_rt, stem_stride_cv:

  Stride in the stem conv layer (default 4 and 1, reflecting the typical
  30:1 RT:CV pixel ratio).
