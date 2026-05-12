# GC-DMS decoder (mirror of encoder, transposed convolutions)

Used during pretraining to reconstruct synthetic inputs through the
encoder bottleneck. Crops output to the exact `target_H` x `target_W` if
transposed-conv strides produce slightly larger output.

## Usage

``` r
dms_decoder(target_H, target_W, stem_stride_rt = 4L, stem_stride_cv = 1L)
```

## Arguments

- target_H, target_W:

  Input image dimensions.

- stem_stride_rt, stem_stride_cv:

  Must match encoder configuration.
