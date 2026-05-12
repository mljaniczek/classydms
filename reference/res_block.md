# Residual block used by encoder and tiny ResNet variants

Two 3x3 convs with batch norm and a skip connection; optionally
downsamples via a 1x1 conv when `stride != 1` or channel counts differ.

## Usage

``` r
res_block(in_ch, out_ch, stride = 1)
```
