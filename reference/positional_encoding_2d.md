# 2D sinusoidal positional encoding

Returns a `(d, h, w)` tensor: the first half of the channels encode the
row (RT) position, the second half the column (CV) position, using the
standard transformer sin/cos frequency ladder.

## Usage

``` r
positional_encoding_2d(h, w, d, device = "cpu")
```
