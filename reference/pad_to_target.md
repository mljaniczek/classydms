# Center-pad a matrix to target dimensions

Pads `Z` with `pad_value` so the result has dimensions
`target_H x target_W`. The original matrix is placed in the center of
the padded canvas. Errors if `Z` is already larger than the target along
either axis.

## Usage

``` r
pad_to_target(Z, target_H, target_W, pad_value = 0)
```

## Arguments

- Z:

  A 2-D numeric matrix.

- target_H, target_W:

  Target row and column counts.

- pad_value:

  Fill value for the padded margins (default 0).

## Value

A `target_H x target_W` matrix.
