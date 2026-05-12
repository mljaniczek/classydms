# Compute padding targets across a list of matrices

Returns target H and W (rounded up to the nearest multiple of
`multiple`) so that all matrices in `Z_list` can be padded to a common
size. The `multiple = 32` default ensures clean compatibility with the
strided convolutions in the encoder (total stride along RT is 32).

## Usage

``` r
compute_pad_targets(Z_list, multiple = 32L)
```

## Arguments

- Z_list:

  List of 2-D numeric matrices.

- multiple:

  Integer multiple to round target dims up to (default 32).

## Value

List with `H` and `W` (the target dims), plus `raw_max_H`, `raw_max_W`
(un-rounded maxes) and per-axis range vectors for diagnostics.
