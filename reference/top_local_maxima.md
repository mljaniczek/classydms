# Non-max-suppressed local maxima of a matrix, ranked by value

Very light wrapper on the same greedy NMS used for peak detection — pick
the brightest pixel, then the next-brightest that is more than
`min_sep_*` away from every accepted pixel on at least one axis, up to
`top_k` picks.

## Usage

``` r
top_local_maxima(M, top_k = 10L, min_sep_rt = 8L, min_sep_cv = 2L)
```

## Arguments

- M:

  numeric matrix (typically a saliency map).

- top_k:

  number of local maxima to return.

- min_sep_rt, min_sep_cv:

  separation on each axis (in pixels).

## Value

integer matrix with columns `r`, `c` and `val`, ranked by `val`.
