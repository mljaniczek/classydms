# Greedily select top-K peak centers with spatial separation

Returns up to `top_k` local-maximum pixels, ranked by intensity, subject
to a minimum separation constraint to prevent double-counting
overlapping peaks.

## Usage

``` r
pick_peak_centers(Z, top_k = 30, eps = NULL, min_sep_rt = 8, min_sep_cv = 2)
```

## Arguments

- Z:

  Intensity matrix.

- top_k:

  Maximum number of peaks to return.

- eps:

  Threshold below which pixels are not considered peak candidates. If
  `NULL`, uses
  [`robust_eps()`](https://mljaniczek.github.io/classydms/reference/robust_eps.md).

- min_sep_rt, min_sep_cv:

  Minimum pixel separation between peaks.

## Value

Integer matrix with two columns (`r`, `c`) giving peak indices.
