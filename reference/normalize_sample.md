# Log-quantile normalize a single sample

Applies `log(1 + x) / Q_q(log(1 + x))`, where `Q_q` is the q-th quantile
of the log-transformed image. Falls back to dividing by 1 if the
quantile is below `eps` or non-finite. Per-sample normalization makes
intensity comparisons robust to instrument/subject-level absolute level
differences.

## Usage

``` r
normalize_sample(Z, q = 0.95, eps = 1e-08)
```

## Arguments

- Z:

  A non-negative numeric matrix.

- q:

  Quantile used as the normalization denominator (default 0.95).

- eps:

  Numerical floor on the denominator (default 1e-8).

## Value

A normalized matrix of the same shape as `Z`.
