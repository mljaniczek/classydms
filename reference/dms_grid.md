# Define a common physical grid for GC-DMS images

Define a common physical grid for GC-DMS images

## Usage

``` r
dms_grid(time, cv, tol = 1e-06)
```

## Arguments

- time:

  Regular numeric vector of RT positions (seconds), one per image row.

- cv:

  Regular numeric vector of CV positions (volts), one per image column,
  increasing.

- tol:

  Relative tolerance for the regular-spacing check.

## Value

An object of class `dms_grid`: a list with `time`, `cv`, `H`, `W`, `dt`,
`dcv`.

## See also

[`make_dms_grid()`](https://mljaniczek.github.io/classydms/reference/make_dms_grid.md),
[`resample_to_grid()`](https://mljaniczek.github.io/classydms/reference/resample_to_grid.md)
