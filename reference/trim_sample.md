# Apply trim bounds to a sample (clamps to actual size if smaller)

Apply trim bounds to a sample (clamps to actual size if smaller)

## Usage

``` r
trim_sample(sample, bounds)
```

## Arguments

- sample:

  A sample list (output of
  [`read_dms_file()`](https://mljaniczek.github.io/classydms/reference/read_dms_file.md)
  or
  [`process_one_sample()`](https://mljaniczek.github.io/classydms/reference/process_one_sample.md)),
  with elements `Z`, `time`, `cv`, and optionally `Z_pretrim`.

- bounds:

  A bounding box from
  [`trim_bounds_from_occupancy()`](https://mljaniczek.github.io/classydms/reference/trim_bounds_from_occupancy.md)
  or
  [`pooled_trim_bounds()`](https://mljaniczek.github.io/classydms/reference/pooled_trim_bounds.md),
  with `rt_start`, `rt_end`, `cv_start`, `cv_end`.

## Value

The sample with its `Z` (and `Z_pretrim`, if present) trimmed to
`bounds`, along with the corresponding subset of the `time` and `cv`
axes. Indices are clamped to the actual matrix dimensions so a sample
smaller than the trim bounds is not an error.
