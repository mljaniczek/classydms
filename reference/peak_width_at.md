# Measure peak width at a fraction of peak height (FWHM-style)

Scans outward from a peak center along rows and columns until the
intensity drops below `frac_height` times the center value, returning
the width in pixels along each axis. Bounded by `max_radius_*` to
prevent walking into adjacent peaks.

## Usage

``` r
peak_width_at(
  Z,
  r0,
  c0,
  frac_height = 0.5,
  eps = 0,
  max_radius_rt = Inf,
  max_radius_cv = Inf
)
```

## Arguments

- Z:

  Intensity matrix.

- r0, c0:

  Peak center indices (1-based).

- frac_height:

  Fraction of peak height defining the width boundary (default 0.5,
  i.e., FWHM).

- eps:

  Floor below which intensities are treated as zero.

- max_radius_rt, max_radius_cv:

  Maximum walk distance per axis.

## Value

Named numeric vector with `rt` and `cv` widths in pixels.
