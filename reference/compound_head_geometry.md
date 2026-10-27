# Catalog geometry for a compound-anchored head

Collects each compound's position on the `H x W` canvas, expressed as a
fraction in the same convention as the head's token centres (pixel `p`
sits at `(p - 0.5) / H`), and a positional spread used as the width of
the head's location prior. Positions come from
[`catalog_pixel_positions()`](https://mljaniczek.github.io/classydms/reference/catalog_pixel_positions.md),
which maps every observation onto the canvas the way
[`pad_to_target()`](https://mljaniczek.github.io/classydms/reference/pad_to_target.md)
centre-pads real samples, so the prior is centred where the peak sits in
both synthetic and real images. The spread is the SD of the
observations' canvas positions combined with any extra displacement the
pretraining augmentation can introduce (affine shift, RT warp), so the
prior stays wide enough to find a shifted peak. (The head additionally
floors the width at one feature-map cell.)

## Usage

``` r
compound_head_geometry(
  catalog,
  H,
  W,
  compound_ids = NULL,
  extra_rt_px = 0,
  extra_cv_px = 0
)
```

## Arguments

- catalog:

  A `peak_catalog`.

- H, W:

  Image dimensions in pixels.

- compound_ids:

  Optional subset of `compound_id`s. Default: every compound with a
  finite position.

- extra_rt_px, extra_cv_px:

  Additional displacement (pixels) to fold into the prior width, e.g.
  the maximum affine shift.

## Value

List with `compound_id`, `rt_frac`, `cv_frac`, `sd_rt_frac`,
`sd_cv_frac`, all as fractions of `H` / `W`.
