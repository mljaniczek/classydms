# Catalog compound positions on a padded canvas

Maps every compound's observations onto an `H x W` canvas. When the
catalog records a
[`dms_grid()`](https://mljaniczek.github.io/classydms/reference/dms_grid.md),
observations are placed from their seconds / volts through the grid
(centred if the canvas is larger). Otherwise they are laid out the way
[`pad_to_target()`](https://mljaniczek.github.io/classydms/reference/pad_to_target.md)
lays out real samples: each observation's fraction position is converted
with its own sample's unpadded data dims
(`catalog$parameters$data_dims`) and centre-padding offset. This is the
frame the synthetic generators render in and the frame the encoder sees
for real samples.

## Usage

``` r
catalog_pixel_positions(catalog, H, W)
```

## Arguments

- catalog:

  A `peak_catalog`.

- H, W:

  Canvas dimensions in pixels.

## Value

A tibble with one row per compound: `compound_id`, `rt_px`, `cv_px`
(median 1-based canvas pixel over the compound's observations, or the
centroid when it has none), and `rt_px_sd`, `cv_px_sd` (SD of the
observations' canvas pixels; 0 with fewer than two).

## Details

Catalogs built before data dims were recorded fall back to
`fraction * H` (with a once-per-session warning).
