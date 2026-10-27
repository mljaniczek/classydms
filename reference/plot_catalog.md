# Plot the spatial distribution of a peak catalog

Compound positions are shown as red points sized by prevalence, overlaid
on a hexbin density map of the underlying pool peaks. Small clusters and
singletons are optionally overlaid as blue points for context.

## Usage

``` r
plot_catalog(catalog, peak_params, show_small = TRUE, show_singletons = FALSE)
```

## Arguments

- catalog:

  A `peak_catalog` object.

- peak_params:

  The `peak_params` used to build the catalog.

- show_small, show_singletons:

  Overlay the below-threshold clusters and singletons? Default `TRUE`
  for `small`, `FALSE` for `singletons` (there are usually many).

## Value

A ggplot object.

## Details

Requires `peak_params` (the object used to build the catalog) so we can
draw the underlying peak density.
