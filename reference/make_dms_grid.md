# Build a dms_grid from start, step and size

Build a dms_grid from start, step and size

## Usage

``` r
make_dms_grid(H, W, t0 = 0, dt = 1, cv0, dcv = 0.55)
```

## Arguments

- H, W:

  Grid rows / columns. Choose them as the model input size (H a multiple
  of 32, W a multiple of 4 for
  [`dms_encoder()`](https://mljaniczek.github.io/classydms/reference/dms_encoder.md));
  rows past a sample's last RT are filled with zeros by
  [`resample_to_grid()`](https://mljaniczek.github.io/classydms/reference/resample_to_grid.md).

- t0, dt:

  First RT (seconds) and RT step.

- cv0, dcv:

  First CV (volts) and CV step. Use
  [`choose_cv_window()`](https://mljaniczek.github.io/classydms/reference/choose_cv_window.md)
  to pick `cv0` from pooled signal occupancy.

## Value

A `dms_grid`.
