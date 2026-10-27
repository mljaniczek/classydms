# Estimate peak parameter distributions from real GC-DMS samples

Pools all samples (cases and controls together; no labels used) and
extracts empirical distributions of: number of peaks per sample, peak
locations (as fractions of image dimensions), peak widths (sigma_RT,
sigma_CV) measured at `frac_height` height and converted to Gaussian
sigmas, peak intensities, and background noise mean / SD. Widths are
fitted to log-normal distributions.

## Usage

``` r
estimate_peak_params(
  Z_list,
  Z_pretrim_list = NULL,
  time_axes = NULL,
  cv_axes = NULL,
  top_k = 150,
  frac_height = 0.25,
  min_sep_rt = 8,
  min_sep_cv = 2,
  grid = NULL,
  verbose = TRUE
)
```

## Arguments

- Z_list:

  List of trimmed, normalized intensity matrices (real data).

- Z_pretrim_list:

  Optional list of pre-dust-threshold matrices.

- top_k:

  Maximum peaks to detect per sample (default 150).

- frac_height:

  Fraction of peak height defining the width (default 0.25; 0.5
  corresponds to FWHM).

- min_sep_rt, min_sep_cv:

  Minimum pixel separation between peaks.

- grid:

  Optional
  [`dms_grid()`](https://mljaniczek.github.io/classydms/reference/dms_grid.md)
  every matrix in `Z_list` was resampled onto
  ([`resample_to_grid()`](https://mljaniczek.github.io/classydms/reference/resample_to_grid.md)).
  Supplies `time_axes` / `cv_axes` for every sample, requires each
  matrix to be exactly `grid$H x grid$W`, and is recorded on the result
  so catalogs and generators place peaks from seconds / volts through
  it.

- verbose:

  Whether to print diagnostic messages (default TRUE).

## Value

A list of fitted distributions and raw observation vectors, suitable for
passing to
[`generate_one_synthetic()`](https://mljaniczek.github.io/classydms/reference/generate_one_synthetic.md)
or
[`generate_synthetic_dataset()`](https://mljaniczek.github.io/classydms/reference/generate_synthetic_dataset.md).
Fields `rt_loc_raw`, `cv_loc_raw`, `sigma_rt_raw`, `sigma_cv_raw`, and
`intensity_raw` are parallel per-peak vectors; `sample_idx_raw` gives
the 1-based index of the source sample for each peak (parallel to the
other `*_raw` vectors), and `sample_names` is `names(Z_list)` for
lookup. `data_dims` is an `n_samples x 2` matrix (`rows`, `cols`) of
each sample's matrix dims, which the generators use to place fraction
positions on a centre-padded canvas. These enable per-sample analyses
(e.g. sample-level shift/clustering) without changing existing callers.

## Details

If `Z_pretrim_list` is provided (pre-dust-threshold matrices), peak
width measurement uses the pre-thresholded data as a fallback to
preserve peak flanks that the dust thresholding zeroed out.
