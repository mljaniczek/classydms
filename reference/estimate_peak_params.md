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
  top_k = 150,
  frac_height = 0.25,
  min_sep_rt = 8,
  min_sep_cv = 2,
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

- verbose:

  Whether to print diagnostic messages (default TRUE).

## Value

A list of fitted distributions and raw observation vectors, suitable for
passing to
[`generate_one_synthetic()`](https://mljaniczek.github.io/classydms/reference/generate_one_synthetic.md)
or
[`generate_synthetic_dataset()`](https://mljaniczek.github.io/classydms/reference/generate_synthetic_dataset.md).

## Details

If `Z_pretrim_list` is provided (pre-dust-threshold matrices), peak
width measurement uses the pre-thresholded data as a fallback to
preserve peak flanks that the dust thresholding zeroed out.
