# Generate one synthetic sample from a peak catalog

For each compound in the catalog, a `Bernoulli(prevalence)` draw decides
whether the compound "fires" in this sample. If it fires, a 2-D Gaussian
peak is placed at `compound$rt_loc, compound$cv_loc` with morphology
bootstrapped from the compound's observed
`(sigma_rt, sigma_cv, intensity)` triples. Per-peak location and size
jitter apply on top. Folded-Normal noise is added at the end when
`add_noise = TRUE`.

## Usage

``` r
generate_one_synthetic_from_catalog(
  catalog,
  H,
  W,
  add_noise = TRUE,
  size_jitter = 0.15,
  location_jitter_rt = 2,
  location_jitter_cv = 1,
  noise_scale = 1,
  prevalence_override = NULL,
  intensity_mult = NULL,
  compound_ids = NULL,
  peak_shape = c("gaussian", "emg", "mix"),
  emg_tau_range = c(0.2, 0.8),
  emg_mix_prob = 0.3
)
```

## Arguments

- catalog:

  A `peak_catalog` object from
  [`build_peak_catalog()`](https://mljaniczek.github.io/classydms/reference/build_peak_catalog.md).

- H, W:

  Output image dimensions in pixels.

- add_noise:

  Whether to add folded-Normal background noise.

- size_jitter:

  Per-peak log-Normal size scale factor SD.

- location_jitter_rt, location_jitter_cv:

  Per-peak location jitter SD (pixels).

- noise_scale:

  Multiplier on `catalog$noise$sd`. Default `1.0`.

- prevalence_override:

  Optional named numeric vector mapping `compound_id` (as character) to
  a replacement prevalence. Not present compounds use the catalog's own
  prevalence.

- intensity_mult:

  Optional named numeric vector mapping `compound_id` (as character) to
  an amplitude multiplier applied when the compound fires. Defaults to
  `1.0` for unlisted compounds.

- compound_ids:

  Optional integer vector of compound_ids to include; `NULL` (default)
  uses every catalog compound.

## Value

List with `clean` (H x W), `noisy` (H x W, `= clean` when
`add_noise = FALSE`), and `fired` (logical vector of length
`nrow(catalog$compounds)` indicating which compounds actually appeared
in this sample).

## Details

Compared to
[`generate_one_synthetic()`](https://mljaniczek.github.io/classydms/reference/generate_one_synthetic.md),
the catalog version reflects the real-data property that the same
compound positions appear across many samples. Which specific compounds
appear per sample is random (via Bernoulli), but the compound
*identities* and positions are shared across all synthetic samples.
