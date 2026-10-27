# Flag samples that fail QC checks

Applies threshold-based and IQR-based criteria to a QC table from
[`qc_features()`](https://mljaniczek.github.io/classydms/reference/qc_features.md)
and returns the table with boolean flag columns added. A sample is
marked `flagged = TRUE` if it triggers any of the active criteria.
Thresholds are tunable; pass `NA` or `Inf` to disable a specific check.

## Usage

``` r
flag_junk_samples(
  qc_tbl,
  min_prop_nonzero = 0.1,
  max_saturation = Inf,
  min_dynamic_range = 5,
  iqr_mult = 3,
  min_cohort_correlation = 0.3,
  features_for_iqr = c("total_signal", "mean_positive", "n_peaks", "prop_nonzero",
    "max_intensity", "dynamic_range", "bbox_size")
)
```

## Arguments

- qc_tbl:

  Output of
  [`qc_features()`](https://mljaniczek.github.io/classydms/reference/qc_features.md).

- min_prop_nonzero:

  Threshold below which a sample is flagged as low-signal (default
  0.10).

- max_saturation:

  Threshold at or above which a sample is flagged for saturation
  (default `Inf` = disabled; set to e.g. 65535 for 16-bit detectors, or
  whatever the instrument's max readout is).

- min_dynamic_range:

  Threshold below which dynamic range is considered too low (default 5).

- iqr_mult:

  Multiplier for IQR-based outlier detection across univariate features
  (default 3).

- min_cohort_correlation:

  Threshold below which a sample's correlation with the cohort median
  spectrum triggers a flag (default 0.3; only applied if
  `cohort_correlation` is in the table).

- features_for_iqr:

  Character vector of QC feature names to include in IQR outlier
  detection. Default uses the package-recommended set; pass NULL to
  disable IQR checks entirely.

## Value

The input tibble with additional columns: one `*_outlier` column per
IQR-checked feature, plus `low_signal`, `saturation`,
`low_dynamic_range`, `low_cohort_correlation`, and a final `flagged`
column that is the OR of all criteria.

## Details

Default criteria (Tier 1 + Tier 2 from the package design):

- **Low signal**: `prop_nonzero` below `min_prop_nonzero`.

- **Saturation**: `max_intensity` at or above `max_saturation`.

- **Low dynamic range**: `dynamic_range` below `min_dynamic_range`.

- **Univariate IQR outliers** on `total_signal`, `mean_positive`,
  `n_peaks`, `prop_nonzero`, `max_intensity`, `dynamic_range`, and
  `bbox_size` (any feature more than `iqr_mult` IQR from the cohort
  median triggers).

- **Low cohort correlation** (only if `cohort_correlation` is in the
  table): correlation below `min_cohort_correlation`.
