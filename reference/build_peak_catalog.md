# Build a peak catalog from a peak_params object

Greedy clustering: peaks are sorted by intensity descending, each
unassigned peak becomes a new cluster seed, and all remaining unassigned
peaks within `(eps_rt, eps_cv)` of the seed are added to that cluster.
Clusters are then classified into three groups:

- compounds:

  clusters with `>= min_cluster_size` observations AND spanning
  `>= min_support_frac` of samples. These are the "reliable positions"
  the catalog is built around.

- small_clusters:

  clusters with \>= 2 observations but below the sample-support
  threshold. Might be sporadic contamination, small compounds, or peaks
  from just a few subjects.

- singletons:

  clusters of size 1 — peaks that didn't co-locate with any other pool
  peak. Often noise, but worth inspecting.

## Usage

``` r
build_peak_catalog(
  peak_params,
  eps_rt = 0.01,
  eps_cv = 0.05,
  min_support_frac = 0.1,
  min_cluster_size = 2L,
  coord_mode = c("fraction", "physical"),
  align = FALSE
)
```

## Arguments

- peak_params:

  Output of
  [`estimate_peak_params()`](https://mljaniczek.github.io/classydms/reference/estimate_peak_params.md)
  containing `rt_loc_raw`, `cv_loc_raw`, `sigma_rt_raw`, `sigma_cv_raw`,
  `intensity_raw`, `sample_idx_raw`, and `n_samples`.

- eps_rt, eps_cv:

  Clustering radius in fraction-of-image units. Defaults `(0.01, 0.05)`
  roughly match the peak-detector's `min_sep_rt = 8, min_sep_cv = 2` on
  a 1400 x 40 native image. Tighter → more, smaller clusters. Looser →
  fewer, larger clusters. In `coord_mode = "physical"` they are in
  seconds and volts.

- min_support_frac:

  Fraction of samples a cluster must span to count as a compound.
  Default `0.10`: appear in at least 10% of samples.

- min_cluster_size:

  Minimum number of pool peaks in a cluster for compound status. Default
  `2L` (excludes singletons from compounds regardless of support).

- coord_mode:

  `"fraction"` (default) clusters on `rt_loc_raw` / `cv_loc_raw`;
  `"physical"` clusters on `rt_seconds_raw` / `cv_volts_raw` (requires
  [`estimate_peak_params()`](https://mljaniczek.github.io/classydms/reference/estimate_peak_params.md)
  to have been given `time_axes` / `cv_axes`). Either way the catalog
  stores fraction positions (`rt_frac`, `cv_frac`, `*_frac_obs`) and,
  when available, physical ones (`rt_seconds`, `cv_volts`,
  `rt_seconds_obs`, `cv_volts_obs`), all from the same member peaks.

- align:

  Optional per-sample drift correction applied to `peak_params` before
  clustering. `FALSE` (default) skips alignment. `TRUE` calls
  [`align_via_catalog()`](https://mljaniczek.github.io/classydms/reference/align_via_catalog.md)
  with sensible defaults (`eps_rt_loose = 3 * eps_rt`,
  `eps_cv_loose = 3 * eps_cv`, `min_support_frac = min_support_frac`,
  `iterate = TRUE`, `max_iter = 3L`). Pass a named list to override any
  of those. Alignment is idempotent and safe to leave off if drift is
  small; see the compactness-delta table in `compare_all_catalogs.Rmd`
  for how to decide whether it helps on your cohort.

## Value

An object of class `peak_catalog` with fields:

- compounds:

  tibble of compound entries — one row per compound with columns
  `compound_id, rt_loc, cv_loc, location_sd_rt, location_sd_cv, n_observations, n_samples, prevalence`,
  and list-columns
  `sigma_rt_obs, sigma_cv_obs, intensity_obs, rt_frac_obs, cv_frac_obs, member_indices`,
  plus `rt_seconds_obs, cv_volts_obs` when physical positions exist and
  `sample_idx_obs` (source sample of each observation) when data dims
  are known.

- small_clusters:

  same schema, for clusters below thresholds.

- singletons:

  same schema, for clusters of size 1.

- parameters:

  the settings used to build the catalog, plus `data_dims`: the
  `n_samples x 2` matrix of each sample's unpadded `rows, cols` (from
  `peak_params$data_dims`, or the axis lengths in `time_ranges` /
  `cv_ranges`), which generators use to place observations on a
  centre-padded canvas
  ([`catalog_pixel_positions()`](https://mljaniczek.github.io/classydms/reference/catalog_pixel_positions.md));
  and `grid`, the
  [`dms_grid()`](https://mljaniczek.github.io/classydms/reference/dms_grid.md)
  the samples were resampled onto (if recorded), through which
  generators map `rt_seconds_obs` / `cv_volts_obs` to pixels.

## Details

Uses a sorted-index optimization so runtime stays `O(n log n + n·k)`
rather than `O(n^2)`, where `k` is the average number of candidates per
seed. Practical for pool sizes up to millions.
