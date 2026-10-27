# Evaluate biomarker recovery from a saliency map

Two complementary metrics on the same saliency map and ground-truth
biomarker spec:

## Usage

``` r
evaluate_biomarker_recovery(
  saliency_map,
  biomarkers,
  H = NULL,
  W = NULL,
  match_radius_rt = 8L,
  match_radius_cv = 2L,
  top_k = 10L
)
```

## Arguments

- saliency_map:

  numeric `H` x `W` matrix in the same coordinate frame as the biomarker
  locations (i.e. after any padding / trimming used for the encoder
  input).

- biomarkers:

  ground-truth spec: tibble with at least `rt_loc` and `cv_loc` columns
  (as fractions of image dimensions).

- H, W:

  dimensions of `saliency_map`. If `NULL`, taken from the matrix.

- match_radius_rt, match_radius_cv:

  match radius in pixels. Default 8 RT / 2 CV, matching the
  peak-detector's default `min_sep_rt`, `min_sep_cv`.

- top_k:

  number of local maxima to extract for precision/recall.

## Value

list(auc_per_pixel, precision_at_k, recall_at_k, matched_biomarkers,
missed_biomarkers, top_hotspots)

## Details

1.  **Per-pixel AUC**: label each pixel `1` if within
    `(match_radius_rt, match_radius_cv)` of any biomarker center, else
    `0`; score with the saliency value; report ROC AUC. This is "does
    saliency track biomarker regions on average?"

2.  **Top-K precision / recall**: extract the top-K non-max-suppressed
    local maxima of the saliency map; for each, check whether it falls
    within the match radius of any biomarker center. Precision = matched
    / K. Recall = unique biomarkers matched / n_biomarkers. This is
    "does saliency's hotspot list correspond to the biomarkers?"
