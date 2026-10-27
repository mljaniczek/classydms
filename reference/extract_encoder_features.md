# Extract encoder features for all samples

Runs each sample through the (frozen) encoder and reduces its feature
map to a vector, stacking the results into a matrix suitable for a
classical classifier.

## Usage

``` r
extract_encoder_features(
  encoder,
  Z_list,
  device = "cpu",
  readout = c("gap", "presence"),
  presence_head = NULL
)
```

## Arguments

- encoder:

  A pre-trained `dms_encoder`.

- Z_list:

  List of preprocessed, padded intensity matrices.

- device:

  "cpu" or "cuda".

- readout:

  `"gap"` or `"presence"`.

- presence_head:

  Required for `readout = "presence"`.

## Value

Numeric matrix with one row per sample and named columns.

## Details

Readouts:

- `"gap"` (default): global average pooling -\> `(N, latent_dim)`,
  columns `encoder_ch_1..`.

- `"presence"`: per-compound presence logits from a pretrained
  [`dms_presence_head()`](https://mljaniczek.github.io/classydms/reference/dms_presence_head.md)
  -\> `(N, K)`, columns `cmp_<compound_id>`. Each column is evidence
  that one catalog compound is present, so an elastic net on these can
  sum evidence across biomarkers and its coefficients name compounds
  directly.
