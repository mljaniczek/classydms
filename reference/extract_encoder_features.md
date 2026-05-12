# Extract 64-dim encoder features for all samples (GAP-pooled)

Runs each sample through the (frozen) encoder, applies global average
pooling to its feature map, and stacks the result into an `(N, 64)`
matrix suitable for use as input to a classical classifier.

## Usage

``` r
extract_encoder_features(encoder, Z_list, device = "cpu")
```

## Arguments

- encoder:

  A pre-trained `dms_encoder`.

- Z_list:

  List of preprocessed, padded intensity matrices.

- device:

  "cpu" or "cuda".

## Value

Numeric matrix with one row per sample, 64 columns.
