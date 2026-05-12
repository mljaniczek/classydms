# Compute a Class Activation Map (CAM) for one sample

Weights each encoder channel's spatial activation by its elastic net
coefficient (`coefs`), summing across channels to produce a 2-D saliency
map. Because the elastic net classifier head is linear in the
global-average-pooled features, this construction is mathematically
equivalent to Class Activation Mapping using elastic net coefficients as
head weights (and similar in spirit to Grad-CAM, without requiring
backpropagation).

## Usage

``` r
compute_saliency_map(encoder, Z, coefs, device = "cpu")
```

## Arguments

- encoder:

  A pre-trained `dms_encoder`.

- Z:

  One sample's preprocessed, padded matrix.

- coefs:

  Numeric vector of length 64 (elastic net coefficients).

- device:

  "cpu" or "cuda".

## Value

List with `raw_cam` and `upsampled_cam`.

## Details

Returns the saliency at the encoder output resolution (raw_cam) and
bilinearly upsampled to the input dimensions (upsampled_cam).
