# Aggregate saliency maps across samples, masking padded regions

Computes a single saliency map averaged over `sample_indices`. For each
sample, only the original (non-padded) region contributes; the divisor
at each pixel is the number of samples that actually had data there.
Prevents data/padding boundary artifacts from dominating the aggregate
map.

## Usage

``` r
aggregate_saliency_masked(
  encoder,
  Z_list,
  coefs,
  sample_indices,
  orig_dims,
  device = "cpu"
)
```

## Arguments

- encoder:

  Pre-trained `dms_encoder`.

- Z_list:

  List of padded matrices.

- coefs:

  Elastic net coefficients (length 64).

- sample_indices:

  Which samples to aggregate over.

- orig_dims:

  List of `c(H, W)` pairs giving each sample's pre-padding dimensions.

- device:

  "cpu" or "cuda".

## Value

Matrix of the same shape as the padded image.
