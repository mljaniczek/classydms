# Elastic net classification on 64-dim encoder features

Performs stratified k-fold CV: extract features once via the frozen
encoder, then fit `cv.glmnet` independently on each fold. Returns
out-of-fold predictions, AUC, the feature matrix, and elastic net
coefficients (which map directly to encoder channels and are the
starting point for spatial saliency interpretation).

## Usage

``` r
cv_encoder_elastic_net(
  encoder,
  Z_list,
  y,
  k = 5L,
  alpha = 0.5,
  nlambda = 100L,
  seed = 1L,
  device = "cpu"
)
```

## Arguments

- encoder:

  A pre-trained `dms_encoder`.

- Z_list:

  List of preprocessed, padded intensity matrices.

- y:

  Binary labels.

- k:

  Number of CV folds.

- alpha:

  Elastic net mixing (1 = LASSO, 0 = ridge).

- nlambda:

  Number of lambda values for `cv.glmnet`.

- seed:

  RNG seed for fold assignment.

- device:

  "cpu" or "cuda".

## Value

List with `probs`, `auc`, `features`, and `coefs`.
