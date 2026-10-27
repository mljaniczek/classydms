# Stratified K-fold partition of binary labels

Always pass `seed` when folds must be reproducible. torch's lazy
first-use initialisation (e.g. the first `torch_load()` or tensor
creation in a session) consumes draws from R's RNG, so calling
[`set.seed()`](https://rdrr.io/r/base/Random.html) before torch work
does not pin folds drawn afterwards from the global stream.

## Usage

``` r
make_stratified_folds(y, k = 5, seed = NULL)
```
