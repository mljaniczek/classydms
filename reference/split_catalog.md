# Split a catalog into two disjoint sub-catalogs by compound

Randomly partitions the catalog's compounds into two sets. Useful for
encoder-transfer experiments: pretrain the encoder on synthetic data
generated from sub-catalog A, then simulate case/control data from
sub-catalog B — biomarker compounds and background compounds in B are
all guaranteed to have been unseen during pretraining.

## Usage

``` r
split_catalog(catalog, frac = 0.5, seed = NULL)
```

## Arguments

- catalog:

  A `peak_catalog` object.

- frac:

  Fraction of compounds to place in the first split. Default `0.5`.

- seed:

  Optional RNG seed for reproducible partitioning.

## Value

A list with `a` and `b`, each a valid `peak_catalog` with its own
compound subset. `a` also carries any `small_clusters` and `singletons`
from the original.

## Details

Small clusters and singletons are left untouched (attached to the first
split only) since they're not typically used for synthetic generation.
