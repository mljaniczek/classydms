# Learned-query attention pooling readout

Replacement for GAP in a downstream classifier: `n_queries` learned
queries attend over the encoder feature map (keys carry a fixed 2D
sinusoidal positional encoding, so queries can become location-
specific). Output is the concatenation of the `n_queries` pooled feature
vectors, `(B, n_queries * C)`. Trained end-to-end with the classifier;
keep `n_queries` small when labeled n is small.

## Usage

``` r
attention_pool_readout(in_channels, n_queries = 4L, d_attn = 32L)
```

## Arguments

- in_channels:

  Encoder output channels.

- n_queries:

  Number of learned queries.

- d_attn:

  Query/key dimension (multiple of 4).
