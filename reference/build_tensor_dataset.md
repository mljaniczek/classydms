# Build a torch dataset from a list of intensity matrices and labels

Coerces matrices to a common (H, W) layout, stacks into a 4-D tensor of
shape `(N, 1, H, W)`, and returns a torch dataset whose items are
`list(x, y)` with one-channel images and integer labels.

## Usage

``` r
build_tensor_dataset(Z_list, y_vec)
```

## Arguments

- Z_list:

  List of 2-D numeric matrices, all of the same dimensions (or
  transposes that will be auto-corrected).

- y_vec:

  Integer or factor labels.

## Value

An instantiated torch dataset.
