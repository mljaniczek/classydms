# Read all GC-DMS files in a directory and assemble a samples list

Convenience wrapper that lists files matching `pattern` in `disease_dir`
and `control_dir`, reads each with
[`read_dms_file()`](https://mljaniczek.github.io/classydms/reference/read_dms_file.md)
(skipping any that error), assigns each sample a name based on the file
basename so downstream results can be traced back to a specific subject,
and returns a parallel pair `samples` and `y`.

## Usage

``` r
load_dms_directory(
  disease_dir,
  control_dir,
  pattern = "POS.*\\.(xls|txt|tsv|csv)$",
  verbose = TRUE
)
```

## Arguments

- disease_dir, control_dir:

  Paths to directories.

- pattern:

  Regex matching the data files. Default selects positive-channel files;
  change to `"NEG.*\\.(xls|txt|tsv|csv)$"` for negative-channel-only.

- verbose:

  Whether to print progress messages.

## Value

List with `samples` (named list of sample lists) and `y` (named integer
vector with 2 = disease, 1 = control).

## Details

Default `pattern` selects the positive ion channel only (POS files),
which is what the classification pipeline trains on. The negative
channel is not used.

Sample names are set to the file basename with extension stripped (e.g.,
`"241115-S1994-S1-DMS1_POS"`). This means you can index a sample by name
later — e.g., `samples[["241115-S1994-S1-DMS1_POS"]]` — and so can `y`
if you set `names(y) <- names(samples)`. The function explicitly
verifies that no name is the empty string, which would otherwise
propagate through [`lapply()`](https://rdrr.io/r/base/lapply.html) and
cause "attempt to use zero-length variable name" errors downstream.
