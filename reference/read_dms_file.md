# Read one raw GC-DMS data file

Expects a tab-delimited file with this layout:

- Line 1: literal "Vc"

- Line 2: tab + compensation-voltage (CV) values, tab-separated

- Line 3: header row (typically "Time Stamp" followed by per-CV labels)

- Lines 4+: each row begins with a retention time (in seconds), followed
  by intensity values for each CV

## Usage

``` r
read_dms_file(path)
```

## Arguments

- path:

  File path to one GC-DMS file.

## Value

List with `path`, `time` (numeric vector), `cv` (numeric vector), and
`Z` (numeric matrix of intensity, rows = RT, cols = CV).

## Details

Implemented with base R only (no readr / stringr dependency) so it works
regardless of which packages the user has attached.
