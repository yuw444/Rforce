# Save an Rforce Object to Disk

`saveRforce()` saves a trained `Rforce` object to disk, including both
the R-side metadata (e.g., model settings, variable importance,
parameters) and the full C-based forest structure. The external C
pointer is handled separately and excluded from the saved `.rds` file to
ensure portability and session safety.

## Usage

``` r
saveRforce(forest, file, ...)
```

## Arguments

- forest:

  An object of class `"Rforce"` representing the fitted random survival
  forest.

- file:

  A character string specifying the directory path where the object
  should be saved. The directory will be created by the underlying C
  routine if it does not already exist.

- ...:

  Additional arguments passed to
  [`saveRDS()`](https://rdrr.io/r/base/readRDS.html), such as `compress`
  or `version`.

## Value

Invisibly returns `NULL`. A message is not printed, but files are
written to the specified directory.

## Details

The function performs the following actions:

1.  Calls the native routine `R_SaveRforce` to serialize and write the
    forest structure (all trees) to the specified directory.

2.  Removes the external C pointer (`_external_forest_C_Ptr`) from the R
    object to prevent session-specific pointer issues.

3.  Saves the remaining R-side object to `Rforce.rds` using
    [`saveRDS()`](https://rdrr.io/r/base/readRDS.html).

The resulting directory structure is compatible with
[`loadRforce()`](https://yuw444.github.io/Rforce/reference/loadRforce.md),
which restores both the R object and its native forest pointer.

## See also

[`loadRforce`](https://yuw444.github.io/Rforce/reference/loadRforce.md)
to reload a saved Rforce object.

## Examples

``` r
if (FALSE) { # \dontrun{
# Fit an Rforce model and then save it
saveRforce(rf_model, file = "saved_model")
} # }
```
