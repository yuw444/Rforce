# Load an Rforce Object from Disk

`loadRforce()` restores a previously saved `Rforce` object from a
specified directory. It reconstructs both the R-side object (stored as
`Rforce.rds`) and the associated C++ forest structure via the
`R_LoadRforce` routine, ensuring full functionality.

## Usage

``` r
loadRforce(path)
```

## Arguments

- path:

  A character string specifying the directory containing the saved
  `Rforce` object. The directory must include `Rforce.rds` and a single
  forest data subdirectory created by
  [`saveRforce()`](https://yuw444.github.io/Rforce/reference/saveRforce.md).

## Value

A fully functional `Rforce` object with a reconstructed external forest
pointer, ready for prediction, visualization, or further analysis.

## Details

The function performs the following tasks:

1.  Reads the `Rforce.rds` file to reconstruct the R-side object.

2.  Validates that the restored object inherits class `"Rforce"`.

3.  Locates the associated forest data directory.

4.  Calls the native routine `R_LoadRforce` to rebuild the C-allocated
    forest structure, restoring the external pointer at
    `_external_forest_C_Ptr`.

If the directory structure is invalid or the C pointer fails to load, an
error is raised.

## Examples

``` r
if (FALSE) { # \dontrun{
# Load a previously saved Rforce model
rf <- loadRforce("path/to/saved/model")
rf
} # }
```
