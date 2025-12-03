# Dummy-encode factor/character covariates in CPIU object

This function takes a CPIU object and dummy-encodes any factor or
character covariates in the design matrix for the outcome variable Y. It
removes the reference level for each factor to avoid multicollinearity.

## Usage

``` r
cpius_to_dummy(object, cols_to_keep = NULL)
```

## Arguments

- object:

  A CPIU object.

- cols_to_keep:

  Optional vector of column names to keep and consider for dummy
  encoding. If NULL, all columns will be considered for dummy-encoded.

## Value

A CPIU object with dummy-encoded covariates in the `designMatrixY`.
