# Dummy-encode factor/character covariates in CPIU object

This function takes a CPIU object and dummy-encodes any factor or
character covariates in the design matrix for the outcome variable Y. It
removes the reference level for each factor to avoid multicollinearity.

## Usage

``` r
cpius_to_dummy(object, cols_to_dummy = NULL)
```

## Arguments

- object:

  A CPIU object.

- cols_to_dummy:

  Optional vector of column names to dummy-encode. If NULL, all
  factor/character columns will be dummy-encoded.

## Value

A CPIU object with dummy-encoded covariates in the `designMatrixY`.
