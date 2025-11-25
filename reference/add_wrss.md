# Calculate the WRSS

Calculate the WRSS

## Usage

``` r
add_wrss(df_test, weights_by_status, lambda_pred, length_cpius, time_point)
```

## Arguments

- df_test:

  a data frame contains composite event observations

- weights_by_status:

  a vector, the weights assign to each event type

- lambda_pred:

  a matrix contains the predicated hazard rate for each subject at each
  interval

- length_cpius:

  a vector, the length of each CPIU

- time_point:

  a scalar, the time point to evaluate WRSS

## Value

a data frame contains WRSS for each subject

## Details

Calculate the WRSS(Gerds T. 2006) to evaluate the performance of
`lambda_pred` when the true Y is unknown
