# Calculate the predicted number of events at given time points

Calculate the predicted number of events at given time points

## Usage

``` r
add_Y_hat(lambda_pred, length_cpius, time_to_evaluate)
```

## Arguments

- lambda_pred:

  a matrix of predicted hazard rates at each interval for multiple
  subjects

- time_to_evaluate:

  a scalar, time point to evaluate \\\hat Y\\

- interval_cpius:

  a vector of lengths for each CPIU, same length with
  `ncol(lambda_pred)`

## Value

a matrix of predicted number of events at each `time_points` for
multiple subjects

## Examples

``` r
lambdas <- matrix(1:10, nrow = 2, byrow = TRUE)
intervals <- rep(3,5)
add_Y_hat(lambdas, intervals, c(3, 6, 9))
#> Error: length(time_to_evaluate) not equal to 1
```
