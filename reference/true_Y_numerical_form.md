# Numerical form to calculate the true Y given the true hazard and other parameters

Numerical form to calculate the true Y given the true hazard and other
parameters

## Usage

``` r
true_Y_numerical_form(
  t,
  constant_baseline_hazard = FALSE,
  baseline_hazard = 1,
  a_shape_weibull,
  sigma_scale_weibull,
  sigma_scale_gamma,
  lambdaZ,
  lambda
)
```

## Arguments

- t:

  the time point to check

- constant_baseline_hazard:

  logical; whether constant baseline hazard is used

- baseline_hazard:

  a scalar, the constant baseline hazard

- a_shape_weibull:

  the parameter of weibull distribution when simulate the data

- sigma_scale_weibull:

  the parameter of weibull distribution when simulate the data

- sigma_scale_gamma:

  the parameter of gamma distribution when simulate the frality term

- lambdaZ:

  a vector, recurrent event hazard rate

- lambda:

  a scalar, the hazard rate of the stopping time

## Value

a scalar, the mean number of event at time `t`
