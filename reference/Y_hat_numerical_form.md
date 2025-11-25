# Numerical form to calculate the true Y given the true hazard and other parameters

Numerical form to calculate the true Y given the true hazard and other
parameters

## Usage

``` r
Y_hat_numerical_form(
  t,
  constant_baseline_hazard,
  baseline_hazard = 1,
  mu_0,
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

- sigma_scale_gamma:

  the parameter of gamma distribution when simulate the frality term

- lambdaZ:

  a vector, recurrent event hazard rate

- lambda:

  a scalar, the hazard rate of the stopping time

- a_shape_weibull:

  the parameter of weibull distribution when simulate the data

- sigma_scale_weibull:

  the parameter of weibull distribution when simulate the data

## Value

a scalar, the mean number of event at time `t`
