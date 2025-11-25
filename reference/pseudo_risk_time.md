# calculate pseudo risk time at the given interval

calculate pseudo risk time at the given interval

## Usage

``` r
pseudo_risk_time(Gt, X, Status = 0, low = 1, high = 3)
```

## Arguments

- Gt::

  a list of unique time and survival probability of Censoring

- X::

  the time point associated with the Status

- Status::

  the event status at the time point X

- low::

  the lower bound of the interval

- high::

  the upper bound of the interval

## Value

: the pseudo risk time at the given interval

## Examples

``` r
# example code
 Gt <- list(unique_time = c(0, 1, 1.8, 3.5, 4),
            surv_prop = c(1, 0.8,0.6, 0.3, 0.1))
pseudo_risk_time(Gt, X =1.5 , Status = 1, low =0.5, high=1.9)
#> Error in pseudo_risk_time(Gt, X = 1.5, Status = 1, low = 0.5, high = 1.9): could not find function "pseudo_risk_time"
```
