# calculate observed risk time at the given interval

calculate observed risk time at the given interval

## Usage

``` r
observed_risk_time(Gt, X, Status = 0, low = 1, high = 3)
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

: the observed risk time at the given interval
