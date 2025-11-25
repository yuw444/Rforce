# convert the recorded event time per patient to the number of events per interval per patient

convert the recorded event time per patient to the number of events per
interval per patient

## Usage

``` r
patients_to_cpius(
  data_to_convert,
  units_of_cpiu,
  weights_by_status = c(0, 1, 1),
  pseudo_risk = TRUE,
  wide_format = TRUE
)
```

## Arguments

- data_to_convert::

  a data frame with columns of Id, X, Status

- units_of_cpiu::

  a vector of units of CPIU

- weights_by_status::

  a vector of weights by status, default is c(0,1,1) for censoring
  (status: 0), terminal event(status: 1) and recurrent event(status: 1)

- pseudo_risk::

  a boolean value indicating whether to use pseudo risk time

## Value

: a dataframe with number of events in each interval by id, row is id,
column is interval
