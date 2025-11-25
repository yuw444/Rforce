# convert the recorded event time per patient to the number of events and wt at the given time point per patient

convert the recorded event time per patient to the number of events and
wt at the given time point per patient

## Usage

``` r
add_wt(data_to_convert, weights_by_status, time_to_evaluate)
```

## Arguments

- time_to_evaluate:

  a scalar, or vector with length equal to
  `length(unique(data_to_convert$Id))` the time point to calculate the
  number of events and wt for each patient;

- data_to_convert::

  a data frame with columns of Id, X, Status

- weights_by_status::

  a vector of weights by status, default is `c(0,1,1)` for censoring
  (status: 0), terminal event(status: 1) and recurrent event(status: 1)

## Value

a data.frame with number of events and wt at the given time point per
patient
