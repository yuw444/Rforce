# calculate the number of events in each interval

calculate the number of events in each interval

## Usage

``` r
counts_by_interval_and_id(
  event_times,
  ids,
  status,
  weights_by_status,
  interval
)
```

## Arguments

- event_times::

  a vector of event times

- ids::

  a vector of ids

- status::

  a vector of status

- weights_by_status::

  a vector of weights by status

- interval::

  a vector of intervals

## Value

a matrix of number of events in each interval by id, row is id, column
is interval
