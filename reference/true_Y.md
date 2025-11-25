# Calculate the mean number events for each subject at different time points in the simulated dataset

Calculate the mean number events for each subject at different time
points in the simulated dataset

## Usage

``` r
true_Y(compo_sim_list, time_points)
```

## Arguments

- compo_sim_list:

  the object from the `compo_sim` or `compo_sim_mao`

- time_points:

  the time points that need to be evaluated

## Value

a data frame (n_subject x length(time_point))
