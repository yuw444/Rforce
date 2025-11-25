# Empirical Mean Number of Events at the Different Time Point

Empirical Mean Number of Events at the Different Time Point

## Usage

``` r
empirical_Y(
  compo_sim_list,
  weibull_baseline = TRUE,
  x,
  time_points,
  n_sims,
  n_size
)
```

## Arguments

- compo_sim_list:

  simulation list either from `compo_sim` or `compo_sim_mao`

- weibull_baseline:

  which simulation scheme is used, default is `compo_sim`

- x:

  patient characteristics to exam

- time_points:

  time point to exam

- n_sims:

  number of simulations to run

- n_size:

  number of patients at each simulation

## Details

Calculate the empirical mean number of event at the different time
points by using data simulation
