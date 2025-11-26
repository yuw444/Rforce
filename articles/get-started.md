# Get Started

## Generating Composite Endpoint Data

- `Rforce` package provides a direct function to generate composite
  endpoints survival data with options
  - Constant baseline hazard
  - Non-constant() baseline hazard
  - Proportional hazard
  - Non-proportional hazard
  - see details for these options in
    `Articles/Generating Composite Endpoints Survival Data`
- By default, we genereate a non-constant baseline hazard with

``` r
library(Rforce)
data_list <- compo_sim(n_patients = 500, verbose = FALSE)
```
