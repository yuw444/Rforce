# apply administrative censoring in the given dataset

censoring the given dataset based on the given quantile of observed time

## Usage

``` r
admin_censoring(data_to_convert, censoring_quantile)
```

## Arguments

- data_to_convert:

  the dataset to be converted; including the columns of `Id`, `Status`
  and `Time`

- censoring_quantile:

  the quantile of observed time to be censored

## Value

the censored dataset

## Examples

``` r
list_data_to_convert <- compo_sim()
df_converted <- manual_censoring(
    list_data_to_convert$dataset,
    0.9
)
#> Error in manual_censoring(list_data_to_convert$dataset, 0.9): could not find function "manual_censoring"
str(df_converted)
#> Error: object 'df_converted' not found
```
