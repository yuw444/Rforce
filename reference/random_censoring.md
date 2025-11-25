# apply random censoring in the given dataset

censoring the given dataset based on the given event rate

## Usage

``` r
random_censoring(data_to_convert, event_rate, seed = 926, verbose = FALSE)
```

## Arguments

- data_to_convert:

  the dataset to be converted; including the columns of `Id`, `Status`
  and `Time`

- event_rate:

  the desired event rate

- seed:

  seed for RNG, default is 926

- verbose:

  print out details

## Value

the censored dataset

## Examples

``` r
list_data_to_convert <- compo_sim()
df_converted <- random_censoring(
    list_data_to_convert$dataset,
    0.8
)
str(df_converted)
#> 'data.frame':    3197 obs. of  13 variables:
#>  $ Id          : int  1 2 3 3 3 3 4 5 6 7 ...
#>  $ X           : num  0.00123 0.54805 0.05301 0.10663 1.24427 ...
#>  $ Status      : num  1 1 2 2 2 1 1 0 1 2 ...
#>  $ binary1     : num  0 1 1 1 1 1 0 1 0 1 ...
#>  $ binary2     : num  1 0 1 1 1 1 0 0 1 1 ...
#>  $ binary3     : num  0 1 1 1 1 1 0 1 0 0 ...
#>  $ binary4     : num  0 0 1 1 1 1 0 1 0 0 ...
#>  $ binary5     : num  0 1 1 1 1 1 0 1 0 0 ...
#>  $ binary6     : num  0 0 1 1 1 1 0 0 0 0 ...
#>  $ continuous7 : num  0.228 0.798 0.463 0.463 0.463 ...
#>  $ continuous8 : num  0.568 0.714 0.0105 0.0105 0.0105 ...
#>  $ continuous9 : num  0.7125 0.2585 0.0927 0.0927 0.0927 ...
#>  $ continuous10: num  0.738 0.447 0.483 0.483 0.483 ...
```
