# get-started

## Rforce

**Composite Endpoint** data type that has gained momentum in the field,
due to its fully consideration on the all observated events for each
patients, in comparison to the **Competing Risk** setting that utilize
only on the first observed event for each patient.

A parametric approach, like **Wcompo** Mao and Lin (2016), has proposed
to addressed the inference challenge while imposing the proportional
hazard assumption. Inspired by **Random Survival Forest** Ishwaran et
al. (2008) and **Counting Process Information Unit(CPIU)** Wongvibulsin,
Wu, and Zeger (2019), we here porposed a non-parametric ensememble
method, **Random Froest for Composite Endpoints*(**Rforce**)*, to
further relax the proportional hazard assumption.

Below, we like to demonstrate a simulated data example to show off the
simplistic pipeline of **Rforce**.

### Generating Composite Endpoint Data

- `Rforce` package provides functions to generate composite endpoints
  survival data with options
  - Constant baseline hazard
  - Non-constant() baseline hazard
  - Proportional hazard
  - Non-proportional hazard
  - see details for these options in
    `Articles/Generating Composite Endpoints/Survival Data`
- By default, we generate a non-constant(Weibull) baseline hazard with
  proportional hazard assumption.

``` r
knitr::opts_chunk$set(echo = TRUE)

library(Rforce)
data_list <- compo_sim(n_patients = 500, verbose = FALSE)
data_list$true_beta
#>  [1] 0.0 0.0 0.0 0.6 0.0 0.0 0.8 0.0 0.7 0.0
data_list$hazard_function
#> function (x) 
#> {
#>     exp(x %*% true_beta)
#> }
#> <bytecode: 0x5557a74d9638>
#> <environment: 0x5557a74dcdd8>
head(data_list$dataset)
#>   Id        Time Status binary1 binary2 binary3 binary4 binary5 binary6
#> 1  1 0.003487709      2       0       1       0       0       1       0
#> 2  1 0.047510838      2       0       1       0       0       1       0
#> 3  1 1.646181137      1       0       1       0       0       1       0
#> 4  2 0.074786350      2       0       0       1       0       1       0
#> 5  2 0.266915324      2       0       0       1       0       1       0
#> 6  2 0.833510030      2       0       0       1       0       1       0
#>   continuous7 continuous8 continuous9 continuous10
#> 1   0.4147548   0.6333024   0.4278975   0.86040861
#> 2   0.4147548   0.6333024   0.4278975   0.86040861
#> 3   0.4147548   0.6333024   0.4278975   0.86040861
#> 4   0.6814544   0.5258878   0.6166269   0.03045993
#> 5   0.6814544   0.5258878   0.6166269   0.03045993
#> 6   0.6814544   0.5258878   0.6166269   0.03045993
```

- Each patient’s observation can span across multiple rows, each row
  representing one of event of interest, ordered by time of occurrence.
- For completeness, each patient has complete observations up to
  terminal(fatal) event.
- However,
  [`compo_sim()`](https://yuw444.github.io/Rforce/reference/compo_sim.md)
  only generates non time-varying covariates although the proposed
  methodology can handle time-varying covariates.
- Here, we generate 6 binary covariates and 4 continuous covariates by
  default.
- The hazard function consist of baseline hazard $\lambda_{0}(t)$ and
  parametric part $g(\mathbf{Z})$, i.e.,
  $$\lambda(t) = \lambda_{0}(t)g(\mathbf{Z})$$ where
  $\lambda_{0}(t) \sim Weibull(1,1)$ and
  $g(\mathbf{Z}) = \exp(\beta\mathbf{Z})$.
- Only two covariates have non-zero effects on the hazard function,
  $\beta$ is indicated by `data_list$true_beta`.

### Implementating Random Censoring

- In real world clinical trials, patients may be censored before
  experiencing terminal event.
- `Rforce` package provides
  [`random_censoring()`](https://yuw444.github.io/Rforce/reference/random_censoring.md)
  function to implement random censoring mechanism on the generated
  complete data.

``` r
df <- random_censoring(
  data = data_list$dataset,
  event_rate = 0.5
)
head(df)
#>   Id           X Status binary1 binary2 binary3 binary4 binary5 binary6
#> 1  1 0.003487709      2       0       1       0       0       1       0
#> 2  1 0.047510838      2       0       1       0       0       1       0
#> 3  1 1.646181137      1       0       1       0       0       1       0
#> 4  2 0.074786350      2       0       0       1       0       1       0
#> 5  2 0.266915324      2       0       0       1       0       1       0
#> 6  2 0.539628591      0       0       0       1       0       1       0
#>   continuous7 continuous8 continuous9 continuous10
#> 1   0.4147548   0.6333024   0.4278975   0.86040861
#> 2   0.4147548   0.6333024   0.4278975   0.86040861
#> 3   0.4147548   0.6333024   0.4278975   0.86040861
#> 4   0.6814544   0.5258878   0.6166269   0.03045993
#> 5   0.6814544   0.5258878   0.6166269   0.03045993
#> 6   0.6814544   0.5258878   0.6166269   0.03045993
```

### Step 1: Converting CPIU-wide Format

- Evidently, the hazard function $\lambda(t)$ is not a constant of time.
  Inspired by **Counting Process Information Unit(CPIU)**, we convert
  the observed data into *CPIU-wide* format to better capture the hazard
  pattern for the population.
- The interval of units `units_of_cpiu` is defined such that each
  interval contains approximately equal number of events, e.g. 10
  intervals.
- [`patients_to_cpius()`](https://yuw444.github.io/Rforce/reference/patients_to_cpius.md)
  is customized for these purpose.

``` r
units_of_cpiu <- diff(c(0, quantile(df$X, 1 / 10 * 1:10)))
cpiu_wide <- patients_to_cpius(
  data_to_convert = df,
  units_of_cpiu = units_of_cpiu
)
names(cpiu_wide)
#>  [1] "data"                   "unitsOfCPIUs"           "nIntervals"            
#>  [4] "eventTypes"             "eventWeights"           "pseudoRisk"            
#>  [7] "wideFormat"             "designMatrix_Y"         "auxiliaryFeatures"     
#> [10] "nPatients"              "variableNameOrignal"    "variableSummaryOrignal"
#> [13] "isDummy"                "variableUsed"           "variableIds"           
#> [16] "formula"
class(cpiu_wide)
#> [1] "CPIU"
```

- `CPIU` is S3 class defined in **Rforce** package to better encapsulate
  `CPIU-wide` formatted data, along with the metadata.
- Two main `data.frame`’s are included,
  - `designMatrix_Y`: covariates + number of events in each defined
    interval
  - `auxilaryFeatures`: `Id` + `pseudo risk time` in each defined
    interval + `X` + `Status`

### Step 2: Encoding Factor/Categorical Variables

- In real world, there may have factor or categorical variable in the
  data. Internally, **Rforce** encodes those variables into dummy
  variables, to avoid the collinearlity, the number of encoded dummy
  varaibles for each categorical `variable` always equal to
  `nlevels(variable) - 1`.

- Worthy to metion, **Rforce** don’t evaluate the variable importance on
  individual levels of these categorical variable. Through the carefully
  design permutation framework, **Rforce** direct output the variable
  importance on the variable itself after the model fitting.

``` r
cpiu_wide <- cpius_to_dummy(cpiu_wide)
```

### Step 3: Fit Random Forest Model

- [`Rforce()`](https://yuw444.github.io/Rforce/reference/Rforce.md) does
  the final model fitting
- Internally, there are various of settings for the hyperparameters for
  growing random forest with some smart default. For example,
  - `n_trees`
  - `mtry`
  - `n_splits`
  - `min_gain`
  - `min_node_size`
  - `max_depth`
  - `seed`.
  - see details for the default in `Reference/Rforce()`.
- Same time, the default `split_rule` is using the proposed
  quasi-likelihood ratio `Rforce-QLR`, it is more computational
  efficient than other two proposed rules,
  - `Rforce-GEE`
  - `Rforce-GEEInter`
  - see details for the default in `Reference/Rforce()`.

``` r
fit <- Rforce(
  cpius = cpiu_wide,
  split_rule = "Rforce-QLR"
)
#> Fitting Rforce-QLR forest...
```

### 

Ishwaran, Hemant, Udaya B. Kogalur, Eugene H. Blackstone, and Michael S.
Lauer. 2008. “Random Survival Forests.” *The Annals of Applied
Statistics* 2 (3). <https://doi.org/10.1214/08-aoas169>.

Mao, Lu, and DY Lin. 2016. “Semiparametric Regression for the Weighted
Composite Endpoint of Recurrent and Terminal Events.” *Biostatistics* 17
(2): 390–403.

Wongvibulsin, Shannon, Katherine C. Wu, and Scott L. Zeger. 2019.
“Clinical Risk Prediction with Random Forests for Survival,
Longitudinal, and Multivariate (RF-SLAM) Data Analysis.” *BMC Medical
Research Methodology* 20 (December).
<https://doi.org/10.1186/s12874-019-0863-0>.
