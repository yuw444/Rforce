# Fit a random forest for composite endpoints using CPIU-wide data

`Rforce()` fits a random forest model for composite endpoints—including
recurrent non-fatal events and a terminal event—using the CPIU
(composite pseudo–at-risk interval)–wide representation as the working
data structure.

The function can either:

- take raw patient-level long-formatted `data` and a
  `Surv(id, time, status) ~ covariates` `formula`, internally construct
  the CPIU-wide design matrix and auxiliary features, and then fit the
  forest; or

- directly accept a pre-computed `cpius` object (typically created by
  [`patients_to_cpius`](https://yuw444.github.io/Rforce/reference/patients_to_cpius.md)
  and
  [`cpius_to_dummy`](https://yuw444.github.io/Rforce/reference/cpius_to_dummy.md))
  that already contains the CPIU-wide matrices and meta-data.

## Usage

``` r
Rforce(
  data = NULL,
  formula = NULL,
  n_intervals = 10,
  weights_by_status = NULL,
  cpius = NULL,
  split_rule = c("Rforce-QLR", "Rforce-GEE", "Rforce-GEEIntr", "RF-SLAM"),
  n_trees = 300,
  mtry = NULL,
  n_splits = NULL,
  min_gain = NULL,
  min_node_size = NULL,
  max_depth = 8,
  seed = 926
)
```

## Arguments

- data:

  Optional `data.frame` containing patient-level data. Must include at
  least the patient ID, event time, and event status columns referenced
  on the left-hand side of `formula`. Ignored if `cpius` is supplied.

- formula:

  A [`Surv`](https://rdrr.io/pkg/survival/man/Surv.html) formula of the
  form `Surv(id, time, status) ~ covariate_1 + covariate_2 + ...`. The
  left-hand side must contain three variables in order: patient ID,
  follow-up time, and event status. If a character string is provided,
  it is internally converted via `as.formula`. Used only when `cpius` is
  `NULL`.

- n_intervals:

  Integer; number of follow-up intervals used to build the CPIU-wide
  representation when `cpius` is not supplied. Defaults to `10`. The
  default break points are based on empirical quantiles of the follow-up
  time.

- weights_by_status:

  Optional numeric vector of non-negative weights, one per unique event
  status in the original patient-level data (including censoring). If
  supplied, its length must match the number of unique values in the
  status variable. If `NULL`, the default is weight 0 for censoring
  (status = 0) and weight 1 for all event types.

- cpius:

  Optional list containing pre-computed CPIU-wide data, typically the
  output of
  [`patients_to_cpius`](https://yuw444.github.io/Rforce/reference/patients_to_cpius.md)
  (optionally processed by
  [`cpius_to_dummy`](https://yuw444.github.io/Rforce/reference/cpius_to_dummy.md)).
  It is expected to contain at least:

  - `design_matrix_Y`: numeric matrix of CPIU-wide outcomes and
    covariates;

  - `auxiliary_features`: numeric matrix of auxiliary CPIU-level
    features (e.g., pseudo–at-risk durations, event indicators, weights)
    aligned with `design_matrix_Y`;

  - `variableIds`: integer vector of group IDs for dummy-encoded
    covariates;

  - `unitsOfCPIUs`: numeric vector giving interval widths;

  - `isDummy`: logical flag indicating whether covariates are already
    dummy-encoded.

  If `cpius` is supplied, the internal data-to-CPIU conversion step is
  skipped.

- split_rule:

  Character string specifying the splitting rule:

  - `"Rforce-QLR"` — quasi-likelihood ratio–based splitting;

  - `"Rforce-GEE"` — GEE-based splitting without interaction;

  - `"Rforce-GEEIntr"` — GEE-based splitting with interaction;

  - `"RF-SLAM"` — RF-SLAM–style splitting without pseudo–at-risk
    durations.

  Partial matching is not supported; the value must match one of these
  strings exactly.

- n_trees:

  Integer; number of trees in the forest. Default is `300`.

- mtry:

  Integer; number of candidate variables drawn at random at each split.
  If `NULL`, a default is chosen internally based on the number of
  covariates.

- n_splits:

  Integer; number of candidate split points considered per variable at
  each node. If `NULL`, a default is chosen internally.

- min_gain:

  Numeric; minimum improvement in the split criterion required to accept
  a split. For `"Rforce-QLR"` and `"RF-SLAM"` this is interpreted as the
  percentage of quasi-likelihood increase (e.g. default `0.01`
  corresponds to 1\\ `"Rforce-GEEIntr"`, it is treated as a threshold on
  an adjusted \\p\\-value (e.g. default 0.2).

- min_node_size:

  Integer; minimum number of patients required in a terminal node (or to
  attempt further splitting). If `NULL`, a default is chosen internally
  based on the sample size.

- max_depth:

  Integer; maximum depth of each tree. Use `NULL` or a large value to
  allow essentially full growth, subject to `min_node_size` and
  `min_gain`. The default is `8`.

- seed:

  Integer random seed used to initialize the C back end for reproducible
  forests. Defaults to `926`.

## Value

An object of class `"Rforce"`, a list with elements including (but not
limited to):

- `data`, `formula` — the original data and model formula used to fit
  the forest (if available);

- `cpius` — the CPIU-wide list containing `design_matrix_Y`,
  `auxiliary_features`, `variableIds`, `unitsOfCPIUs`, and related
  meta-data;

- `bagMatrix` — bootstrap sample indicators for each tree;

- `treePhi` — estimates of the mean structure for each interval in each
  tree;

- `forestMatrix` — a compact summary of the forest structure in matrix
  form (not intended to be interpreted directly by users);

- `vimpStat` — forest structure and variable importance summaries
  returned from the C back end;

- `predicted`, `oobPredicted` — fitted and out-of-bag predictions on the
  CPIU scale;

- hyperparameters such as `nTrees`, `maxDepth`, `minNodeSize`,
  `minGain`, `mtry`, `nsplits`, and `seed`;

- `_external_forest_C_Ptr` — an external pointer storing the forest from
  the Rforce C API. It is valid only for the current R session unless
  serialized and restored via
  [`saveRforce`](https://yuw444.github.io/Rforce/reference/saveRforce.md)
  and
  [`loadRforce`](https://yuw444.github.io/Rforce/reference/loadRforce.md).

## Details

When `cpius` is `NULL`, `Rforce()`:

1.  checks that `data` and `formula` are supplied and that the left-hand
    side of `formula` is of the form `Surv(id, time, status)`;

2.  validates that all variables appearing on the left- and right-hand
    sides of `formula` are present in `data`;

3.  computes CPIU interval widths from the empirical distribution of the
    follow-up time using `n_intervals` empirical quantiles;

4.  constructs a long-format working dataset combining covariates and
    the `Id`, `X`, and `Status` response variables;

5.  calls
    [`patients_to_cpius`](https://yuw444.github.io/Rforce/reference/patients_to_cpius.md)
    to build a CPIU-wide design matrix and auxiliary features, with
    pseudo–at-risk durations used for `"Rforce-QLR"`, `"Rforce-GEE"`,
    and `"Rforce-GEEIntr"`, and disabled for `"RF-SLAM"`;

6.  applies
    [`cpius_to_dummy`](https://yuw444.github.io/Rforce/reference/cpius_to_dummy.md)
    to dummy-encode factor and character covariates if they are not
    already dummy-encoded;

7.  passes the resulting CPIU-wide matrices, variable grouping, and
    interval widths to the C routine `R_Rforce`, together with the
    specified splitting rule and hyperparameters.

The fitted forest is returned as an S3 object of class `"Rforce"` with
components including the CPIU-wide data structure, forest structure, and
out-of-bag predictions, along with various meta-data to facilitate
downstream methods such as `print.Rforce()`, `summary.Rforce()`, and
[`predict.Rforce()`](https://yuw444.github.io/Rforce/reference/predict.Rforce.md).

## See also

[`compo_sim`](https://yuw444.github.io/Rforce/reference/compo_sim.md),
[`random_censoring`](https://yuw444.github.io/Rforce/reference/random_censoring.md),
[`patients_to_cpius`](https://yuw444.github.io/Rforce/reference/patients_to_cpius.md),
[`cpius_to_dummy`](https://yuw444.github.io/Rforce/reference/cpius_to_dummy.md),
[`Surv`](https://rdrr.io/pkg/survival/man/Surv.html),
[`predict.Rforce`](https://yuw444.github.io/Rforce/reference/predict.Rforce.md),
`summary.Rforce`,
[`vimp.Rforce`](https://yuw444.github.io/Rforce/reference/vimp.Rforce.md),
[`saveRforce`](https://yuw444.github.io/Rforce/reference/saveRforce.md),
[`loadRforce`](https://yuw444.github.io/Rforce/reference/loadRforce.md).

## Examples

``` r
## Example using simulated data
library(Rforce)

set.seed(926)

data_list <- compo_sim()

df_censored <- random_censoring(
  data       = data_list$dataset,
  event_rate = 0.5
)

fit <- Rforce(
  data         = df_censored,
  formula      = Surv(Id, X, Status) ~ .,
  n_intervals  = 5,
  n_trees      = 50,
  max_depth    = 10,
  min_node_size = 20,
  min_gain     = 0,
  mtry         = 2,
  n_splits     = 5
)
#> Error: The data frame must contain the column: .

fit
#> Error: object 'fit' not found
```
