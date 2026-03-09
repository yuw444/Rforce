# Rforce

**Rforce** implements the methodology described in [Rforce: Random
Forests for Composite
Endpoints](https://onlinelibrary.wiley.com/doi/10.1002/sim.70413), which
models composite endpoints consisting of **non-fatal events and terminal
events**.

The method builds random forests using **generalized estimating
equations (GEE)** and handles dependent censoring caused by terminal
events using the concept of **pseudo-at-risk duration**.

This work received the **2024 Student Paper Competition Award** from the
American Statistical Association (ASA), jointly from the [Section on
Statistical Computing and Section on Statistical
Graphics](https://community.amstat.org/jointscsg-section/awards/student-paper-competition).

The paper is published in *Statistics in Medicine*:

- PMID: [41640374](https://pubmed.ncbi.nlm.nih.gov/41640374/)
- DOI:
  [10.1002/sim.70413](https://onlinelibrary.wiley.com/doi/10.1002/sim.70413)

The software provides both:

- **R API**
- **C API**

Key features include:

- High computational and memory efficiency
- Parallel computation using [OpenMP](https://www.openmp.org/)
- Reproducible results (see the reproducibility example
  [here](https://yuw444.github.io/Rforce/articles/get-started.html#in-1-equivalent-command))

------------------------------------------------------------------------

# Installation

## Dependencies

- `cmake >= 3.16.0` – build system for the C API
- `OpenMP` – parallel computing
- `R >= 4.3.3` – R interface

------------------------------------------------------------------------

## Install R API

``` r
# install.packages("devtools")
devtools::install_github("yuw444/Rforce")
```

### C API

``` bash
git clone https://github.com/yuw444/Rforce.git
cd Rforce
mkdir build
cd build
cmake ..
make
```

A `CMakeLists.txt` file is provided in the repository.

------------------------------------------------------------------------

## Usage

### R Examples

- Examples: [Get
  Started](https://yuw444.github.io/Rforce/articles/get-started.html).

------------------------------------------------------------------------

### Shell Scripts

``` bash
./Rforce [subcommands] <options>
```

### Available Subcommands

- `train` — Train a composite endpoint forest
- `predict` — Predict using a trained composite endpoint forest and new
  observations
- `-h, --help` — Show help message and exit

------------------------------------------------------------------------

## C API Subcommands

### Train

Train a composite endpoint forest model.

``` bash
Rforce train <options>
```

**Options:**

| Option                      | Description                                                                                | Required/Optional | Default                                   |
|-----------------------------|--------------------------------------------------------------------------------------------|-------------------|-------------------------------------------|
| `-d, --designMatrixY=<str>` | Path to design matrix                                                                      | **Required**      |                                           |
| `-a, --auxiliary=<str>`     | Path to auxiliary features                                                                 | **Required**      |                                           |
| `-u, --unitsOfCPIU=<str>`   | Path to unitsOfCPIU file                                                                   | **Required**      |                                           |
| `-o, --out=<str>`           | Path to output directory                                                                   | Optional          | Current working directory                 |
| `-v, --verbose=<int>`       | Verbosity level (0–3)                                                                      | Optional          | 0                                         |
| `-m, --maxDepth=<int>`      | Maximum tree depth                                                                         | Optional          | 10                                        |
| `-n, --minNodeSize=<int>`   | Minimum node size                                                                          | Optional          | 2 × len(unitsOfCPIU) - 1                  |
| `-g, --gain=<float>`        | Minimum gain for split                                                                     | Optional          | 0.0 (likelihood-based) or 1.3 (GEE-based) |
| `-t, --mtry=<int>`          | Number of variables to try during splitting                                                | Optional          | √(number of variables)                    |
| `-s, --nsplits=<int>`       | Number of splits to try per variable                                                       | Optional          | 10                                        |
| `-r, --nTrees=<int>`        | Number of trees                                                                            | Optional          | 200                                       |
| `-e, --seed=<int>`          | Random seed                                                                                | Optional          | 926                                       |
| `-p, --nPerms=<int>`        | Number of permutations for variable importance                                             | Optional          | 10                                        |
| `-u, --nVars=<int>`         | Number of variables in the design matrix                                                   | Optional          | Number of columns                         |
| `-i, --pathVarIds=<str>`    | Variable IDs (categorical variables supported via repeated IDs)                            | Optional          |                                           |
| `-x, --iDot`                | Output tree DOT files                                                                      | Optional          | False                                     |
| `-k, --k=<int>`             | Bayesian estimator parameter for leaf output                                               | Optional          | 4                                         |
| `-L, --long`                | Use multiple rows per patient (RF-SLAM style)                                              | Optional          |                                           |
| `-N, --nopseudo`            | Do not estimate pseudo risk time                                                           | Optional          |                                           |
| `-P, --pseudorisk1`         | Use original pseudo-risk time (population level)                                           | Optional          |                                           |
| `-B, --pseudorisk2`         | Recalculate pseudo-risk time at each tree (default)                                        | Optional          |                                           |
| `-D, --dynamicrisk`         | Dynamically estimate pseudo-risk time at each split                                        | Optional          |                                           |
| `-F, --nophi`               | Fix φ = 1, do not estimate φ                                                               | Optional          |                                           |
| `-P, --phi1`                | Estimate φ at population level                                                             | Optional          |                                           |
| `-H, --phi2`                | Estimate φ at tree level (default)                                                         | Optional          |                                           |
| `-Y, --dynamicphi`          | Dynamically estimate φ at each split                                                       | Optional          |                                           |
| `-G, --gee`                 | Use GEE approach                                                                           | Optional          |                                           |
| `-A, --padjust=<str>`       | p-value adjustment method (`bonferroni`, `holm`, `hochberg`, `hommel`, `BH`, `BY`, `none`) | Optional          | `BH`                                      |
| `-I, --interaction`         | Add interaction terms for GEE                                                              | Optional          | NULL                                      |
| `-S, --asym`                | Use asymptotic approach                                                                    | Optional          |                                           |
| `-T, --threads=<int>`       | Number of parallel computing threads                                                       | Optional          | 8                                         |

------------------------------------------------------------------------

### Predict

Predict using a trained model and test data.

``` bash
Rforce predict <options>
```

**Options:**

| Option              | Description              | Required/Optional | Default                   |
|---------------------|--------------------------|-------------------|---------------------------|
| `-m, --model=<str>` | Path to trained model    | **Required**      |                           |
| `-t, --test=<str>`  | Path to test data        | **Required**      |                           |
| `-o, --out=<str>`   | Path to output directory | Optional          | Current working directory |

------------------------------------------------------------------------

## Examples

Train a model:

``` bash
Rforce train -d design_matrix.csv -a auxiliary_features.csv -u unitsOfCPIU.txt -o output_folder -v 1
```

Predict with a trained model:

``` bash
Rforce predict -m output_folder/model.rforce -t test_data.csv -o prediction_results/
```

------------------------------------------------------------------------

## Notes

- By default, pseudo-risk time and φ (phi) are re-estimated at each tree
  level.
- Dynamic options (`--dynamicrisk`, `--dynamicphi`) allow estimates at
  each split for more flexibility.
- Parallel computation is supported via the `--threads` option.
- GEE-based splitting with p-value adjustment is available.
- An R API is currently actively developing which includes:
  - Classical survival data generation
  - Composite endpoint data generation
  - [`Wcompo`](https://cran.r-project.org/web/packages/Wcompo/index.html)
    methodology realization
  - An R interface to **Rforce**
