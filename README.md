# Rforce

We introduce random forests for composite endpoints(**Rforce**) consisting of non-fatal composite events and terminal events. It utilizes generalized estimating equations to build trees and handles the dependent censoring due to the terminal events with the concept of pseudo-at-risk duration.

## Installation

(Instructions for installation if needed, or write something like:  
"Download the binary and add it to your PATH.")

## Usage

```bash
Rforce [subcommands] <options>
```

### Available Subcommands

- `train` — Train a composite endpoint forest
- `predict` — Predict using a trained composite endpoint forest and new observations
- `-h, --help` — Show help message and exit

---

## Subcommand Details

### Train

Train a composite endpoint forest model.

```bash
Rforce train <options>
```

**Options:**

| Option | Description | Required/Optional | Default |
|:------|:------------|:------------------|:--------|
| `-d, --designMatrixY=<str>` | Path to design matrix | **Required** | |
| `-a, --auxiliary=<str>` | Path to auxiliary features | **Required** | |
| `-u, --unitsOfCPIU=<str>` | Path to unitsOfCPIU file | **Required** | |
| `-o, --out=<str>` | Path to output directory | Optional | Current working directory |
| `-v, --verbose=<int>` | Verbosity level (0–3) | Optional | 0 |
| `-m, --maxDepth=<int>` | Maximum tree depth | Optional | 10 |
| `-n, --minNodeSize=<int>` | Minimum node size | Optional | 2 × len(unitsOfCPIU) - 1 |
| `-g, --gain=<float>` | Minimum gain for split | Optional | 0.0 (likelihood-based) or 1.3 (GEE-based) |
| `-t, --mtry=<int>` | Number of variables to try during splitting | Optional | √(number of variables) |
| `-s, --nsplits=<int>` | Number of splits to try per variable | Optional | 10 |
| `-r, --nTrees=<int>` | Number of trees | Optional | 200 |
| `-e, --seed=<int>` | Random seed | Optional | 926 |
| `-p, --nPerms=<int>` | Number of permutations for variable importance | Optional | 10 |
| `-u, --nVars=<int>` | Number of variables in the design matrix | Optional | Number of columns |
| `-i, --pathVarIds=<str>` | Variable IDs (categorical variables supported via repeated IDs) | Optional | |
| `-x, --iDot` | Output tree DOT files | Optional | False |
| `-k, --k=<int>` | Bayesian estimator parameter for leaf output | Optional | 4 |
| `-L, --long` | Use multiple rows per patient (RF-SLAM style) | Optional | |
| `-N, --nopseudo` | Do not estimate pseudo risk time | Optional | |
| `-P, --pseudorisk1` | Use original pseudo-risk time (population level) | Optional | |
| `-B, --pseudorisk2` | Recalculate pseudo-risk time at each tree (default) | Optional | |
| `-D, --dynamicrisk` | Dynamically estimate pseudo-risk time at each split | Optional | |
| `-F, --nophi` | Fix φ = 1, do not estimate φ | Optional | |
| `-P, --phi1` | Estimate φ at population level | Optional | |
| `-H, --phi2` | Estimate φ at tree level (default) | Optional | |
| `-Y, --dynamicphi` | Dynamically estimate φ at each split | Optional | |
| `-G, --gee` | Use GEE approach | Optional | |
| `-A, --padjust=<str>` | p-value adjustment method (`bonferroni`, `holm`, `hochberg`, `hommel`, `BH`, `BY`, `none`) | Optional | `BH` |
| `-I, --interaction` | Add interaction terms for GEE | Optional | NULL |
| `-S, --asym` | Use asymptotic approach | Optional | |
| `-T, --threads=<int>` | Number of parallel computing threads | Optional | 8 |

---

### Predict

Predict using a trained model and test data.

```bash
Rforce predict <options>
```

**Options:**

| Option | Description | Required/Optional | Default |
|:------|:------------|:------------------|:--------|
| `-m, --model=<str>` | Path to trained model | **Required** | |
| `-t, --test=<str>` | Path to test data | **Required** | |
| `-o, --out=<str>` | Path to output directory | Optional | Current working directory |

---

## Examples

Train a model:

```bash
Rforce train -d design_matrix.csv -a auxiliary_features.csv -u unitsOfCPIU.txt -o output_folder -v 1
```

Predict with a trained model:

```bash
Rforce predict -m output_folder/model.rforce -t test_data.csv -o prediction_results/
```

---

## Notes

- By default, pseudo-risk time and φ (phi) are re-estimated at each tree level.
- Dynamic options (`--dynamicrisk`, `--dynamicphi`) allow estimates at each split for more flexibility.
- Parallel computation is supported via the `--threads` option.
- GEE-based splitting with p-value adjustment is available.

