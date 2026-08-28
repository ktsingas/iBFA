# iBFA

**interval Bayesian Factor Analysis** for incomplete multi-omics data.

> ktsingas · github.com/ktsingas/iBFA

---

<!-- Insert graphical abstract here once ready -->
<!-- ![](manuscript/graphical_abstract.pdf) -->

---

## What it does

iBFA learns shared latent structure across mixed data types (continuous, binary, count) that are typically collected from the same subjects across multiple omics platforms. Key aspects:

- **Mixed likelihoods.** Gaussian, Bernoulli, and Binomial features are handled jointly via Pólya-Gamma data augmentation — no pre-processing or separating the data required.
- **Graph-structured shrinkage.** A precision-matrix prior on the factor loadings is informed by a user-supplied feature adjacency graph (e.g. from pathway databases), encouraging features connected in the graph to have correlated loadings and improving interpretability.
- **Missing data.** Missing entries are imputed jointly with the latent factors inside the sampler (data augmentation), so uncertainty about missing values propagates through all estimates automatically.
- **Credible intervals.** 95% posterior credible intervals on every loading identify features with non-zero contribution to each factor, and sparse factors are pruned automatically after burn-in.

---

## Setup

Not on CRAN — source the scripts directly. Install dependencies once:

```r
install.packages(c("MASS", "statmod", "Matrix", "matrixcalc", "corpcor", "BayesLogit"))
```

Then from inside the `iBFA/` directory:

```r
source("R/iBFA.R")
```

---

## Quick example

```r
setwd("path/to/iBFA")
source("examples/example_mixed.R")
```

Or on your own data:

```r
source("R/iBFA.R")

# X   : p x n matrix (rows = features, cols = subjects; NA = missing)
# dt  : integer vector of length p  (0 = Gaussian, 1 = Bernoulli, 2 = Binomial)
# graph: p x p binary adjacency matrix from pathway database (optional)

fit <- iBFA(X, dt, graph = graph, L = 8, T = 2000, burn_in = 1000)

W <- fit$w_est   # p x L  factor loading matrix
Z <- fit$z_est   # L x n  factor score matrix
```

Progress is printed each iteration as `"Task Progress: 0.001"` through `"Task Progress: 1"`.

---

## Building a graph prior

```r
# Two omics blocks, each with 3 pathways of 10 features (G2 = within-pathway edges)
graph <- working_graph(
  x            = 2,
  pathway_list = list(rep(10, 3), rep(10, 3)),
  H            = 2,
  data_dim     = c(30, 30),
  ind_s        = c(1, 31),
  ind_e        = c(30, 60)
)
```

`x = 1`: no graph (zero matrix). `x = 2`: star + random within-pathway edges (G2). `x = 3`: G2 plus sparse across-pathway edges (G3).

---

## `iBFA()` — main function

```r
fit <- iBFA(X, dt, graph = NULL, L = 5, T = 1000, burn_in = 500,
            nu_1 = 1, nu_2 = 0.5, eta = 10, eps = 0.2,
            trials = NULL, missing = NULL)
```

### Inputs

| Argument | Description |
|----------|-------------|
| `X` | p × n data matrix. NA entries are treated as missing. |
| `dt` | Integer vector (length p): `0` = Gaussian, `1` = Bernoulli, `2` = Binomial. |
| `graph` | p × p symmetric binary adjacency matrix. `NULL` = no graph prior (default). |
| `L` | Number of latent factors. Use `iBFA_tune()` if unknown. |
| `T` | Total MCMC iterations. ≥ 2000 for publication; 200–500 for exploration. |
| `burn_in` | Iterations to discard (must be < T). Posterior uses `(burn_in+1):T`. |
| `nu_1` | Prior mean for the log-shrinkage parameter α. Larger = less sparse loadings. Default `1`. |
| `nu_2` | Prior variance for α. Smaller = tighter graph-guided correlation. Default `0.5`. |
| `eta`, `eps` | Shape/rate for the Wishart-like prior on Ω. Defaults `10`, `0.2`. |
| `trials` | p × n integer matrix of binomial trial counts. Required if any `dt[j] == 2`. Defaults to all-ones. |
| `missing` | `TRUE` to impute NAs inside the sampler. Inferred from `is.na(X)` if `NULL`. |

### Outputs (named list)

| Field | Dimensions | Description |
|-------|------------|-------------|
| `w_est` | p × L | Posterior mean factor loading matrix. |
| `z_est` | L × n | Posterior mean latent factor scores. |
| `mu_est` | p × n | Posterior mean linear predictor (WZ + m). Gaussian = data scale; Bernoulli/Binomial = logit scale. |
| `x_est` | vector | Posterior mean imputed values for each NA in X (column-major of mask). Empty if no missing data. |
| `rho_est` | p × n | Posterior mean precision (Gaussian: noise precision; others: PG augmentation). |
| `interval` | 2 × (p·L) | 95% credible interval bounds for each W entry (row-major). |
| `df` | scalar | Effective df: count of W entries with CI excluding zero. |
| `bic_1` | scalar | BIC (all observed + imputed cells). Used by `iBFA_tune()`. |
| `bic_1obs` | scalar | BIC on observed cells only — preferred for selection when data has many missing entries. |
| `keep` | integer vector | Indices of retained factors (after pruning sparse ones). |
| `L`, `L_old` | scalar | Factors after / before pruning. |
| `nu_1`, `nu_2` | scalar | Hyperparameters used. |
| `mask` | p × n logical | TRUE where X was missing. |

### Extracting outputs

```r
L <- fit$L
W <- fit$w_est          # p x L
Z <- fit$z_est          # L x n
mu <- fit$mu_est        # p x n

# Credible-interval matrix for factor l  (p x 2)
cols <- ((l-1)*p + 1):(l*p)
CI_l <- t(fit$interval[, cols])
sig  <- CI_l[, 1] * CI_l[, 2] > 0   # loadings with CI not crossing zero
```

---

## `iBFA_tune()` — hyperparameter selection

```r
source("R/iBFA_tune.R")

tuned <- iBFA_tune(X, dt, graph = graph,
                   grid_L    = c(5, 8, 10),
                   grid_nu_1 = c(-1, 0, 1),
                   grid_nu_2 = c(0.25, 0.5),
                   T = 2000, burn_in = 1000,
                   n_cores = 4)

best <- tuned$best    # full iBFA result for the best grid point
tuned$grid            # data frame: all (L, nu_1, nu_2, bic) combinations
```

Uses `mclapply` on Unix/macOS and `parLapply` on Windows. BIC selection uses `bic_1` by default; switch to `bic_1obs` in `iBFA_tune.R` if data is heavily missing.

---

## R session info

```
R version 4.3.2 (2023-10-31) — Windows 11 x64

MASS 7.3-60 | statmod 1.5.0 | Matrix 1.6-4
matrixcalc 1.0-6 | corpcor 1.6.10 | BayesLogit 2.1
parallel 4.3.2 (base; iBFA_tune only)
```

---

## Files

```
iBFA/
├── R/
│   ├── iBFA.R          # working_graph(), iBFA_mcmc(), iBFA()
│   └── iBFA_tune.R     # BIC grid search: iBFA_tune()
├── examples/
│   └── example_mixed.R # Gaussian + Bernoulli, G2 graph, 20% MAR missingness
└── README.md
```

---

## Citation

> Tsingas K et al. (2026). *Interval Bayesian Factor Analysis for incomplete multi-omics data with graph-structured priors.* [Journal TBD].
