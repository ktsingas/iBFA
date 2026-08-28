## iBFA_tune: BIC grid search over L, nu_1, nu_2
## Author: ktsingas

if (!exists("iBFA")) source(file.path(getwd(), "R", "iBFA.R"))
library(parallel)


# Grid search over L, nu_1, nu_2 by BIC.  Returns the best iBFA fit and the full grid.
#
# grid_nu_1 / grid_nu_2 accept two formats:
#   scalar vector  c(-1, 0, 1)         — each value applied uniformly to all modalities
#   list of vectors list(c(1,0,-1), …) — each element is a per-modality/type vector
#
# modality: integer vector length p (1..H).  When provided, nu_1/nu_2 are indexed by
#   modality block rather than data type.  NULL = index by data type (dt).
#
# n_cores > 1 uses mclapply on Unix/macOS and parLapply on Windows.
iBFA_tune <- function(X, dt, graph = NULL,
                       grid_L    = c(3, 5, 8),
                       grid_nu_1 = c(0, 1),
                       grid_nu_2 = c(0.25, 0.5),
                       eta = 10, eps = 0.2,
                       T = 1000, burn_in = 500,
                       trials = NULL, missing = NULL,
                       modality = NULL,
                       n_cores = 1,
                       print_every = 100) {

  stopifnot(burn_in < T)

  if (!is.list(grid_nu_1)) grid_nu_1 <- as.list(grid_nu_1)
  if (!is.list(grid_nu_2)) grid_nu_2 <- as.list(grid_nu_2)

  idx_grid <- expand.grid(
    L      = grid_L,
    nu1_i  = seq_along(grid_nu_1),
    nu2_i  = seq_along(grid_nu_2),
    stringsAsFactors = FALSE
  )

  nu1_labels <- sapply(grid_nu_1, function(v) paste(round(v, 3), collapse = "/"))
  nu2_labels <- sapply(grid_nu_2, function(v) paste(round(v, 3), collapse = "/"))
  idx_grid$nu_1_label <- nu1_labels[idx_grid$nu1_i]
  idx_grid$nu_2_label <- nu2_labels[idx_grid$nu2_i]

  .run_one <- function(i) {
    iBFA(X        = X,
         dt       = dt,
         graph    = graph,
         L        = idx_grid$L[i],
         T        = T,
         burn_in  = burn_in,
         nu_1     = grid_nu_1[[idx_grid$nu1_i[i]]],
         nu_2     = grid_nu_2[[idx_grid$nu2_i[i]]],
         eta      = eta,
         eps      = eps,
         trials   = trials,
         missing  = missing,
         modality = modality,
         print_every = print_every)
  }

  if (n_cores > 1 && .Platform$OS.type == "unix") {
    results <- mclapply(seq_len(nrow(idx_grid)), .run_one, mc.cores = n_cores)

  } else if (n_cores > 1) {
    cl <- makeCluster(n_cores)
    on.exit(stopCluster(cl), add = TRUE)
    clusterExport(cl,
                  varlist = c("iBFA", "iBFA_mcmc",
                              ".empty_chain", ".p_density_mmh", ".sample_quantile",
                              "X", "dt", "graph", "idx_grid",
                              "grid_nu_1", "grid_nu_2",
                              "T", "burn_in", "eta", "eps",
                              "trials", "missing", "modality", "print_every"),
                  envir = environment())
    clusterEvalQ(cl, {
      library(MASS); library(statmod); library(Matrix)
      library(matrixcalc); library(corpcor); library(BayesLogit)
    })
    results <- parLapply(cl, seq_len(nrow(idx_grid)), .run_one)

  } else {
    results <- lapply(seq_len(nrow(idx_grid)), .run_one)
  }

  ## select by bic_1; switch to bic_1obs when data is heavily missing
  bics     <- sapply(results, function(r) r$bic_1)
  best_idx <- which.min(bics)

  grid <- data.frame(
    L    = idx_grid$L,
    nu_1 = idx_grid$nu_1_label,
    nu_2 = idx_grid$nu_2_label,
    bic  = bics
  )

  cat(sprintf("Best: L=%d  nu_1=%s  nu_2=%s  BIC=%.2f\n",
              grid$L[best_idx], grid$nu_1[best_idx],
              grid$nu_2[best_idx], bics[best_idx]))

  list(best = results[[best_idx]], grid = grid, best_idx = best_idx)
}
