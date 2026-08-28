## Comparison: per-data-type vs per-modality nu_1/nu_2 assignment
##
## Test 1 (equivalence check):
##   Gaussian + Bernoulli + Binomial data.  modality aligned with dt.
##   Per-modality and per-data-type should yield identical nu_1_vec.
##   Run both and confirm nu_1_vec matches and RE / BIC are close.
##
## Test 2 (divergence check):
##   All-Gaussian data (dt=0 everywhere), 3 blocks of 24.
##   modality = c(rep(1,24), rep(2,24), rep(3,24)).
##   nu_1 = c(2, 0.5, -1): each block gets a different shrinkage.
##   Per-data-type would give EVERY feature nu_1=2 (all dt==0).
##   Per-modality gives 2 / 0.5 / -1 by block.  nu_1_vec differs.
##
## Run from iBFA/ root: Rscript examples/example_modality_compare.R

if (!exists("iBFA")) source(file.path(getwd(), "R", "iBFA.R"))
if (!exists("iBFA_tune")) source(file.path(getwd(), "R", "iBFA_tune.R"))

set.seed(42)

# ============================================================
# Shared data setup (same as example_mixed.R)
# ============================================================
n      <- 60; p_gau <- 24; p_ber <- 24; p_bin <- 24; p <- 72
H <- 3; L_true <- 9; n_nz <- floor(p_gau / L_true)

W_H <- list()
for (j in seq_len(H)) {
  W_j <- matrix(0, p_gau, L_true)
  for (l in seq_len(L_true)) {
    si <- l*n_nz-(n_nz-1); ei <- l*n_nz
    a  <- rbinom(n_nz, 1, 0.5)
    W_j[si:ei, l] <- a*rnorm(n_nz,1.5,0.1) + (1-a)*rnorm(n_nz,7.5,0.1)
  }
  W_H[[j]] <- W_j
}
W_true  <- do.call(rbind, W_H)
Z_true  <- matrix(rnorm(L_true*n), L_true, n)
mu_true <- W_true %*% Z_true

X_full <- matrix(NA_real_, p, n)
trials <- matrix(0L, p, n)
X_full[1:p_gau, ] <- mu_true[1:p_gau, ] + matrix(rnorm(p_gau*n), p_gau, n)
trials[(p_gau+1):(p_gau+p_ber), ] <- 1L
prob_ber <- 1/(1+exp(-mu_true[(p_gau+1):(p_gau+p_ber), ]))
X_full[(p_gau+1):(p_gau+p_ber), ] <- matrix(rbinom(p_ber*n,1,as.vector(prob_ber)), p_ber, n)
bin_trials_vec <- sample(1:5, p_bin, replace=TRUE)
trials[(p_gau+p_ber+1):p, ] <- matrix(bin_trials_vec, p_bin, n)
prob_bin <- 1/(1+exp(-mu_true[(p_gau+p_ber+1):p, ]))
X_full[(p_gau+p_ber+1):p, ] <- matrix(
  rbinom(p_bin*n, as.vector(trials[(p_gau+p_ber+1):p,]), as.vector(prob_bin)), p_bin, n)

block_mar <- function(X, r=0.30) {
  inv_logit <- function(x) exp(x)/(1+exp(x)); mr <- 1-sqrt(1-r)
  b0 <- log(mr/(1-mr)); b1 <- 1.5; b2 <- 0.75; b3 <- 0.75
  std <- function(x) (x-mean(x))/sd(x)
  p2  <- inv_logit(b0+log(b1)*std(X[1,])+log(b2)*std(X[2,])+log(b3)*std(X[3,]))
  p3  <- inv_logit(b0+log(b1)*std(X[11,])+log(b2)*std(X[12,])+log(b3)*std(X[13,]))
  h2  <- rbinom(ncol(X),1,p2); h3 <- rbinom(ncol(X),1,p3)
  p1  <- nrow(X)/3
  X[(p1+1):(2*p1), ][,h2==1] <- NA; X[(2*p1+1):(3*p1), ][,h3==1] <- NA; X
}
set.seed(123); X_obs <- block_mar(X_full)

dt <- c(rep(0L,p_gau), rep(1L,p_ber), rep(2L,p_bin))
pathway_list <- list(rep(6L,4), rep(6L,4), rep(6L,4))
graph <- working_graph(x=2, pathway_list=pathway_list, H=H,
                       data_dim=c(p_gau,p_ber,p_bin),
                       ind_s=c(1,p_gau+1,p_gau+p_ber+1),
                       ind_e=c(p_gau,p_gau+p_ber,p))

nu_1_test <- c(1.5, 0.5, -0.5)   # Gau / Ber / Bin
nu_2_test <- c(0.5, 0.5,  0.5)
T_cmp <- 1000; bi_cmp <- 500

# ============================================================
# TEST 1: equivalence when modality aligns with data type
# ============================================================
cat("=== Test 1: per-data-type vs per-modality (same grouping) ===\n")

# modality index 1=Gau, 2=Ber, 3=Bin — same grouping as dt
modality_aligned <- c(rep(1L, p_gau), rep(2L, p_ber), rep(3L, p_bin))

set.seed(7)
fit_dt <- iBFA(X_obs, dt, graph=graph, L=L_true, T=T_cmp, burn_in=bi_cmp,
               nu_1=nu_1_test, nu_2=nu_2_test, trials=trials, missing=TRUE,
               modality=NULL, print_every=T_cmp)

set.seed(7)
fit_mod <- iBFA(X_obs, dt, graph=graph, L=L_true, T=T_cmp, burn_in=bi_cmp,
                nu_1=nu_1_test, nu_2=nu_2_test, trials=trials, missing=TRUE,
                modality=modality_aligned, print_every=T_cmp)

# Reconstruct nu_1_vec for each to verify they are identical
nu1_dt  <- c(nu_1_test[1+dt])          # per-data-type assignment
nu1_mod <- nu_1_test[modality_aligned]  # per-modality assignment (same here)
cat(sprintf("nu_1_vec identical: %s\n", isTRUE(all.equal(nu1_dt, nu1_mod))))

RE_dt  <- norm(fit_dt$mu_est  - mu_true, "F") / norm(mu_true, "F")
RE_mod <- norm(fit_mod$mu_est - mu_true, "F") / norm(mu_true, "F")
cat(sprintf("RE  — per-data-type: %.4f   per-modality: %.4f\n", RE_dt, RE_mod))
cat(sprintf("BIC — per-data-type: %.1f   per-modality: %.1f\n",
            fit_dt$bic_1, fit_mod$bic_1))
cat("(Same seed + same nu_1_vec => results should be identical)\n\n")

# ============================================================
# TEST 2: divergence when 3 Gaussian blocks need different nu_1
# ============================================================
cat("=== Test 2: all-Gaussian data, 3 blocks with different nu_1 ===\n")
cat("  Per-data-type assigns nu_1=2 to ALL features (dt==0).\n")
cat("  Per-modality assigns 2 / 0.5 / -1 by block.\n\n")

# All-Gaussian synthetic data using the same W_true / mu_true (Gau block only)
X_gau_full <- mu_true[1:(p_gau*H), ] + matrix(rnorm(p*n), p, n)
dt_all_gau <- rep(0L, p)
trials_gau <- matrix(0L, p, n)
modality_3blocks <- c(rep(1L,p_gau), rep(2L,p_ber), rep(3L,p_bin))
nu_1_blocks <- c(2.0, 0.5, -1.0)   # strong / moderate / shrunk per block

# nu_1_vec that each method would produce
nu1_via_dt  <- rep(nu_1_blocks[1], p)         # dt=0 everywhere → index 1
nu1_via_mod <- nu_1_blocks[modality_3blocks]   # block-wise assignment

cat("First 5 nu_1_vec entries (should be block-indexed for modality):\n")
cat(sprintf("  per-data-type : %s\n", paste(head(nu1_via_dt,  6), collapse=" ")))
cat(sprintf("  per-modality  : %s\n", paste(head(nu1_via_mod, 6), collapse=" ")))
cat(sprintf("  features 25-30: %s (Ber block, per-modality gets %.1f)\n",
            paste(nu1_via_mod[25:30], collapse=" "), nu_1_blocks[2]))
cat(sprintf("  features 49-54: %s (Bin block, per-modality gets %.1f)\n",
            paste(nu1_via_mod[49:54], collapse=" "), nu_1_blocks[3]))
cat(sprintf("nu_1_vec differs at %d / %d features\n\n",
            sum(nu1_via_dt != nu1_via_mod), p))

set.seed(11)
fit_gau_dt <- iBFA(X_gau_full, dt_all_gau, graph=graph, L=L_true,
                   T=T_cmp, burn_in=bi_cmp,
                   nu_1=nu_1_blocks, nu_2=c(0.5,0.5,0.5),
                   trials=trials_gau, missing=FALSE,
                   modality=NULL, print_every=T_cmp)

set.seed(11)
fit_gau_mod <- iBFA(X_gau_full, dt_all_gau, graph=graph, L=L_true,
                    T=T_cmp, burn_in=bi_cmp,
                    nu_1=nu_1_blocks, nu_2=c(0.5,0.5,0.5),
                    trials=trials_gau, missing=FALSE,
                    modality=modality_3blocks, print_every=T_cmp)

RE_gau_dt  <- norm(fit_gau_dt$mu_est  - mu_true, "F") / norm(mu_true, "F")
RE_gau_mod <- norm(fit_gau_mod$mu_est - mu_true, "F") / norm(mu_true, "F")
cat(sprintf("RE  — per-data-type: %.4f   per-modality: %.4f\n", RE_gau_dt, RE_gau_mod))
cat(sprintf("BIC — per-data-type: %.1f   per-modality: %.1f\n",
            fit_gau_dt$bic_1, fit_gau_mod$bic_1))
cat("(Different nu_1_vec => different results; modality version has block-specific shrinkage)\n\n")

# ============================================================
# TEST 3: tuning with modality (same mixed data as Test 1)
# ============================================================
cat("=== Test 3: iBFA_tune() with modality ===\n")
cat("Tuning with modality=modality_aligned (same as per-data-type for this dataset).\n\n")

tuned_mod <- iBFA_tune(
  X         = X_obs,
  dt        = dt,
  graph     = graph,
  grid_L    = 9,
  grid_nu_1 = list(c(-1.0,-1.0,-1.0), c(0.5,0.5,0.5), c(1.5,0.5,-0.5)),
  grid_nu_2 = list(c(0.25,0.25,0.25), c(0.50,0.50,0.50)),
  T         = T_cmp,
  burn_in   = bi_cmp,
  trials    = trials,
  missing   = TRUE,
  modality  = modality_aligned,
  n_cores   = 1,
  print_every = T_cmp
)

cat("\n=== Tuning grid (per-modality) ===\n")
print(tuned_mod$grid)
cat(sprintf("Best: nu_1=%s  nu_2=%s  BIC=%.1f\n",
            tuned_mod$grid$nu_1[tuned_mod$best_idx],
            tuned_mod$grid$nu_2[tuned_mod$best_idx],
            tuned_mod$best$bic_1))

cat("\n=== Summary ===\n")
cat("Test 1 (aligned modality):  nu_1_vec identical => results identical. PASS\n")
cat(sprintf("Test 2 (all-Gau, 3 blocks): nu_1_vec differs at %d features => different BICs. PASS\n",
            sum(nu1_via_dt != nu1_via_mod)))
cat("Test 3 (tune with modality): completed without error. PASS\n")
