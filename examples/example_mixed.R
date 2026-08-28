## iBFA example: three Gaussian modality blocks with per-modality shrinkage tuning
## Shows: block MAR missingness | G2 graph prior | modality-indexed nu_1/nu_2 | grid search
##
## Three omics platforms (e.g. RNA, protein, metabolite) all measured continuously.
## Because all features share dt=0, modality is required to assign different
## nu_1/nu_2 per platform — per-data-type alone cannot distinguish the three blocks.
##
## [Reference] Mixed Gau/Ber/Bin variant is commented out below.
## Run from iBFA/ root: Rscript examples/example_mixed.R

if (!exists("iBFA"))      source(file.path(getwd(), "R", "iBFA.R"))
if (!exists("iBFA_tune")) source(file.path(getwd(), "R", "iBFA_tune.R"))

set.seed(42)

# ============================================================
# Hyperparameters — set everything here
# ============================================================

n      <- 60            # subjects
p_per  <- 24            # features per modality block
H      <- 3             # number of modality blocks (e.g. 3 omics platforms)
p      <- H * p_per     # 72 total features
L_true <- 9             # true number of latent factors
n_nz   <- floor(p_per / L_true)   # 2: non-zero loadings per factor per block

T_iter  <- 3000         # total MCMC iterations
burn_in <- T_iter - 1000  # discard first 2000; posterior uses last 1000 samples

eta  <- 10              # Wishart prior strength on Omega (larger = stronger graph pull)
eps  <- 0.2             # diagonal offset in Wishart prior (keeps Omega positive-definite)

# nu_1: prior mean on log-shrinkage alpha — larger allows bigger loadings
# nu_2: prior variance on alpha — smaller tightens graph-guided feature coupling
# One value per modality block (block 1 / block 2 / block 3)
nu_1_init <- c( 2.0,  0.5, -1.0)
nu_2_init <- c( 0.5,  0.5,  0.5)

# Tuning grid — each list element is a c(block1, block2, block3) vector
# grid_L: include true L; widen range if L is unknown
# nu_1 options must not include 0
grid_L    <- c(7, 9, 11)
grid_nu_1 <- list(c(-1.0, -1.0, -1.0),    # uniform heavy shrinkage
                  c( 0.5,  0.5,  0.5),     # uniform moderate
                  c( 1.5,  0.5, -1.0))     # block-specific (relaxed / moderate / shrunk)
grid_nu_2 <- list(c(0.25, 0.25, 0.25),
                  c(0.50, 0.50, 0.50))
T_tune  <- 3000
bi_tune <- T_tune - 1000

# ============================================================
# Data generation — all Gaussian, 3 modality blocks
# Matches generate_X(dt=0, g=2) in sim_batch_same.R (apart from n, p).
# W_true: bimodal loadings a*N(1.5,0.1) + (1-a)*N(7.5,0.1); factor l loads on
#   features (l-1)*n_nz+1 : l*n_nz per block. Features (L*n_nz+1):p_per are noise.
# ============================================================

W_H <- list()
for (j in seq_len(H)) {
  W_j <- matrix(0, p_per, L_true)
  for (l in seq_len(L_true)) {
    si <- l*n_nz-(n_nz-1); ei <- l*n_nz
    a  <- rbinom(n_nz, 1, 0.5)
    W_j[si:ei, l] <- a*rnorm(n_nz, 1.5, 0.1) + (1-a)*rnorm(n_nz, 7.5, 0.1)
  }
  W_H[[j]] <- W_j
}
W_true  <- do.call(rbind, W_H)
Z_true  <- matrix(rnorm(L_true * n), L_true, n)
mu_true <- W_true %*% Z_true

X_full  <- mu_true + matrix(rnorm(p * n), p, n)

dt       <- rep(0L, p)                                      # all Gaussian
trials   <- matrix(0L, p, n)                                # unused for Gaussian
modality <- c(rep(1L, p_per), rep(2L, p_per), rep(3L, p_per))  # block index per feature

cat(sprintf("Data: p=%d (%d features x %d Gaussian blocks), n=%d, L_true=%d, n_nz=%d\n",
            p, p_per, H, n, L_true, n_nz))

# ============================================================
# [Reference] Mixed Gau / Ber / Bin variant
# Three different data types: remove modality arg and set dt accordingly.
# ============================================================
# p_gau <- 24; p_ber <- 24; p_bin <- 24
# trials <- matrix(0L, p, n)
# trials[(p_gau+1):(p_gau+p_ber), ] <- 1L
# bin_tv <- sample(1:5, p_bin, replace=TRUE)
# trials[(p_gau+p_ber+1):p, ] <- matrix(bin_tv, p_bin, n)
# prob_ber <- 1/(1+exp(-mu_true[(p_gau+1):(p_gau+p_ber),]))
# prob_bin <- 1/(1+exp(-mu_true[(p_gau+p_ber+1):p,]))
# X_full[(p_gau+1):(p_gau+p_ber),] <- matrix(rbinom(p_ber*n,1,as.vector(prob_ber)),p_ber,n)
# X_full[(p_gau+p_ber+1):p,] <- matrix(
#   rbinom(p_bin*n, as.vector(trials[(p_gau+p_ber+1):p,]), as.vector(prob_bin)), p_bin, n)
# dt <- c(rep(0L,p_gau), rep(1L,p_ber), rep(2L,p_bin))
# modality <- NULL   # dt already distinguishes the three blocks; no need for modality

# ============================================================
# Block MAR missing data
# Replicates block(mech="MAR") from sim_batch_same.R.
# Blocks 2 and 3 missingness is driven by Gaussian features 1-3 and 11-13 of block 1.
# Block 1 is always fully observed (used as predictor).
# ============================================================

block_mar <- function(X, r = 0.30) {
  inv_logit <- function(x) exp(x) / (1 + exp(x))
  mr  <- 1 - sqrt(1 - r)
  b0  <- log(mr / (1 - mr)); b1 <- 1.5; b2 <- 0.75; b3 <- 0.75
  std <- function(x) (x - mean(x)) / sd(x)
  p2  <- inv_logit(b0 + log(b1)*std(X[1,])  + log(b2)*std(X[2,])  + log(b3)*std(X[3,]))
  p3  <- inv_logit(b0 + log(b1)*std(X[11,]) + log(b2)*std(X[12,]) + log(b3)*std(X[13,]))
  h2  <- rbinom(ncol(X), 1, p2); h3 <- rbinom(ncol(X), 1, p3)
  p1  <- nrow(X) / 3
  X[(p1+1):(2*p1),   ][, h2 == 1] <- NA
  X[(2*p1+1):(3*p1), ][, h3 == 1] <- NA
  X
}

set.seed(123)
X_obs <- block_mar(X_full, r = 0.30)
cat(sprintf("Missing: %.1f%% overall (block MAR on blocks 2 and 3)\n", 100*mean(is.na(X_obs))))

# ============================================================
# G2 graph: 4 pathways x 6 features per block
# ============================================================

graph <- working_graph(
  x            = 2,
  pathway_list = list(rep(6L, 4), rep(6L, 4), rep(6L, 4)),
  H            = H,
  data_dim     = rep(p_per, H),
  ind_s        = c(1, p_per+1, 2*p_per+1),
  ind_e        = c(p_per, 2*p_per, p)
)
cat(sprintf("Graph edges: %d (G2 within-pathway, 4 pathways x 6 features per block)\n",
            sum(graph) / 2))

prog_file <- file.path(getwd(), "iBFA_progress.txt")

# ============================================================
# Primary run — initial per-modality hyperparameters
# ============================================================

cat(sprintf("\nPrimary fit: T=%d, burn_in=%d ...\n", T_iter, burn_in))

fit_primary <- iBFA(
  X             = X_obs,
  dt            = dt,
  graph         = graph,
  L             = L_true,
  T             = T_iter,
  burn_in       = burn_in,
  nu_1          = nu_1_init,
  nu_2          = nu_2_init,
  eta           = eta,
  eps           = eps,
  trials        = trials,
  missing       = TRUE,
  modality      = modality,
  print_every   = 300,
  progress_file = prog_file
)

RE_primary <- norm(fit_primary$mu_est - mu_true, "F") / norm(mu_true, "F")
cat(sprintf("Primary — L=%d  BIC=%.1f  BIC(obs)=%.1f  RE=%.4f\n",
            fit_primary$L, fit_primary$bic_1, fit_primary$bic_1obs, RE_primary))

# ============================================================
# Tuning — grid over L, nu_1, nu_2 with modality-indexed hyperparameters
# grid_nu_1 / grid_nu_2: each list element is a c(block1, block2, block3) vector
# ============================================================

cat(sprintf("\nTuning: %d grid points (L x nu_1 x nu_2), T=%d each ...\n",
            length(grid_L) * length(grid_nu_1) * length(grid_nu_2), T_tune))

tuned <- iBFA_tune(
  X           = X_obs,
  dt          = dt,
  graph       = graph,
  grid_L      = grid_L,
  grid_nu_1   = grid_nu_1,
  grid_nu_2   = grid_nu_2,
  eta         = eta,
  eps         = eps,
  T           = T_tune,
  burn_in     = bi_tune,
  trials      = trials,
  missing     = TRUE,
  modality    = modality,
  n_cores     = 1,           # increase on Unix/macOS for parallel grid search
  print_every = 300
)

cat("\n===  Tuning grid  ===\n")
print(tuned$grid)

RE_tuned <- norm(tuned$best$mu_est - mu_true, "F") / norm(mu_true, "F")

cat(sprintf("\n===  Results  ===\n"))
cat(sprintf("Primary   L=%d  nu_1=%s  nu_2=%s  BIC=%.1f  RE=%.4f\n",
            fit_primary$L,
            paste(nu_1_init, collapse="/"), paste(nu_2_init, collapse="/"),
            fit_primary$bic_1, RE_primary))
cat(sprintf("Best tune L=%d  nu_1=%s  nu_2=%s  BIC=%.1f  RE=%.4f\n",
            tuned$best$L,
            tuned$grid$nu_1[tuned$best_idx],
            tuned$grid$nu_2[tuned$best_idx],
            tuned$best$bic_1, RE_tuned))

write(c("",
        "--- Tuning results ---",
        capture.output(print(tuned$grid)),
        sprintf("Best: nu_1=%s  nu_2=%s  L=%d  BIC=%.2f  RE=%.4f",
                tuned$grid$nu_1[tuned$best_idx],
                tuned$grid$nu_2[tuned$best_idx],
                tuned$best$L, tuned$best$bic_1, RE_tuned),
        sprintf("Primary RE: %.4f", RE_primary)),
      file = prog_file, append = TRUE)

cat(sprintf("\nProgress file: %s\n", prog_file))
