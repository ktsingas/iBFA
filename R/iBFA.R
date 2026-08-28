## iBFA: interval Bayesian Factor Analysis for mixed-type incomplete multi-omics data
## Author: ktsingas
## Deps: MASS, statmod, Matrix, matrixcalc, corpcor, BayesLogit

library(MASS)
library(statmod)
library(Matrix)
library(matrixcalc)
library(corpcor)
library(BayesLogit)


## x=1: no graph; x=2: star + within-pathway random edges (G2); x=3: G2 + across-pathway edges (G3)
working_graph <- function(x, pathway_list, H, data_dim, ind_s, ind_e) {
  p     <- sum(data_dim)
  graph <- matrix(0, p, p)
  if (x == 1) return(graph)

  for (h in 1:H) {
    graph_tmp   <- matrix(0, data_dim[h], data_dim[h])
    pathway_tmp <- pathway_list[[h]]
    node <- c(1)
    for (xx in seq_len(length(pathway_tmp) - 1))
      node <- c(node, 1 + sum(pathway_tmp[1:xx]))

    for (nn in seq_along(node)) {
      r_ind <- node[nn]
      end_n <- sum(pathway_tmp[1:nn])

      if (end_n > r_ind)
        graph_tmp[r_ind, seq(r_ind + 1, end_n)] <- 1

      if (x == 2 && end_n - r_ind > 1) {
        inner <- (r_ind + 1):end_n
        graph_tmp[inner, inner] <- rbinom(length(inner)^2, 1, 0.25)
      }

      if (x == 3) {
        graph_tmp[(r_ind + 1):end_n, (r_ind + 1):end_n] <-
          rbinom((pathway_tmp[nn] - 1)^2, 1, 0)
        if (end_n + 1 <= sum(pathway_tmp)) {
          rest <- (end_n + 1):sum(pathway_tmp)
          graph_tmp[r_ind:end_n, rest] <-
            rbinom(length(r_ind:end_n) * length(rest), 1, 0.1)
        }
      }
    }
    graph[ind_s[h]:ind_e[h], ind_s[h]:ind_e[h]] <- graph_tmp
  }

  graph <- as.matrix(forceSymmetric(graph, uplo = "U"))
  diag(graph) <- 0
  graph
}


## internal helpers ----

.empty_chain <- function(n, p, T, L, rho_ini, w_ini, tau_ini, z_ini, alpha_ini, mask) {
  ch <- list(
    chain_m     = matrix(0, T + 2, p),
    chain_rho   = matrix(0, T + 2, p * n),
    chain_w     = matrix(0, T + 2, p * L),
    chain_tau   = matrix(0, T + 2, p * L),
    chain_z     = matrix(0, T + 2, L * n),
    chain_alpha = matrix(0, T + 2, p * L),
    chain_x     = matrix(0, T + 2, sum(mask))
  )
  ch$chain_rho[1, ]   <- as.vector(t(rho_ini))
  ch$chain_w[1, ]     <- as.vector(t(w_ini))
  ch$chain_tau[1, ]   <- as.vector(t(tau_ini))
  ch$chain_z[1, ]     <- as.vector(t(z_ini))
  ch$chain_alpha[1, ] <- as.vector(t(alpha_ini))
  ch
}

## log-ratio for alpha Metropolis-Hastings accept step
.p_density_mmh <- function(alpha_l, alpha_h, tau_l, alpha, omega, j, nu_1j, nu_2j) {
  ll <- exp(alpha_l); lh <- exp(alpha_h)
  alp_aug <- replace(alpha, j, alpha_l)
  (ll^2 / lh^2) *
    exp(-(ll^2 - lh^2) * tau_l / 2) *
    exp((-1 / (2 * nu_2j)) * (alpha_l - alpha_h) *
          (omega[j, ] %*% (alpha + alp_aug - 2 * nu_1j)))
}

.sample_quantile <- function(x) quantile(x, probs = c(0.025, 0.975))


## MCMC sampler ----

iBFA_mcmc <- function(dt, T, L, p, n, X, trials,
                       nu_1, nu_2, Sigma, Q, eta, eps,
                       rho_temp, tau_temp, omega_temp, inv_omega_temp,
                       w_temp, z_temp, alpha_temp, m_temp,
                       start, end,
                       chain_rho, chain_alpha, chain_tau,
                       chain_w, chain_z, chain_m, chain_x,
                       missing = FALSE,
                       modality = NULL,
                       print_every = 100,
                       progress_file = NULL) {

  no_dt  <- unique(dt)
  dt_gau <- which(dt == 0)
  dt_ber <- which(dt == 1)
  dt_bin <- which(dt == 2)
  dt_nb  <- which(dt == 3)

  ## build per-feature nu vectors
  ## modality (integer 1..H per feature) takes priority over per-data-type grouping
  if (!is.null(modality)) {
    nu_1_vec <- numeric(p); nu_2_vec <- numeric(p)
    for (m in seq_along(nu_1)) {
      idx <- which(modality == m)
      if (length(idx) > 0) { nu_1_vec[idx] <- nu_1[m]; nu_2_vec[idx] <- nu_2[m] }
    }
  } else if (length(nu_1) > 1) {
    nu_1_vec <- numeric(p); nu_2_vec <- numeric(p)
    mod_list <- list(dt_gau, dt_ber, dt_bin, dt_nb)
    for (m in seq_along(nu_1)) {
      idx <- mod_list[[m]]
      if (length(idx) > 0) { nu_1_vec[idx] <- nu_1[m]; nu_2_vec[idx] <- nu_2[m] }
    }
  } else {
    nu_1_vec <- rep(nu_1, p)
    nu_2_vec <- rep(nu_2, p)
  }

  ## initial imputation and sufficient statistics
  if (missing) {
    mask <- is.na(X)
    mask_gau <- if (length(dt_gau) > 0) is.na(X[dt_gau, ]) else matrix(FALSE, 0, n)
    mask_ber <- if (length(dt_ber) > 0) is.na(X[dt_ber, ]) else matrix(FALSE, 0, n)
    mask_bin <- if (length(dt_bin) > 0) is.na(X[dt_bin, ]) else matrix(FALSE, 0, n)

    if (length(dt_gau) > 0)
      X[dt_gau, ] <- t(apply(X[dt_gau, , drop = FALSE], 1, function(row) {
        row[is.na(row)] <- mean(row, na.rm = TRUE); row }))

    if (length(dt_ber) > 0)
      X[dt_ber, ] <- t(apply(X[dt_ber, , drop = FALSE], 1, function(row) {
        mv <- max(row, na.rm = TRUE); sp <- sum(row, na.rm = TRUE) / (mv * sum(!is.na(row)))
        row[is.na(row)] <- rbinom(sum(is.na(row)), mv, sp); row }))

    if (length(dt_bin) > 0)
      X[dt_bin, ] <- t(apply(X[dt_bin, , drop = FALSE], 1, function(row) {
        row[is.na(row)] <- median(row, na.rm = TRUE); row }))

    chain_x[1, ] <- as.vector(t(X[mask]))
  } else {
    mask     <- matrix(FALSE, p, n)
    mask_gau <- matrix(FALSE, length(dt_gau), n)
    mask_ber <- matrix(FALSE, length(dt_ber), n)
    mask_bin <- matrix(FALSE, length(dt_bin), n)
  }

  psi             <- matrix(0, p, n)
  if (length(dt_gau) > 0) psi[dt_gau, ] <- X[dt_gau, ]
  phi             <- as.matrix(rep(1, p))
  b               <- matrix(1, p, n)
  if (length(dt_bin) > 0) b[dt_bin, ] <- trials[dt_bin, ]
  if (length(dt_nb)  > 0) b[dt_nb,  ] <- X[dt_nb, ] + trials[dt_nb, ]
  kappa           <- X - b / 2
  if (length(dt_gau) > 0) kappa[dt_gau, ] <- 0
  rowsum_kappa    <- rowSums(kappa)

  ## MCMC loop
  for (t in 1:T) {

    rowsum_rho  <- rowSums(rho_temp)
    row_rho_psi <- rowSums(rho_temp * psi)

    ## m: intercept, conjugate Gaussian update
    covar_m <- 1 / (rowsum_rho + 1 / Sigma)
    mean_m  <- covar_m * (rowsum_kappa + row_rho_psi -
                            rowSums((w_temp %*% z_temp) * rho_temp))
    m_temp  <- rnorm(p, mean_m, sqrt(covar_m))
    chain_m[t + 1, ] <- as.vector(m_temp)

    ## rho: Gaussian precision via Gamma; non-Gaussian via Polya-Gamma augmentation
    mu_all <- matrix(m_temp, p, n) + w_temp %*% z_temp
    if (length(dt_gau) > 0) {
      r <- phi[dt_gau, 1] + rowSums((X[dt_gau, ] - mu_all[dt_gau, ])^2)
      rho_temp[dt_gau, ] <- matrix(
        rgamma(length(dt_gau), shape = (phi[dt_gau, 1] + n) / 2, rate = r / 2),
        nrow = length(dt_gau), ncol = n)
    }
    for (j in c(dt_ber, dt_bin))
      rho_temp[j, ] <- rpg.devroye(n, b[j, ], mu_all[j, ])
    chain_rho[t + 1, ] <- as.vector(t(rho_temp))

    ## alpha: log-shrinkage, Metropolis-Hastings per feature
    for (l in 1:L) for (j in 1:p) {
      candi <- rnorm(1, alpha_temp[j, l], 1)
      comp  <- as.numeric(.p_density_mmh(candi, alpha_temp[j, l],
                                          tau_temp[j, l], alpha_temp[, l],
                                          omega_temp, j, nu_1_vec[j], nu_2_vec[j]))
      if (is.na(comp)) comp <- 0
      if (runif(1) < min(comp, 1)) alpha_temp[j, l] <- candi
    }
    chain_alpha[t + 1, ] <- as.vector(t(alpha_temp))

    ## Omega: graph-structured precision on alpha columns
    A_resid <- (alpha_temp - nu_1_vec) / sqrt(nu_2_vec)
    A <- eta * (diag(eps, p) + 1) + A_resid %*% t(A_resid)

    for (j in 1:p) {
      A_jj    <- A[j, j]
      inv_o12 <- inv_omega_temp[-j, j]
      inv_o22 <- inv_omega_temp[j, j]
      O_11    <- inv_omega_temp[-j, -j] - outer(inv_o12, inv_o12) / inv_o22
      nz_ind  <- setdiff(which(omega_temp[, j] != 0), j)

      if (length(nz_ind) > 1) {
        prec_nz        <- inv_omega_temp[nz_ind, nz_ind] -
          outer(inv_omega_temp[nz_ind, j], inv_omega_temp[nz_ind, j]) / inv_o22
        prec_nz_sqrt   <- chol(prec_nz)
        mean_nz        <- backsolve(prec_nz_sqrt,
                                    forwardsolve(t(prec_nz_sqrt), (-1/A_jj)*A[nz_ind,j]))
        non_zero_omega <- backsolve(prec_nz_sqrt, rnorm(length(nz_ind), 0, 1/sqrt(A_jj))) + mean_nz
        xi             <- rgamma(1, shape = 1 + (eta*(1+eps)+L)/2, rate = A_jj/2)
        diag_omega     <- xi + sum((prec_nz_sqrt %*% non_zero_omega)^2)

      } else if (length(nz_ind) == 1) {
        prec_nz        <- inv_omega_temp[nz_ind, nz_ind] -
          inv_omega_temp[nz_ind, j]^2 / inv_o22
        non_zero_omega <- rnorm(1, (-1/(A_jj*prec_nz))*A[nz_ind,j], 1/sqrt(prec_nz*A_jj))
        xi             <- rgamma(1, shape = 1 + (eta*(1+eps)+L)/2, rate = A_jj/2)
        diag_omega     <- xi + non_zero_omega^2 * prec_nz

      } else {
        xi             <- rgamma(1, shape = 1 + (eta*(1+eps)+L)/2, rate = A_jj/2)
        diag_omega     <- xi
        non_zero_omega <- numeric(0)
      }

      if (length(nz_ind) > 0) {
        omega_temp[nz_ind, j] <- non_zero_omega
        omega_temp[j, nz_ind] <- non_zero_omega
      }
      omega_temp[j, j] <- diag_omega

      inv_omega_temp[j, j]   <- 1 / xi
      pp                     <- O_11 %*% omega_temp[-j, j]
      inv_omega_11           <- O_11 + (1/xi) * (pp %*% t(pp))
      inv_omega_12           <- -inv_omega_11 %*% omega_temp[-j, j] / diag_omega
      inv_omega_temp[-j, -j] <- inv_omega_11
      inv_omega_temp[-j, j]  <- inv_omega_12
      inv_omega_temp[j, -j]  <- inv_omega_12
    }

    ## Z: latent factor scores
    for (i in 1:n) {
      pt_cov      <- t(w_temp) %*% (w_temp * rho_temp[, i])
      diag(pt_cov) <- diag(pt_cov) + 1
      chol_cov    <- chol(pt_cov)
      pt_mean     <- t(w_temp) %*% (rho_temp[, i] * (psi[, i] - m_temp) + kappa[, i])
      z_temp[, i] <- backsolve(chol_cov, rnorm(L)) +
        backsolve(chol_cov, forwardsolve(t(chol_cov), pt_mean))
    }
    chain_z[t + 1, ] <- as.vector(t(z_temp))

    ## W: loadings; tau: local shrinkage (inverse Gaussian)
    for (j in 1:p) {
      cov_pt      <- z_temp %*% (t(z_temp) * rho_temp[j, ])
      diag(cov_pt) <- diag(cov_pt) + 1 / tau_temp[j, ]
      chol_cov    <- chol(cov_pt)
      pt_mean     <- t(t(z_temp) * rho_temp[j, ]) %*%
        (psi[j, ] - m_temp[j] + kappa[j, ] / rho_temp[j, ])
      w_temp[j, ] <- backsolve(chol_cov, rnorm(L)) +
        backsolve(chol_cov, forwardsolve(t(chol_cov), pt_mean))
      tau_temp[j, ] <- 1 / rinvgauss(L,
                                      abs(exp(alpha_temp[j, ]) / w_temp[j, ]),
                                      exp(2 * alpha_temp[j, ]))
    }
    chain_w[t + 1, ]   <- as.vector(t(w_temp))
    chain_tau[t + 1, ] <- as.vector(t(tau_temp))

    if (t %% print_every == 0 || t == T)
      print(paste0("Task Progress: ", round(t / T, 4)))

    if (!is.null(progress_file) && (t == 1 || t %% 50 == 0 || t == T)) {
      writeLines(c(
        sprintf("iBFA progress — %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
        sprintf("Iteration : %d / %d  (%.1f%%)", t, T, 100 * t / T),
        sprintf("Burn-in   : %d  |  Post-burn-in starts at %d", start - 1, start),
        if (t < T) "Status    : running" else "Status    : DONE"
      ), progress_file)
    }

    ## imputation: resample missing entries from their conditional posterior
    if (missing) {
      mu_temp <- m_temp + w_temp %*% z_temp
      if (length(dt_gau) > 0 && any(mask_gau)) {
        X[dt_gau, ][mask_gau]   <- rnorm(sum(mask_gau),
                                          mu_temp[dt_gau, ][mask_gau],
                                          sqrt(1 / rho_temp[dt_gau, ][mask_gau]))
        psi[dt_gau, ][mask_gau] <- X[dt_gau, ][mask_gau]
      }
      if (length(dt_ber) > 0 && any(mask_ber)) {
        X[dt_ber, ][mask_ber]     <- rbinom(sum(mask_ber), 1,
                                             1 / (1 + exp(-mu_temp[dt_ber, ][mask_ber])))
        kappa[dt_ber, ][mask_ber] <- X[dt_ber, ][mask_ber] - b[dt_ber, ][mask_ber] / 2
      }
      if (length(dt_bin) > 0 && any(mask_bin)) {
        X[dt_bin, ][mask_bin]     <- rbinom(sum(mask_bin),
                                             trials[dt_bin, ][mask_bin],
                                             1 / (1 + exp(-mu_temp[dt_bin, ][mask_bin])))
        kappa[dt_bin, ][mask_bin] <- X[dt_bin, ][mask_bin] - b[dt_bin, ][mask_bin] / 2
      }
      chain_x[t + 1, ] <- as.vector(t(X[mask]))
    }
  }

  ## post-burn summarise: prune factors with <2 credible non-zero loadings
  .summarise <- function(s, e, L_cur) {
    mu_e <- matrix(0, p, n)
    for (i in s:e) {
      w_t <- matrix(as.numeric(chain_w[i, ]), p, L_cur, byrow = TRUE)
      z_t <- matrix(as.numeric(chain_z[i, ]), L_cur, n, byrow = TRUE)
      mu_e <- mu_e + w_t %*% z_t + as.numeric(chain_m[i, ])
    }
    mu_e <- mu_e / (e - s + 1)

    w_e   <- apply(chain_w[s:e, ],   2, mean)
    z_e   <- apply(chain_z[s:e, ],   2, mean)
    m_e   <- apply(chain_m[s:e, ],   2, mean)
    rho_e <- matrix(apply(chain_rho[s:e, ], 2, mean), p, n, byrow = TRUE)
    x_e   <- apply(chain_x[s:e, ],   2, mean)
    intv  <- apply(chain_w[s:e, ], 2, .sample_quantile)
    df_   <- sum(intv[1, ] * intv[2, ] > 0)

    cred_mat <- matrix((intv[1, ] * intv[2, ]) > 0, p, L_cur, byrow = TRUE)
    keep_idx <- which(colSums(cred_mat) >= 2)
    if (length(keep_idx) == 0) keep_idx <- seq_len(L_cur)
    L_old <- L_cur; L_new <- length(keep_idx)

    if (L_new < L_cur) {
      idx_w <- as.vector(sapply(seq_len(p), function(j) (j-1L)*L_cur + keep_idx))
      idx_z <- unlist(lapply(keep_idx, function(l) ((l-1L)*n+1L):(l*n)))
      mu_e  <- matrix(0, p, n)
      for (i in s:e) {
        wk   <- matrix(as.numeric(chain_w[i,]), p, L_cur, byrow=TRUE)[, keep_idx, drop=FALSE]
        zk   <- matrix(as.numeric(chain_z[i,]), L_cur, n, byrow=TRUE)[keep_idx, , drop=FALSE]
        mu_e <- mu_e + wk %*% zk + as.numeric(chain_m[i,])
      }
      mu_e  <- mu_e / (e - s + 1L)
      cw    <- chain_w[s:e, idx_w, drop = FALSE]
      cz    <- chain_z[s:e, idx_z, drop = FALSE]
      w_e   <- apply(cw, 2, mean)
      z_e   <- apply(cz, 2, mean)
      intv  <- apply(cw, 2, .sample_quantile)
      df_   <- sum((intv[1,] * intv[2,]) > 0)
      L_new_ <- L_new
    } else {
      keep_idx <- seq_len(L_cur)
      L_new_   <- L_cur
    }

    list(w_est = w_e, z_est = z_e, mu_est = mu_e, x_est = x_e,
         rho_est = rho_e, m_est = m_e, interval = intv, df = df_,
         keep = keep_idx, L = L_new_, L_old = L_old)
  }

  res <- .summarise(start, end, L)

  ## BIC: bic_1 on all cells; bic_1obs on observed cells only (prefer when data heavily missing)
  X_cur <- X; mu_e <- res$mu_est; rho_e <- res$rho_est
  ll_gau <- if (length(dt_gau) > 0)
    0.5*log(rho_e[dt_gau,]) - 0.5*log(2*pi) -
    0.5*rho_e[dt_gau,]*(X_cur[dt_gau,] - mu_e[dt_gau,])^2 else numeric(0)
  ll_ber <- if (length(dt_ber) > 0)
    log(choose(trials[dt_ber,], X_cur[dt_ber,])) +
    mu_e[dt_ber,]*X_cur[dt_ber,] - log(1+exp(mu_e[dt_ber,]))*trials[dt_ber,] else numeric(0)
  ll_bin <- if (length(dt_bin) > 0)
    log(choose(trials[dt_bin,], X_cur[dt_bin,])) +
    mu_e[dt_bin,]*X_cur[dt_bin,] - log(1+exp(mu_e[dt_bin,]))*trials[dt_bin,] else numeric(0)

  bic_1    <- -2*(sum(ll_gau)+sum(ll_ber)+sum(ll_bin)) + log(n)*res$df
  bic_1obs <- -2*(sum(ll_gau[!mask_gau])+sum(ll_ber[!mask_ber])+sum(ll_bin[!mask_bin])) +
    log(n)*res$df

  out <- c(res, list(bic_1 = bic_1, bic_1obs = bic_1obs))

  if (!is.null(progress_file)) {
    write(c("",
      "--- Final results ---",
      sprintf("Factors   : %d  (started %d, pruned %d)", out$L, out$L_old, out$L_old-out$L),
      sprintf("Eff. df   : %d", out$df),
      sprintf("BIC       : %.2f", bic_1),
      sprintf("BIC (obs) : %.2f", bic_1obs)),
      file = progress_file, append = TRUE)
  }

  out
}


## user-facing wrapper ----

# X        : p x n data matrix (NA = missing)
# dt       : integer vector length p; 0=Gaussian 1=Bernoulli 2=Binomial
# modality : integer vector length p (1..H); overrides dt for nu assignment.
#            Use when multiple blocks share a data type but need separate nu_1/nu_2.
#            If NULL, nu_1/nu_2 are indexed by data type.
# nu_1     : prior mean on log-shrinkage alpha. Scalar, or length-H vector (per modality/type).
# nu_2     : prior variance on alpha. Same shape as nu_1.
# eta, eps : Wishart-like prior on Omega. eta=strength, eps=diagonal offset.
iBFA <- function(X, dt, graph = NULL, L = 5, T = 1000, burn_in = 500,
                  nu_1 = 1, nu_2 = 0.5, eta = 10, eps = 0.2,
                  trials = NULL, missing = NULL,
                  modality = NULL,
                  print_every = 100,
                  progress_file = NULL) {

  p <- nrow(X); n <- ncol(X)
  stopifnot(length(dt) == p, burn_in < T)

  mask    <- is.na(X)
  if (is.null(missing)) missing <- any(mask)
  if (is.null(trials))  trials  <- matrix(1L, p, n)
  if (is.null(graph))   graph   <- matrix(0, p, p)

  ## initialise Omega from graph structure
  omega_ini <- diag(1, p)
  ind <- which(matrix(as.logical(graph), p), arr.ind = TRUE)
  if (nrow(ind) > 0)
    for (i in seq_len(nrow(ind)))
      if (ind[i, 1] > ind[i, 2]) omega_ini[ind[i, 1], ind[i, 2]] <- 0.05
  omega_temp     <- as.matrix(forceSymmetric(omega_ini, uplo = "L"))
  inv_omega_temp <- solve(omega_temp)

  rho_ini <- matrix(1, p, n)
  tau_ini <- matrix(1, p, L)

  ## initialise alpha at its prior mean (per-feature if modality/type given)
  if (!is.null(modality) && length(nu_1) > 1) {
    alpha_ini <- matrix(nu_1[modality], p, L)
  } else if (length(nu_1) > 1) {
    alpha_ini <- matrix(nu_1[dt + 1], p, L)
  } else {
    alpha_ini <- matrix(nu_1, p, L)
  }

  w_ini   <- matrix(rnorm(L * p), p, L)
  z_ini   <- matrix(rnorm(L * n), L, n)
  Sigma   <- rep(0.1, p)
  Q       <- diag(4, p, p)

  box <- .empty_chain(n, p, T, L, rho_ini, w_ini, tau_ini, z_ini, alpha_ini, mask)

  w_temp     <- matrix(as.numeric(box$chain_w[1, ]),     p, L, byrow = TRUE)
  z_temp     <- matrix(as.numeric(box$chain_z[1, ]),     L, n, byrow = TRUE)
  rho_temp   <- matrix(as.numeric(box$chain_rho[1, ]),   p, n, byrow = TRUE)
  alpha_temp <- matrix(as.numeric(box$chain_alpha[1, ]), p, L, byrow = TRUE)
  tau_temp   <- matrix(as.numeric(box$chain_tau[1, ]),   p, L, byrow = TRUE)
  m_temp     <- rep(0, p)

  out <- iBFA_mcmc(
    dt = dt, T = T, L = L, p = p, n = n, X = X, trials = trials,
    nu_1 = nu_1, nu_2 = nu_2, Sigma = Sigma, Q = Q, eta = eta, eps = eps,
    rho_temp = rho_temp, tau_temp = tau_temp,
    omega_temp = omega_temp, inv_omega_temp = inv_omega_temp,
    w_temp = w_temp, z_temp = z_temp, alpha_temp = alpha_temp,
    m_temp = m_temp, start = burn_in + 1, end = T,
    chain_rho   = box$chain_rho,   chain_alpha = box$chain_alpha,
    chain_tau   = box$chain_tau,   chain_w     = box$chain_w,
    chain_z     = box$chain_z,     chain_m     = box$chain_m,
    chain_x     = box$chain_x,     missing     = missing,
    modality    = modality,
    print_every = print_every,     progress_file = progress_file
  )

  L_final   <- out$L
  out$w_est <- matrix(out$w_est, p, L_final, byrow = TRUE)
  out$z_est <- matrix(out$z_est, L_final, n, byrow = TRUE)
  out$nu_1  <- nu_1
  out$nu_2  <- nu_2
  out$mask  <- mask
  out
}
