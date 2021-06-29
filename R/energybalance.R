#Estimates energy balancing weights to match distributions

energybalance <- function(sampleX, Z = NULL, targetX = NULL, sampleW = NULL, targetW = NULL, std = "studentized", improved = TRUE) {

  sampleX <- process_X(sampleX)

  std <- match.arg(std, c("studentized", "mahalanobis", "none"))

  n_s <- nrow(sampleX)

  if (is.null(Z)) {
    Z <- rep(0L, n_s)
    z_levels <- 0L
    nz <- 1L
    improved <- FALSE
  } else {
    z_levels <- unique(Z)
    nz <- length(z_levels)
  }

  sw <- process_w(sampleW, Z, n = n_s)

  if (is.null(targetX)) {
    targetX <- sampleX
    tw <- sw
  }
  else {
    targetX <- process_X(targetX)
  }

  n_t <- nrow(targetX)
  tw <- process_w(targetW, n = n_t)

  #Distance matrices
  #Standardize variables based on target sample
  Sinv <- switch(std,
                 "studentized" = generalized_inverse(wcov(targetX, tw, diag = TRUE)),
                 "mahalanobis" = generalized_inverse(wcov(targetX, tw, diag = FALSE)),
                 "none" = NULL)

  min.w <- 1e-8

  P <- matrix(0, nrow = n_s, ncol = n_s)
  Q <- rep(0, n_s)
  A <- matrix(0, nrow = n_s + nz, ncol = n_s)
  l <- rep(min.w, n_s + nz)
  u <- rep(Inf, n_s + nz)

  diag(A)[seq_len(n_s)] <- 1

  for (i in seq_len(nz)) {
    in_zi <- which(Z == z_levels[i])
    n_zi <- length(in_zi)

    d_zi_zi <- dist_mat(sampleX[in_zi,,drop = FALSE], Sinv = Sinv)

    d_zi_t <- dist_mat(sampleX[in_zi,,drop = FALSE], targetX, Sinv = Sinv)

    sw_zi <- sw[in_zi]/n_zi

    if (improved) {
      P[in_zi, in_zi] <- -nz * xAy(sw_zi, d_zi_zi, sw_zi)
    } else {
      P[in_zi, in_zi] <- -1 * xAy(sw_zi, d_zi_zi, sw_zi)
    }
    Q[in_zi] <- 2 * ((sw_zi * d_zi_t) %*% tw/n_t)

    A[n_s + i, in_zi] <- sw[in_zi]
    l[n_s + i] <- n_zi
    u[n_s + i] <- n_zi
  }

  if (improved) {
    z_combns <- utils::combn(seq_len(nz), 2, simplify = FALSE)

    for (t_ in z_combns) {
      in_z1 <- which(Z == z_levels[t_[1]])
      in_z2 <- which(Z == z_levels[t_[2]])

      n_z1 <- length(in_z1)
      n_z2 <- length(in_z2)

      sw_z1 <- sw[in_z1]/n_z1
      sw_z2 <- sw[in_z2]/n_z2

      d_z1_z2 <- dist_mat(sampleX[in_z1,,drop = FALSE],
                          sampleX[in_z2,,drop = FALSE],
                          Sinv = Sinv)

      P[in_z1, in_z2] <- P[in_z1, in_z2] + xAy(sw_z1, d_z1_z2, sw_z2)
      P[in_z2, in_z1] <- t(P[in_z1, in_z2])
    }
  }

  pars <- osqp::osqpSettings(verbose = FALSE,
                             max_iter = 2000L,
                             eps_abs = 1e-8,
                             eps_rel = 1e-8)
  opt.out <- osqp::solve_osqp(2 * P, Q, A, l, u, pars = pars)
  w <- opt.out$x
  w[w < min.w] <- min.w

  return(w)
}
