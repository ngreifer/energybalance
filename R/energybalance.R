#Estimates energy balancing weights to match distributions

energybalance <- function(sampleX, Z = NULL, targetX = NULL, sampleW = NULL, targetW = NULL, std = "studentized", improved = TRUE,
                          lambda = 0) {

  std <- match.arg(std, c("studentized", "mahalanobis", "none"))

  sampleX <- process_X(sampleX)

  n_s <- nrow(sampleX)

  if (is.null(Z)) {
    Z <- rep(0L, n_s)
    z_levels <- 0L
    nz <- 1L
    improved <- FALSE
  } else {
    Z <- process_Z(Z)
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

  check_lengths(sampleX, Z, sampleW)
  check_lengths(targetX, targetW)
  check_vars(sampleX, targetX)

  #Distance matrices
  #Standardize variables based on target sample
  if (std == "studentized") {
    sds <- apply(targetX, 2, w_sd, w = tw)
    sampleX <- scale(sampleX, scale = sds, center = FALSE)
    targetX <- scale(targetX, scale = sds, center = FALSE)
  }
  else if (std == "mahalanobis") {
    Sinv <- generalized_inverse(wcov(targetX, tw, diag = FALSE))
    ch <- chol2(Sinv)
    sampleX <- tcrossprod(sampleX, ch)
    targetX <- tcrossprod(targetX, ch)
  }

  min.w <- 1e-8

  P <- matrix(0, nrow = n_s, ncol = n_s)
  Q <- rep(0, n_s)

  #Constraint matrix; sum of weights in each Z
  A <- matrix(0, nrow = nz, ncol = n_s)
  eq <- rep(0, nz)

  #LB and LB for weights
  lb <- rep(min.w, n_s)
  ub <- rep(Inf, n_s)

  for (i in seq_len(nz)) {
    in_zi <- which(Z == z_levels[i])
    n_zi <- length(in_zi)

    d_zi_zi <- dist_mat(sampleX[in_zi,,drop = FALSE])

    d_zi_t <- dist_mat(sampleX[in_zi,,drop = FALSE], targetX)

    sw_zi <- sw[in_zi]/n_zi

    if (improved) {
      P[in_zi, in_zi] <- -nz * xAy(sw_zi, d_zi_zi, sw_zi)
    } else {
      P[in_zi, in_zi] <- -1 * xAy(sw_zi, d_zi_zi, sw_zi)
    }
    Q[in_zi] <- 2 * ((sw_zi * d_zi_t) %*% tw/n_t)

    A[i, in_zi] <- sw[in_zi]
    eq[i] <- n_zi
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
                          sampleX[in_z2,,drop = FALSE])

      P[in_z1, in_z2] <- P[in_z1, in_z2] + xAy(sw_z1, d_z1_z2, sw_z2)
      P[in_z2, in_z1] <- t(P[in_z1, in_z2])
    }
  }

  if (length(lambda) == 0) lambda <- 0
  else if (length(lambda) > 1 || !is.numeric(lambda) || lambda < 0) {
    stop("'lambda' must be a single non-negative number.", call. = FALSE)
  }
  if (lambda > 0) {
    P <- P + (lambda/n_s^2) * diag(n_s)
  }

  w <- quad_solve(solver = "osqp", P, Q, A, eq, lb, ub)
  w[w < min.w] <- min.w

  return(w)
}
