#Helper functions

#Generalized inverse; port of MASS::ginv()
generalized_inverse <- function(sigma) {
  sigmasvd <- svd(sigma)
  pos <- sigmasvd$d > max(1e-7 * sigmasvd$d[1L], 0)
  sigma_inv <- sigmasvd$v[, pos, drop = FALSE] %*% (sigmasvd$d[pos]^-1 * t(sigmasvd$u[, pos, drop = FALSE]))
  sigma_inv[upper.tri(sigma_inv)] <- t(sigma_inv)[upper.tri(t(sigma_inv))]
  return(sigma_inv)
}

#Speed-optimized shortcut for diag(x) %*% A %*% diag(y)
xAy <- function(x, A, y) {
  tcrossprod(x, y) * A
}

#Weighted standard deviation
w_sd <- function(x, w) {
  sumw <- sum(w)
  sqrt(sum(w * (x - sum(x*w)/sumw) ^ 2) * sumw / (sumw^2 - sum(w^2)))
}

#Compute distance matrix
dist_mat <- function(X, Y = X, Sinv = NULL) {
  indx <- seq_len(nrow(X))
  indy <- nrow(X) + seq_len(nrow(Y))

  distances::distance_columns(distances::distances(rbind(X, Y), normalize = "none", weights = unname(Sinv)),
                              row_indices = indx, column_indices = indy)
}

wcov <- function(X, W = NULL, diag = FALSE) {
  if (is.null(W)) W <- rep(1, nrow(X))

  if (diag) {
    diag(apply(X, 2, w_sd, w = W)^2)
  }
  else {
    cov.wt(X, W)$cov
  }
}

#Process supplied weights
process_w <- function(W, Z = NULL, n = NULL) {
  if (is.null(Z)) {
    if (is.null(W)) W <- rep(1, n)
    else W <- W/mean(W)
  }
  else if (is.null(W)) W <- rep(1, length(Z))
  else {
    for (i in unique(Z)) W[Z==i] <- W[Z==i]/mean(W[Z==i])
  }

  W
}

#Process covariates
process_X <- function(X) {
  nm <- deparse(substitute(X))
  if (anyNA(X)) stop(paste0("NAs are not allowed in ", nm, "."), call. = FALSE)
  if (!is.matrix(X) && !is.data.frame(X)) stop(paste0(nm, " must be a data frame or matrix."), call. = FALSE)
  if (is.data.frame(X)) {
    if (any(sapply(X, is.factor))) X <- cobalt::splitfactor(X, drop.first = FALSE)
    X[sapply(X, is.logical)] <- lapply(X[sapply(X, is.logical)], as.numeric)
    X <- as.matrix(X)
  }
  if (is.character(X)) stop(paste0(nm, " must be a numeric matrix or data frame."), call. = FALSE)

  return(X)
}