#Helper functions

#Generalized inverse; port of MASS::ginv()
generalized_inverse <- function(sigma) {
  sigmasvd <- svd(sigma)
  pos <- sigmasvd$d > max(1e-9 * sigmasvd$d[1L], 0)
  sigma_inv <- sigmasvd$v[, pos, drop = FALSE] %*% (sigmasvd$d[pos]^-1 * t(sigmasvd$u[, pos, drop = FALSE]))
  sigma_inv[upper.tri(sigma_inv)] <- t(sigma_inv)[upper.tri(t(sigma_inv))]

  dimnames(sigma_inv) <- dimnames(sigma)
  return(sigma_inv)
}

#Choleski decomp, works for nonpositive definite matrices too
chol2 <- function(Sinv) {
  ch <- suppressWarnings(chol(Sinv, pivot = TRUE))
  p <- order(attr(ch, "pivot"))
  return(ch[,p])
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

#Computes v1 %*% M %*% v2 efficiently
quad_mult <- function(v1, M, v2 = v1) {
  tcrossprod(v1, tcrossprod(v2, M))
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

#Process covariates (and M)
process_X <- function(X, std = "none") {
  nm <- paste0("'", deparse(substitute(X)), "'")
  if (is.atomic(X) && length(dim(X)) == 0) X <- data.frame(X)
  else if (!is.matrix(X) && !is.data.frame(X)) stop(paste0(nm, " must be a data frame or matrix."), call. = FALSE)
  if (anyNA(X)) stop(paste0("NAs are not allowed in ", nm, "."), call. = FALSE)
  if (is.data.frame(X)) {
    logicals <- which(vapply(X, is.logical, logical(1L)))
    X[logicals] <- lapply(X[logicals], as.numeric)
    if (any(sapply(X, function(x) is.factor(x) || is.character(x))))
      X <- cobalt::splitfactor(X, drop.first = "if2")
    X <- as.matrix(X)
  }
  if (is.character(X)) stop(paste0(nm, " must be a numeric matrix or data frame."), call. = FALSE)

  return(X)
}

#Process Z
process_Z <- function(Z, bin = FALSE) {
  nm <- paste0("'", deparse(substitute(Z)), "'")
  if (!is.atomic(Z) || length(dim(Z)) > 0) stop(paste0(nm, " must be an atomic vector."), call. = FALSE)
  if (anyNA(Z)) stop(paste0("NAs are not allowed in ", nm, "."), call. = FALSE)

  if (bin) {
    if (length(unique(Z)) != 2) stop(paste0(nm, " must have exactly 2 unique values."), call. = FALSE)
    if (is.logical(Z)) Z <- as.numeric(Z)
    else if (is.factor(Z)) Z <- as.numeric(Z == levels(Z)[2])
    else Z <- as.numeric(Z == max(Z))
  }

  return(Z)
}

#Process focal for eb_att
process_focal <- function(focal = NULL, Z) {
  if (is.null(focal)) {
    if (is.factor(Z) && levels(Z)[nlevels(Z)] %in% Z) focal <- levels(Z)[nlevels(Z)]
    else focal <- max(Z)
  }
  else {
    if (length(focal) != 1) {
      stop("'focal' must be a single value.", call. = FALSE)
    }
    if (!focal %in% Z) {
      stop("'focal' must be a value in Z.", call. = FALSE)
    }
  }

  return(focal)
}

#Check to ensure lengths are correct
check_lengths <- function(...) {
  nm <- unlist(lapply(substitute(list(...))[-1], deparse1))

  lens <- unlist(lapply(list(...), function(x) {
    if (length(x) == 0) 0L
    else if (length(dim(x)) == 0) length(x)
    else nrow(x)
  }))

  nm <- nm[lens > 0]
  lens <- lens[lens > 0]
  if (length(lens) > 1) {
    bad_lens <- which(lens != lens[1])
    if (length(bad_lens) > 0) {
      stop(paste0("The following argument(s) must have the same number of units as '", nm[1],"': ", paste0(nm[bad_lens], collapse = ", ")),
           call. = FALSE)
    }
  }
}

#Check to ensure same variables are used
check_vars <- function(sampleX, targetX) {

  vn <- lapply(list(sampleX, targetX), function(x) {
    n <- colnames(x)
    if (length(n) == 0) n <- seq_len(ncol(x))
    n
  })

  vn_len <- lengths(vn)

  if (any(vn_len != vn_len[1])) {
    stop("'sampleX' and 'targetX' must have the same number of columns.", call. = FALSE)
  }

  char_nm <- which(vapply(vn, is.character, logical(1L)))

  if (length(char_nm) >= 2) {
    if (!all(vn[[char_nm[1]]] %in% unlist(vn[char_nm]))) {
      stop("'sampleX' and 'targetX' must have the same variable names.", call. = FALSE)
    }
  }

}

#Solve the qudratic program using a solver
quad_solve <- function(solver = "osqp", P, Q, A, eq, lb, ub) {
  if (solver == "osqp") {
    A <- rbind(A, diag(length(lb)))
    l <- c(eq, lb)
    u <- c(eq, ub)

    pars <- osqp::osqpSettings(verbose = FALSE,
                               max_iter = 2000L,
                               eps_abs = 1e-8,
                               eps_rel = 1e-8)
    for (i in 1:15) {
      #May need to run multiple times to get a solution
      opt.out <- osqp::solve_osqp(P = 2 * P, q = Q, A = A, l = l, u = u, pars = pars)
      w <- opt.out$x

      if (abs(max(w) - min(w)) > sqrt(.Machine$double.eps)) {
        break
      }
      else if (i == 15) {
        warning("The optimization failed to converge. Try running it again.", call. = FALSE)
        w[] <- 1
      }
    }
  }

  return(w)
}

if (FALSE) {
  for (i in dir("R/")) source(paste0("R/", i))
}