#Compute the (weighted) energy distance between two distributions

edist <- function(sampleX, targetX = NULL, sampleW = NULL, targetW = NULL, Z = NULL, std = "studentized") {

  std <- match.arg(std, c("studentized", "mahalanobis", "none"))

  sampleX <- process_X(sampleX, std)

  if (is.null(targetX)) targetX <- sampleX
  else targetX <- process_X(targetX, std)

  check_vars(sampleX, targetX)

  if (is.null(Z)) {
    n1 <- nrow(sampleX)
    n0 <- nrow(targetX)

    w1 <- process_w(sampleW, n = n1)
    w0 <- process_w(targetW, n = n0)

    check_lengths(sampleX, sampleW)
    check_lengths(targetX, targetW)

    if (std == "studentized") {
      sds <- apply(targetX, 2, w_sd, w = w0)
      sampleX <- scale(sampleX, scale = sds, center = FALSE)
      targetX <- scale(targetX, scale = sds, center = FALSE)
    }
    else if (std == "mahalanobis") {
      Sinv <- generalized_inverse(wcov(targetX, w0, diag = FALSE))
      ch <- chol2(Sinv)
      sampleX <- tcrossprod(sampleX, ch)
      targetX <- tcrossprod(targetX, ch)
    }

    d_1_1 <- dist_mat(sampleX)   #sample x sample
    d_0_0 <- dist_mat(targetX) #target x target
    d_1_0 <- dist_mat(sampleX, targetX) #sample x target

  }
  else {
    Z <- process_Z(Z, bin = TRUE)

    n1 <- sum(Z==1)
    n0 <- sum(Z==0)

    w <- process_w(sampleW, Z = Z)
    w1 <- w[Z==1]
    w0 <- w[Z==0]

    check_lengths(sampleX, sampleW, Z, targetW)

    if (is.null(targetW)) targetW <- sampleW

    tw <- process_w(targetW, n = nrow(targetX))

    if (std == "studentized") {
      sds <- apply(targetX, 2, w_sd, w = tw)
      sampleX <- scale(sampleX, scale = sds, center = FALSE)
    }
    else if (std == "mahalanobis") {
      Sinv <- generalized_inverse(wcov(targetX, tw, diag = FALSE))
      ch <- chol2(Sinv)
      sampleX <- tcrossprod(sampleX, ch)
    }

    sample1 <- sampleX[Z==1,,drop = FALSE]
    sample0 <- sampleX[Z==0,,drop = FALSE]

    d_1_1 <- dist_mat(sample1)   #sample x sample
    d_0_0 <- dist_mat(sample0) #target x target
    d_1_0 <- dist_mat(sample1, sample0) #sample x target
  }
  #Formula for energy distance
  drop(
    (2/(n1*n0)) * (w1 %*% d_1_0 %*% w0) - (1/n1^2) * (w1 %*% d_1_1 %*% w1) - (1/n0^2) * (w0 %*% d_0_0 %*% w0)
  )
}
