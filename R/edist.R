#Compute the (weighted) energy distance between two distributions

edist <- function(sampleX, targetX = sampleX, sampleW = NULL, targetW = NULL, Z = NULL, std = "studentized") {

  sampleX <- process_X(sampleX)
  targetX <- process_X(targetX)

  std <- match.arg(std, c("studentized", "mahalanobis", "none"))

  if (is.null(Z)) {
    n1 <- nrow(sampleX)
    n0 <- nrow(targetX)

    w1 <- process_w(sampleW, n = n1)
    w0 <- process_w(targetW, n = n0)

    Sinv <- switch(std,
                   "studentized" = generalized_inverse(wcov(targetX, w0, diag = TRUE)),
                   "mahalanobis" = generalized_inverse(wcov(targetX, w0, diag = FALSE)),
                   "none" = NULL)

    d_1_1 <- dist_mat(sampleX, Sinv = Sinv)   #sample x sample
    d_0_0 <- dist_mat(targetX, Sinv = Sinv) #target x target
    d_1_0 <- dist_mat(sampleX, targetX, Sinv = Sinv) #sample x target

  }
  else {

    if (length(unique(Z)) != 2) {
      stop("'Z' must have exactly 2 unique values.", call. = FALSE)
    }

    i <- Z[1]
    n1 <- sum(Z==i)
    n0 <- sum(Z!=i)

    w <- process_w(sampleW, Z = Z)
    w1 <- w[Z==i]
    w0 <- w[Z!=i]

    if (is.null(targetW)) targetW <- sampleW

    tw <- process_w(targetW, n = nrow(targetX))

    Sinv <- switch(std,
                   "studentized" = generalized_inverse(wcov(targetX, tw, diag = TRUE)),
                   "mahalanobis" = generalized_inverse(wcov(targetX, tw, diag = FALSE)),
                   "none" = NULL)

    sample1 <- sampleX[Z==i,,drop = FALSE]
    sample0 <- sampleX[Z!=i,,drop = FALSE]

    d_1_1 <- dist_mat(sample1, Sinv = Sinv)   #sample x sample
    d_0_0 <- dist_mat(sample0, Sinv = Sinv) #target x target
    d_1_0 <- dist_mat(sample1, sample0, Sinv = Sinv) #sample x target
  }

  #Formula for energy distance
  drop(
    (2/(n1*n0)) * (w1 %*% d_1_0 %*% w0) - (1/n1^2) * (w1 %*% d_1_1 %*% w1) - (1/n0^2) * (w0 %*% d_0_0 %*% w0)
  )
}
