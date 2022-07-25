#Compute the (weighted) energy distance between two distributions

edist <- function(sampleX, targetX = NULL, sampleW = NULL, targetW = NULL, Z = NULL, std = "studentized") {

  std <- match.arg(std, c("studentized", "mahalanobis", "none"))

  sampleX <- process_X(sampleX)

  if (is.null(targetX)) targetX <- sampleX
  else targetX <- process_X(targetX)

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
  # drop(
  #   (2/(n1*n0)) * (w1 %*% d_1_0 %*% w0) - (1/n1^2) * (w1 %*% d_1_1 %*% w1) - (1/n0^2) * (w0 %*% d_0_0 %*% w0)
  # )

  2 * quad_mult(w1, d_1_0, w0)/(n1*n0) - quad_mult(w1, d_1_1, w1)/(n1*n1) - quad_mult(w0, d_0_0, w0)/(n0*n0)
}

edist_ate <- function(sampleX, Z, sampleW = NULL, targetW = NULL, std = "studentized", improved = TRUE) {

  std <- match.arg(std, c("studentized", "mahalanobis", "none"))

  sampleX <- process_X(sampleX)

  n_s <- nrow(sampleX)

  Z <- process_Z(Z)
  z_levels <- sort(unique(Z))
  nz <- length(z_levels)

  sampleW <- process_w(sampleW, Z, n = n_s)

  targetX <- sampleX

  n_t <- nrow(targetX)
  targetW <- process_w(targetW, n = n_t)

  check_lengths(sampleX, Z, sampleW, targetW)

  sw <- sampleW
  tw <- targetW

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

  d_t <- dist_mat(targetX)

  ed <- setNames(numeric(nz), paste(z_levels, "vs. full"))

  for (i in seq_len(nz)) {
    in_zi <- which(Z == z_levels[i])
    n_zi <- length(in_zi)

    d_zi_zi <- dist_mat(sampleX[in_zi,,drop = FALSE])

    d_zi_t <- dist_mat(sampleX[in_zi,,drop = FALSE], targetX)

    sw_zi <- sw[in_zi]

    ed[i] <- 2 * quad_mult(sw_zi, d_zi_t, tw)/(n_zi*n_t) - quad_mult(sw_zi, d_zi_zi, sw_zi)/(n_zi*n_zi) - quad_mult(tw, d_t, tw)/(n_t*n_t)
  }

  if (improved) {
    z_combns <- utils::combn(seq_len(nz), 2, simplify = FALSE)
    ed_imp <- setNames(numeric(length(z_combns)), sapply(z_combns, function(t_) paste(z_levels[t_], collapse = " vs. ")))

    d_z1_z1 <- d_z2_z2 <- NULL
    for (i_ in seq_along(z_combns)) {
      t_ <- z_combns[[i_]]
      in_z1 <- which(Z == z_levels[t_[1]])
      in_z2 <- which(Z == z_levels[t_[2]])

      n_z1 <- length(in_z1)
      n_z2 <- length(in_z2)

      sw_z1 <- sw[in_z1]
      sw_z2 <- sw[in_z2]

      d_z1_z2 <- dist_mat(sampleX[in_z1,,drop = FALSE],
                          sampleX[in_z2,,drop = FALSE])

      if (i_ > 1) {
        if (t_[1] == z_combns[[i_-1]][1]) d_z1_z1 <- d_z1_z1
        else if (t_[1] == z_combns[[i_-1]][2]) d_z1_z1 <- d_z2_z2
        else d_z1_z1 <- dist_mat(sampleX[in_z1,,drop = FALSE])

        if (t_[2] == z_combns[[i_-1]][1]) d_z2_z2 <- d_z1_z1
        else if (t_[1] == z_combns[[i_-1]][2]) d_z2_z2 <- d_z2_z2
        else d_z2_z2 <- dist_mat(sampleX[in_z2,,drop = FALSE])
      }
      else {
        d_z1_z1 <- dist_mat(sampleX[in_z1,,drop = FALSE])
        d_z2_z2 <- dist_mat(sampleX[in_z2,,drop = FALSE])
      }

      ed_imp[i_] <- 2 * quad_mult(sw_z1, d_z1_z2, sw_z2)/(n_z1*n_z2) - quad_mult(sw_z1, d_z1_z1, sw_z1)/(n_z1*n_z1) - quad_mult(sw_z2, d_z2_z2, sw_z2)/(n_z2*n_z2)
    }
    ed <- c(ed, ed_imp)
  }

  total_edist <- sum(ed)
  attr(total_edist, "components") <- ed

  return(total_edist)

}