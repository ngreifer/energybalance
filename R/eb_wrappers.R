eb_ate <- function(sampleX, Z, sampleW = NULL, std = "studentized", improved = TRUE, lambda = 0) {
  energybalance(sampleX, Z = Z, sampleW = sampleW, std = std, improved = improved,
                lambda = lambda)
}

eb_target <- function(sampleX, targetX, sampleW = NULL, targetW = NULL, std = "studentized", lambda = 0) {
  energybalance(sampleX, targetX = targetX, sampleW = sampleW, targetW = targetW,
                std = std, lambda = lambda)
}

eb_mediation <- function(sampleX, Z, M, sampleW = NULL, std = "studentized", lambda = 0) {
  Z <- process_Z(Z, bin = TRUE)
  M <- process_X(M)

  sampleW <- process_w(sampleW, Z)

  check_lengths(Z, sampleX, M, sampleW)

  w_ate_c <- eb_target(sampleX[Z==0,], targetX = sampleX, sampleW = sampleW[Z==0],
                       targetW = sampleW, std = std, lambda = lambda)

  sampleXM <- cbind(sampleX, M)

  eb_target(sampleXM[Z==1,], targetX = sampleXM[Z==0,], sampleW = sampleW[Z==1],
            targetW = w_ate_c * sampleW[Z==0], std = std, lambda = lambda)
}

eb_att <- function(sampleX, Z, sampleW = NULL, std = "studentized", focal = NULL, improved = FALSE, lambda = 0) {
  Z <- process_Z(Z)
  focal <- process_focal(focal, Z)

  sampleW <- process_w(sampleW, Z)

  sampleX <- process_X(sampleX)

  w <- rep(1, length(Z))

  e <- energybalance(sampleX[Z != focal,,drop = FALSE], Z = Z[Z != focal], targetX = sampleX[Z == focal,,drop = FALSE],
                     sampleW = sampleW[Z != focal], targetW = sampleW[Z == focal],
                     std = std, improved = improved, lambda = lambda)

  w[Z != focal] <- e

  return(w)
}