eb_ate <- function(sampleX, Z, sampleW = NULL, std = "studentized", improved = TRUE) {
  energybalance(sampleX, Z = Z, sampleW = sampleW, std = std, improved = improved)
}

eb_target <- function(sampleX, targetX, sampleW = NULL, targetW = NULL, std = "studentized") {
  energybalance(sampleX, targetX = targetX, sampleW = sampleW, targetW = targetW, std = std)
}

eb_mediation <- function(sampleX, Z, M, sampleW = NULL, std = "studentized") {
  sampleWc <- process_w(sampleW, Z)[Z==0]

  w_ate_c <- eb_target(sampleX[Z==0,], targetX = sampleX, sampleW = sampleWc,
                       targetW = sampleWc, std = std)

  sampleXM <- cbind(sampleX, M)

  eb_target(sampleXM[Z==1,], targetX = sampleXM[Z==0,], sampleW = sampleW,
            targetW = w_ate_c * sampleWc, std = std)
}