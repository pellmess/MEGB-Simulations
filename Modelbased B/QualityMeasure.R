


QualityMeasure <- function(True_mean, Est_mean, MSE, MSETF) {
  ResultsAnaModelBasedMSE <- function (True_mean, Est_mean, MSE) {
    warning("MSE Version noch nicht korrigiert")
    m <- dim(True_mean)[1]
    NoSim <- dim(True_mean)[2]
    True_mean <- t(True_mean)
    Est_mean <- t(Est_mean)
    MSE <- t(MSE)
    
    min_runs <- min(nrow(True_mean), nrow(Est_mean)) # Bestimme wie viele runs
    # minimal funktionierten, um nur anhand dieser den MSE zu berechnen
    True_mean <- True_mean[1:min_runs, ]
    MSE <- MSE[1:min_runs, ]
    
    assign("m", m, pos = 1)
    RB <- rep(0, m)
    Bias <- rep(0, m)
    RRMSE <- rep(0, m)
    
    True.MSE <- rep(0, m)
    Est.MSE <- rep(0, m)
    CV.MSE <- rep(0, m)
    BIAS.MSE <- rep(0, m)
    
    for (i in 1:m) {
      RB[i] <- mean(((Est_mean[, i]) - (True_mean[, i])) / (True_mean[, i]), na.rm = TRUE)
      Bias[i] <- mean(((Est_mean[, i]) - (True_mean[, i])), na.rm = TRUE)
      RRMSE[i] <-
        sqrt(mean((((Est_mean[, i]) - (True_mean[, i])
        ) / (True_mean[, i])) ^ 2, na.rm = TRUE))
      True.MSE[i] <- mean((Est_mean[, i] - True_mean[, i]) ^ 2, na.rm = TRUE)
      Est.MSE[i] <- mean(MSE[, i], na.rm = TRUE)
      
      CV.MSE[i] <-
        (sqrt(mean((
          sqrt(MSE[, i]) - sqrt(True.MSE[i])
        ) ^ 2, na.rm = TRUE)) / sqrt(True.MSE[i])) * 100
      BIAS.MSE[i] <- mean(MSE[, i] - True.MSE[i], na.rm = TRUE) / True.MSE[i]
    }
    True.RootMSE = sqrt(True.MSE)
    Est.RootMSE = sqrt(Est.MSE)
    RB.estRRMSE = ((Est.RootMSE - True.RootMSE) / True.RootMSE)
    
    CR <- rep(0, m)
    CR_Length <- rep(0, m)
    coverage <- matrix(0, NoSim, m)
    Length <- matrix(0, NoSim, m)
    
    
    for (k in 1:NoSim)
    {
      for (oo in 1:m)
      {
        Length[k, oo] <- 2 * 1.96 * sqrt(MSE[k, oo])
        if ((
          Est_mean[k, oo] - 1.96 * sqrt(MSE[k, oo]) < True_mean[k, oo] &&
          True_mean[k, oo] < Est_mean[k, oo] + 1.96 * sqrt(MSE[k, oo])
        ))
          coverage[k, oo] <- 1
        else
          coverage[k, oo] <- 0
      }
    }
    
    for (i in 1:m) {
      CR[i] <- mean(coverage[, i])
      CR_Length[i] <- mean(Length[, i])
    }
    
    
    list(
      Bias = Bias,
      RB = RB,
      RRMSE = RRMSE,
      True.RMSE = True.RootMSE,
      Est.RMSE = Est.RootMSE,
      RB_RMSE = RB.estRRMSE,
      RRMSE_RMSE = CV.MSE,
      Coverage = CR,
      CI_Length = CR_Length
    )
  }
  
  ResultsAnaModelBased <- function (True_mean, Est_mean) {
    m <- dim(True_mean)[1] # Anzahl an Areas
    NoSim <- dim(True_mean)[2] # Anzahl an Simulationen, basierend auf True Mean
    True_mean <- t(True_mean) # m x NoSim Matrix
    Est_mean <- t(Est_mean)
    
    min_runs <- min(nrow(True_mean), nrow(Est_mean)) # Bestimme wie viele runs
    # minimal funktionierten, um nur anhand dieser den MSE zu berechnen
    True_mean <- True_mean[1:min_runs, ]
    assign("m", m, pos = 1)
    RB <- rep(0, m)
    ARB <- rep(0, m)
    RRMSE <- rep(0, m)
    Bias <- rep(0, m)
    True <- rep(0, m)
    Dir <- rep(0, m)
    
    True.MSE <- rep(0, m)
    Est.MSE <- rep(0, m)
    CV.MSE <- rep(0, m)
    BIAS.MSE <- rep(0, m)
    if(anyNA(Est_mean)){
      
    warning("NA's werden ignoiert, wahrscheinlich nicht alle duchgelaufen")  
    }
    for (i in 1:m) {
      Bias[i] <- mean((Est_mean[, i]) - (True_mean[, i]), na.rm = TRUE)
      True[i] <- mean(True_mean[, i], na.rm = TRUE)
      Dir[i] <- mean(Est_mean[, i], na.rm = TRUE)
      RB[i] <- mean(((Est_mean[, i]) - (True_mean[, i])) / True_mean[, i], na.rm = TRUE)
      ARB[i] <-
        mean(abs(((Est_mean[, i]) - (True_mean[, i]))) / abs(True_mean[, i]), na.rm = TRUE)
      RRMSE[i] <-
        sqrt(mean((((Est_mean[, i]) - (True_mean[, i])
        ) / (True_mean[, i])) ^ 2, na.rm = TRUE))
      True.MSE[i] <- mean((Est_mean[, i] - True_mean[, i]) ^ 2, na.rm = TRUE)
      warning("na.rm ERGÄNZT")
    }
    
    True.RootMSE = sqrt(True.MSE)
    
    
    
    list(
      Bias = Bias,
      True = True,
      Dir = Dir,
      RB = RB,
      RRMSE = RRMSE,
      True.RMSE = True.RootMSE
    )
  }
  
  if (MSETF == TRUE) {
    results <- ResultsAnaModelBasedMSE(True_mean, Est_mean, MSE)
    return(results)
  }
  
  if (MSETF == FALSE) {
    results <- ResultsAnaModelBased(True_mean, Est_mean)
    return(results)
  }
}
################################################################################

# Zum Test (Yeonjoo)
#load("res.Rdata")
#res <- QualityMeasure(True_mean = True_Mean, Est_mean = Mean_Full, MSETF = F)

#Bias_Mean$Bias_Mean[Bias_Mean$Mod == "Full"]
#sqrt(MSE_Mean$MSE_Mean[MSE_Mean$Mod == "Full"])

#all.equal(res$Bias, Bias_Mean$Bias_Mean[Bias_Mean$Mod == "Full"])
#all.equal(res$True.RMSE, sqrt(MSE_Mean$MSE_Mean[MSE_Mean$Mod == "Full"]))
#
# load("Res_EBP.RData")
# res <- QualityMeasure(True_mean = Mean, Est_mean = Mean_Full, MSETF = F)
#
# res$Bias
# rowMeans(Bias_Mean_Full)
# all.equal(res$Bias, rowMeans(Bias_Mean_Full))
# res$True.RMSE
# sqrt(rowMeans(MSE_Mean_Full))
# all.equal(res$True.RMSE, sqrt(rowMeans(MSE_Mean_Full)))
