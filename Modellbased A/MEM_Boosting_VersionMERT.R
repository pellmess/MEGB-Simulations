################################################################################
# MEM boosting - Version MERT
################################################################################


#===============================================================================
# get function for gradient tree boosting

source("gtb.R") 
#===============================================================================


#===============================================================================
# function for checking convergence

mem_boost_gll <- function( Z = NULL, ID = NULL, bhat = NULL,
                           ehat = NULL, UniqueID = NULL, NID = NULL,
                           D = NULL, Sigma2 = NULL )
{
  #- define ll:
  ll <- 0
  #- compute inverse and determinant of D:
  invD <- solve( D )
  logD <- as.numeric( determinant( D, log = TRUE )$modulus )
  #- let's go:
  for ( ii in 1:NID ) {
    #- get index variable:
    idx <- which( ID == UniqueID[ii] ) 
    #- get person data:
    ei <- ehat[idx]
    Zi <- Z[idx,,drop = FALSE]
    ni <- dim( Zi )[1]
    Ri <- diag( as.numeric( Sigma2 ), ni )
    #- compute inverse matrix and determinant:
    InvRi <- solve( Ri )
    logR <- as.numeric( determinant( Ri, log = TRUE )$modulus )
    #- compute ll:
    ll <- ll + logD + logR + t(ei)%*%InvRi%*%(ei) + t(bhat[ii,])%*%invD%*%bhat[ii,]
  }
  return(-ll)
}


#===============================================================================

# MERT Boosting

mem_boost_mert <- function( formula, data = NULL, random = NULL,
                            conv_memboost = 0.001, 
                            maxIter_memboost = 100, minIter_memboost = 0,
                            verbose_memboost = FALSE,
                            n.trees = 100, interaction.depth = 1,
                            shrinkage = 0.1, loss = c( "L2", "L1", "Huber" ),
                            modification = c( "none", "beta"), 
                            minsplit = 20, 
                            bag.fraction = 0.5, 
                            seed = NULL,
                            cp = 0.01 ){
  
  
  if(missing(modification)) { modification <- "none" }

  
  if (loss == "L2" & modification != "none") {
    stop('modification only used for L1 or Huber loss')
  }
  
  #===============================================================================
  # STEP 0: PREPARATION
  #===============================================================================
 
  #- Get X
  PredNames <- attr( stats::terms( formula ), "term.labels" )
  X <- model.frame(terms(reformulate(PredNames)), data = data)
  
  #- Get Y
  OutcomeName <- formula[[2]]
  if( length( OutcomeName ) > 1 ) { OutcomeName <- OutcomeName[3] } 
  Y <- data[ , toString( OutcomeName )]
  
  #- Get ID and Z
  FormulaRandom <- random
  NamesRandom <- attr( stats::terms( FormulaRandom ), "term.labels" )
  NamesRandom <- gsub("\\s", "", NamesRandom ) # delete all spaces 
  HasBar <- grepl("|", NamesRandom, fixed=TRUE) 
  if ( any( !HasBar ) ) { stop("'random' must contain a grouping variable after the | symbol.") }
  FormulaRandomSplit <- strsplit( NamesRandom, "\\|", perl = FALSE )[[1]]
  IdVar <- FormulaRandomSplit[2]
  if ( !( IdVar %in% colnames( data ) ) ) { stop("Level-2 identifier not found.") }
  ID <- data[ , toString( IdVar )]    
  FormulaRandom <- stats::formula( paste0("~", FormulaRandomSplit[1], collapse="") ) 
  Z <- model.matrix( FormulaRandom, data )
  
  #- Some initial specifications
  TotalObs <- dim( data )[1]
  UniqueID <- unique( ID )
  NID <- length( UniqueID )
  p <- dim(Z)[2]
  Dhat <- diag( 1, p ) 
  Sigma2hat <- 1  
  bhat <- matrix( 0, nrow = NID , ncol = p )    
  ehat <- rep( 0, TotalObs )
  
  # Some preparations for saving the means of the transformed outcome,
  # of the random intercept and of the boosting ensemble predictions
  # as well as the estimated covariance matrix of the random effects
  # and the estimated error variance per iteration 
  means.Ystar <- NULL
  means.fhat <- NULL
  means.ranint <- NULL
  DhatList <- list()
  errorVarList <- list()
  
  # Dataframe and initializations for the while-loop:
  newdata <- data
  toIterate <- TRUE
  convWarning <- FALSE 
  noIterations <- 0
  llnew <- 0
  
  #- step 5: Start the while loop    
  while ( toIterate ) {
    
    #- Count number of iterations
    noIterations <- noIterations + 1
    llold <- llnew
    
    #===============================================================================
    # STEP 1a: Get an estimate f
    #===============================================================================
    #----------------------------------------------------
    #- (i): compute the transformed outcome Ystar 
    #----------------------------------------------------
    Ystar <- rep( 0, TotalObs )
    for ( ii in 1:NID ) {
      #- get index variable:
      idx <- which( ID == UniqueID[ii] )
      #- get relevant matrices and vectors:
      Yi <- Y[idx]
      Zi <- Z[idx,,drop = FALSE]
      bi <- bhat[ii,]
      Ystar[idx] <- Yi - Zi%*%bi 
    }
    
    meanYstar <- mean(Ystar)
    means.Ystar <- rbind(means.Ystar, meanYstar)
    
    #----------------------------------------------------
    #- (ii): estimate f via gradient tree boosting 
    #----------------------------------------------------
    newdata[,"Ystar"] <- Ystar
    formula <- update.formula(formula, as.formula('Ystar ~ .'))
    tmpGTB <- gtb( formula, data = newdata,
                   n.trees = n.trees, interaction.depth = interaction.depth,
                   loss = loss,
                   shrinkage = shrinkage, minsplit = minsplit,
                   bag.fraction = bag.fraction, seed = seed, cp = cp )
    
    #- get the boosting ensemble predictions
    fhat <- predict.gtb( tmpGTB, newdata[, PredNames], n.trees = n.trees )
    
    means.fhat <- rbind( means.fhat, mean(fhat) )
    
    
    #===============================================================================
    # STEP 1b and STEP 2: Update the random effects and variance components
    #===============================================================================
    
    if (modification == "beta") { # in case of Huber and L1 loss, use modified updating 
      
      #----------------------------------------------------
      #- (iii): compute fixed effects
      #----------------------------------------------------

      XTmp <- model.matrix(X, data = dfsTrain) 
      
      #- to compute new betahat:
      betahat_sum1 <- 0
      betahat_sum2 <- 0
      #- let's go:
      for ( ii in 1:NID ) {
        #- get index variable:
        idx <- which( ID == UniqueID[ii] ) 
        #- get person data:
        Yi <- Y[idx]
        fhati <- fhat[idx]
        ei <- Yi - fhati
        Zi <- Z[idx,,drop = FALSE]
        Xi <- XTmp[idx,,drop = FALSE]
        ni <- dim( Zi )[1]
        Ri <- diag( as.numeric( Sigma2hat ), ni )
        Vi <- Zi%*%Dhat%*%t(Zi) + Ri
        InvVi <- solve( Vi )
        #- compute beta:
        betahat_sum1 <- betahat_sum1 + t(Xi)%*%InvVi%*%Xi 
        betahat_sum2 <- betahat_sum2 + t(Xi)%*%InvVi%*%ei
      }
      betahat <- solve( betahat_sum1 ) %*% betahat_sum2
      
      #----------------------------------------------------
      #- (iv): compute random effects and new epsilon
      # and update Dhat and Sigma2hat
      #----------------------------------------------------
      
      DhatNew <- diag( 0, p ) 
      Sigma2hatNew <- 0
      for ( ii in 1:NID ) {
        idx <- which( ID == UniqueID[ii] ) 
        #- get person data:
        Yi <- Y[idx]
        fhati <- fhat[idx]
        ei <- Yi - fhati
        Zi <- Z[idx,,drop = FALSE]
        Xi <- XTmp[idx,,drop = FALSE]
        ni <- dim( Zi )[1]
        Ri <- diag( as.numeric( Sigma2hat ), ni )
        Vi <- Zi%*%Dhat%*%t(Zi) + Ri
        InvVi <- solve( Vi )
        
        #- compute bhat:
        bhat[ii, ] <- Dhat%*%t(Zi)%*%InvVi %*% (ei - Xi %*% betahat)
        
        #- compute new epsilon:
        ehat[idx] <- ei - Zi%*%bhat[ii,]
        
        #- compute new variance components:
        DhatNew <- DhatNew + ( bhat[ii,]%*%t(bhat[ii,]) + ( Dhat - Dhat%*%t(Zi)%*%InvVi%*%Zi%*%Dhat ) )
        tmpSigma <- as.numeric( Sigma2hat )*( ni - as.numeric( Sigma2hat )*sum( diag( InvVi ) ) ) 
        Sigma2hatNew <- Sigma2hatNew + ( t(ehat[idx])%*%ehat[idx] + tmpSigma )  
      }
      
      #- update matrices:
      Dhat <- DhatNew/NID
      Sigma2hat <- Sigma2hatNew/TotalObs  
      
    } else {
      
      #----------------------------------------------------
      #- (iii): compute new bs and new epsilons
      #----------------------------------------------------
      for ( ii in 1:NID ) {
        #- get index variable:
        idx <- which( ID == UniqueID[ii] )
        #- get relevant matrices and vectors:
        Yi <- Y[idx]
        fhati <- fhat[idx]
        Zi <- Z[idx,,drop = FALSE]
        ni <- dim( Zi )[1]
        Ri <- diag( as.numeric( Sigma2hat ), ni )
        Vi <- Zi%*%Dhat%*%t(Zi) + Ri
        InvVi <- solve( Vi )
        #- compute new bhati and new epsilons:
        bhat[ii,] <- Dhat%*%t(Zi)%*%(InvVi%*%( Yi - fhati ))  
        ehat[idx] <- Yi - fhati - Zi%*%bhat[ii,]
      }
      
      #----------------------------------------------------
      #- (iv): update Dhat and Sigma2hat
      #----------------------------------------------------
      DhatNew <- diag( 0, p ) 
      Sigma2hatNew <- 0
      for ( ii in 1:NID ) {
        #- get index variable:
        idx <- which( ID == UniqueID[ii] )
        #- get relevant matrices and vectors:
        Yi <- Y[idx]
        Zi <- Z[idx,,drop = FALSE]
        ni <- dim( Zi )[1]
        Ri <- diag( as.numeric( Sigma2hat ), ni )
        Vi <- Zi%*%Dhat%*%t(Zi) + Ri
        InvVi <- solve( Vi )
        #- compute new variance components:
        DhatNew <- DhatNew + ( bhat[ii,]%*%t(bhat[ii,]) + ( Dhat - Dhat%*%t(Zi)%*%InvVi%*%Zi%*%Dhat ) )
        tmpSigma <- as.numeric( Sigma2hat )*( ni - as.numeric( Sigma2hat )*sum( diag( InvVi ) ) ) 
        Sigma2hatNew <- Sigma2hatNew + ( t(ehat[idx])%*%ehat[idx] + tmpSigma )     
      }
      #- update matrices:
      Dhat <- DhatNew/NID
      Sigma2hat <- Sigma2hatNew/TotalObs  
      
    }
    
    means.ranint <- rbind( means.ranint, mean(bhat[, 1]) )
    errorVarList[[noIterations]] <- Sigma2hat
    DhatList[[noIterations]] <- Dhat
    
    #- Compute GLLnew
    llnew <- mem_boost_gll( Z = Z, ID = ID, bhat = bhat, ehat = ehat, UniqueID = UniqueID, 
                            NID = NID, D = Dhat, Sigma2 = Sigma2hat )
    
    #- Verbose output:
    if( verbose_memboost ) {
      h1 <- paste0("Loglikelihood: ", round( llnew, 2),
                   " | No. iteration: ", noIterations )
      cat(h1, "\n")
      utils::flush.console()    
    }
    #- Leaving the while loop? 
    absDiffLogLik <- abs( (llold - llnew)/llold )
    
    if (  noIterations > minIter_memboost &
          ( absDiffLogLik < conv_memboost | noIterations >= maxIter_memboost )  ) {   
      toIterate <- FALSE
    }      
  } # while
  
  if(absDiffLogLik >= conv_memboost) {
    warning("EM algorithm did not converge")
    convWarning <- TRUE
  }
  
  if( modification != "beta" ){
    betahat <- NULL
  } else {
    betahat <- as.data.frame(betahat)
    colnames(betahat) <- "betahat"
    rownames(betahat)[1] <- "(Intercept)"
  }
  
  
  #- output:
  out <- list( boosting_ensemble = tmpGTB, 
               var_random_effects = Dhat,
               errorVar = Sigma2hat,
               logLik = llnew,
               raneffs = bhat,
               fhat = fhat,
               noIterations = noIterations,
               convWarning = convWarning,
               means.Ystar = means.Ystar, means.fhat = means.fhat,
               means.ranint = means.ranint, DhatList = DhatList,
               errorVarList = errorVarList,
               modification = modification, 
               fixeffs = betahat )
  class(out) <- "MEMBoost"
  return( out )
}