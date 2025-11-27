################################################################################
# Function for gradient tree boosting
################################################################################

library(rpart)


gtb <- function( formula, data = NULL,
                 loss = c( "L2", "L1", "Huber" ),
                 n.trees = 100,
                 shrinkage = 0.1,
                 interaction.depth = 1,
                 bag.fraction = 0.5,
                 minsplit = 20, 
                 cp = 0.01,
                 seed = NULL ) 
{ 
  
  if (missing(loss)) {
    stop('A loss function needs to be specified')
  }
  
  if (loss != "L1" & loss != "L2" & loss != "Huber") {
    stop('undefined loss function')
  }
  
  #-----------------------------------------------------------------------------
  # Some preparations
  #-----------------------------------------------------------------------------
  #- get X
  PredNames <- attr( stats::terms( formula ), "term.labels" )
  X <- model.frame(terms(reformulate(PredNames)), data = data)
  
  #- get y
  OutcomeName <- formula[[2]]
  y <- data[ , toString( OutcomeName )]
  
  #- get subsample size
  n <- nrow(X)
  subsample.n <- floor(n * bag.fraction)
  
  #- things for saving
  trees <- vector( "list", n.trees )
  inbag <- matrix(0, ncol = n.trees, nrow = n)
  outcome <-  matrix(NA, ncol = n.trees, nrow = n)
  preds <- matrix( NA, ncol = n.trees, nrow = n )
  leaves <-  matrix(NA, ncol = n.trees, nrow = n)
  
  #- create a vector of seeds
  if (!is.null( seed )) {
    seeds <- seq(from = seed, to = seed + n.trees - 1, by = 1)
  }
  
  #-----------------------------------------------------------------------------
  # Step 1: Initialize f
  #-----------------------------------------------------------------------------
  if (loss == "L2") {
    f <- rep( mean(y), n )
  } else {
    f <- rep( median(y), n )
  }
  
  f0 <- unique(f) # save initializing value
  
  #-----------------------------------------------------------------------------
  # Step 2
  #-----------------------------------------------------------------------------
  # Loop over a total of n.trees iterations
  for (m in 1:n.trees) {
    
    #- a: Get a subsample:
    if (!is.null( seed )) {
      set.seed( seeds[m] )
    }
    inbag.idx <- sample( n, replace = F, size = subsample.n )
    ysub <- y[ inbag.idx ]
    Xsub <- data.frame( X[ inbag.idx, ] )
    colnames(Xsub) <- colnames(X) # to make sure that variables are found when
    # using rpart:::pred.rpart( fit, rpart:::rpart.matrix( X ) )
    inbag[ inbag.idx, m ] <- 1 
    
    #- b: Compute the negative gradient of L at predictions of function f_m-1
    #  for squared error loss, the negative gradient is just the residual
    fsub <- f[ inbag.idx ]
    dsub <- ysub - fsub
    if (loss == "L2") {
      r <- dsub
    } else if (loss == "L1") {
      r <- sign( dsub )
    } else {
      delta <- 1.345 * mad( dsub, constant = 1 ) # as in Lutz et al (2008)
      abs_e <- abs(dsub)
      r <- ifelse( abs_e <= delta, dsub, delta * sign( dsub ) ) 
    }
    outcome[ inbag.idx, m ] <- r
    
    #- c: Fit a regression tree to the target r
    fit <- rpart::rpart( r ~ ., data = Xsub, method = "anova", 
                         control = list( cp = cp, xval = 0, 
                                         maxdepth = interaction.depth,
                                         minsplit = minsplit) )
    
    # save leave indices for full training sample 
    leaves[, m ] <- rpart:::pred.rpart( fit, rpart:::rpart.matrix( X ) )
    
    #- d: Compute new predictions per region
    if (loss == "L2") {
      gamma <- as.vector( predict( fit, newdata = X ) )
      
    } else if (loss == "L1") {
      
      # save median per leaf and replace the leaf predictions in the tree
      # by the median values
      inbagLeaves <- factor( leaves[ inbag.idx, m] )
      leavesPreds <- cbind.data.frame( leaves = inbagLeaves,
                                       median = ave( dsub, inbagLeaves, FUN = median ) )
      # use cbind.data.frame instead of data.frame to stop cbind from
      # converting factor values into numerics
      
      adjPreds <- unique(leavesPreds)
      adjPreds[, 1] <- as.numeric(levels(adjPreds[, 1]))[adjPreds[, 1]]
      # transform leaves (factor) to a numeric using the factor values as numbers
      # (factor values correspond to rows in the tree frame)
      
      # replace the predictions in the fit object
      fit$frame[adjPreds[, 1], ]$yval <- adjPreds[,2]
      
      # now get leaf predictions for the full training sample
      gamma <- as.vector( predict( fit, newdata = X ) )
      leaves[, m ] <- rpart:::pred.rpart( fit, rpart:::rpart.matrix( X ) )
      
    } else {
      
      # compute the Huber loss leaf predictions as in Friedman (2001)
      inbagLeaves <- factor( leaves[ inbag.idx, m] )
      leavesPreds <- cbind.data.frame( leaves = inbagLeaves,
                                       predsTMP1 = ave( dsub, inbagLeaves, FUN = median ),
                                       leaveN = ave( dsub, inbagLeaves, FUN = length ) )
      leavesPreds$predsTMP2 <-  sign( dsub - leavesPreds$predsTMP1 ) *
        pmin( delta, abs( dsub - leavesPreds$predsTMP1 ) )
      leavesPreds <- merge( leavesPreds[, 1:3],
                            aggregate( predsTMP2 ~ leaves, data = leavesPreds, FUN = sum) )
      leavesPreds$preds <- with(leavesPreds, 
                                predsTMP1 + 1/leaveN * predsTMP2 ) 
      
      adjPreds <- unique( leavesPreds[, c("leaves", "preds")])  
      adjPreds[, 1] <- as.numeric(levels(adjPreds[, 1]))[adjPreds[, 1]]
      
      # replace the predictions in the fit object
      fit$frame[adjPreds[, 1], ]$yval <- adjPreds[,2]
      
      # get leaf predictions for the full training sample
      gamma <- as.vector( predict( fit, newdata = X ) )
      leaves[, m ] <- rpart:::pred.rpart( fit, rpart:::rpart.matrix( X ) )
    }
    
    trees[[m]] <- fit
    
    #- e: Update overall prediction
    f <- f + shrinkage * gamma
    
    #- some additional output:
    preds[,m] <- shrinkage * gamma
    
  }
  out <- list( trees = trees, fhat = f, initF = f0,
               preds = preds,
               inbag = inbag, outcome = outcome, leaves = leaves,
               loss = loss, mean.y = mean(y), median.y = median(y),
               shrinkage = shrinkage, formula = formula, n.trees = n.trees )
  class(out) <- "gtb"
  return(out)
} 


#===============================================================================
# function for generating predictions from a gtb or MEMBoost object 
# (In case of MEMboost objects, random effect predictions have to be added
# separately)

predict.gtb <- function(object, X, n.trees){
  
  #- check whether object is a MEM boosting or a standard gtb object
  if( class(object) == "MEMBoost" ){
    modification <- object$modification # extract modification argument
    object <- object$boosting_ensemble  # now use only the boosting object to generate predictions 
  } else {
    modification <- FALSE # for not getting an error in the if statement below when class is "gtb"
  }
  
  if (missing(n.trees)) { # per default use all trees in the ensemble for prediction
    n.trees <- object$n.trees
  }
  
  #- get X 
  PredNames <- attr( stats::terms( object$formula ), "term.labels" )
  X <- model.frame(terms(reformulate(PredNames)), data = X)

  shrinkage <- object$shrinkage
  
  #- get predictions from the ntree base learners and multiply with learning rate
  fit_pred_trees <- sapply( 1:n.trees, function(i) {
    shrinkage * predict(object$trees[[i]],
                               newdata = X,
                               type = "vector") } )
  
  #- sum all the predictions
  # the initializing predictions need to be added separately
  f0 <- object$initF
  if (modification == "beta") { # add the estimated intercept to the predictions
    estIntercept <- object$fixeffs[1, ]
    f0 <- f0 + as.vector( estIntercept )
    } 

  
  if ( is.null( dim(fit_pred_trees) ) ){ 
    fit_pred <- f0 + sum(fit_pred_trees)
  } else { 
    fit_pred <- f0 + apply( fit_pred_trees, 1, sum )
  }
  
  
  return(fit_pred)
}