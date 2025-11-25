for (i in r) {
  print(i)
  no <- i
  
  path <- paste0(mod, "/results/",
                 "mboost_l2/")
  dir.create(path, showWarnings = F)
  # gradient process ----------------------------------------------------------
  if(file.exists(paste0(path, i, ".txt"))){
    warning(paste("Run no.", i, "has already been completed, skipping"))
    next
  } else{
    write.table(x = i, paste0(path, i, ".txt"))
    input.elements.no <- input_elements[[no]]
    input.elements.no$pop[,group.id] <- as.factor(input.elements.no$pop[,group.id])
    input.elements.no$sample[,group.id] <- as.factor(input.elements.no$sample[,group.id])
    
    # -------------------------------------------------------------------------------
    # - mboost
    # -------------------------------------------------------------------------------
    model <- as.formula(paste(y.names, "~", paste(x.names, collapse = " + ")))
    random <- as.formula(paste("~ 1 |", group.id))
    
    conv_memboost = 0.001
    maxIter_memboost = 100
    verbose_memboost = F
    
    minIter_memboost = 0
    cp = 0.01
    

    
    #-------------------------------------------------------------------------------
    #- Model-based boosting using simple linear models as base-learners
    #-------------------------------------------------------------------------------
    smp_data <- input.elements.no$sample
    pop_data <- input.elements.no$pop
    dfs_mboost <- list( smp_data = smp_data,
                        pop_data = pop_data)
    dfs_mboost <- lapply( dfs_mboost, function(Y){
      # add intercept variable
      Y$int <- rep(1, nrow(Y))
      # mean-center the predictors
      x.names <- 
      Y[, x.names] <- apply(Y[, x.names], 2, function(x) scale(x, center = T, scale = F))
      return(Y)  } )
    
    
    tictoc::tic("MBOOST_L2")
    s_MBOOST_L2 <- Sys.time()
    fitMBOOST_L2 <- mboost::gamboost( y ~ bols(int, intercept = F) + # intercept
                                        bols(x1, intercept = F) + bols(x2, intercept = F) +
                                        brandom(idD), 
                                      data = dfs_mboost$smp_data)
    cvm <- cvrisk( fitMBOOST_L2 ) # find the optimal number of boosting iterations
    fitMBOOST_L2[ mstop(cvm) ] # choose the model based on optimal number of boosting iterations
    e_MBOOST_L2 <- Sys.time()
    time_MBOOST_L2 <- difftime(e_MBOOST_L2, s_MBOOST_L2, units = "secs")[[1]]
    tictoc::toc(log = TRUE)
    

    pop_data_pred <- predict(fitMBOOST_L2, dfs_mboost$pop_data)
    area_MBOOST_L2 <-
      cbind.data.frame(pop_data_pred, idd = input.elements.no$pop[, group.id]) |>
      group_by (idd) |> dplyr::summarise(mw     = mean(pop_data_pred))
    names(area_MBOOST_L2) <- c("domain", "Mean")
    area_MBOOST_L2$run <- no

    
    csv_path <- paste0(path, "/", sprintf("%03d", no), "_mboost_l2.csv")
    message(paste("csv_path:", csv_path))
    
    write.table(
      area_MBOOST_L2,
      file = csv_path,
      row.names = F,
      col.names = T
    )

    time_path <- paste0(path, "/", sprintf("%03d", no), "_time.txt")
    
    write(time_MBOOST_L2, file = time_path, ncolumns = 1)
    
  }
}

