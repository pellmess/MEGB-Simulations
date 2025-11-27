for (i in r) {
  print(i)
  no <- i
  
  path <- paste0(mod, "/results/",
                 "salditt_l2/")
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
    # - Standard boosting
    # -------------------------------------------------------------------------------
    model <- as.formula(paste(y.names, "~", paste(x.names, collapse = " + ")))
    random <- as.formula(paste("~ 1 |", group.id))
    

    n.trees <- 100
    interaction.depth <- 6
    shrinkage <- .1
    bag.fraction <- 0.5
    minsplit <- 20 
    
    tictoc::tic("Standard Boosting L2")
    s_Boost_L2 <- Sys.time()
    fitBoost_L2 <- gtb(model, data = input.elements.no$sample,
                        n.trees = n.trees,
                        loss = "L2",
                        interaction.depth = interaction.depth,
                        shrinkage = shrinkage,
                        bag.fraction = bag.fraction,
                        minsplit = minsplit)
    e_Boost_L2 <- Sys.time()
    tictoc::toc(log = T)
    time_Boost_L2 <- difftime(e_Boost_L2, s_Boost_L2, units = "secs")[[1]]

    pop_data_pred <- predict.gtb( fitBoost_L2, input.elements.no$pop, n.trees = n.trees )
    area_est_L2 <-
      cbind.data.frame(pop_data_pred, idd = input.elements.no$pop[, group.id]) |>
      group_by (idd) |> dplyr::summarise(mw     = mean(pop_data_pred))
    names(area_est_L2) <- c("domain", "Mean")
    area_est_L2$run <- no

    
    csv_path <- paste0(path, "/", sprintf("%03d", no), "_salditt_l2.csv")
    message(paste("csv_path:", csv_path))
    
    write.table(
      area_est_L2,
      file = csv_path,
      row.names = F,
      col.names = T
    )

    time_path <- paste0(path, "/", sprintf("%03d", no), "_time.txt")
    
    write(time_Boost_L2, file = time_path, ncolumns = 1)
    
  }
}

