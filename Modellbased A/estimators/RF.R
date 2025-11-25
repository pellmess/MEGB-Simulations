for (i in r) {
  print(i)
  no <- i
  
  path <- paste0(mod, "/results/",
                 "rf/")
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
    # - RF
    # -------------------------------------------------------------------------------
    model <- as.formula(paste(y.names, "~", paste(x.names, collapse = " + ")))
    random <- as.formula(paste("~ 1 |", group.id))
    

    tictoc::tic("Standard RF")
    s_RF <- Sys.time()
    fitRF <- randomForest::randomForest(formula = model, data = input.elements.no$sample)
    e_RF <- Sys.time()
    tictoc::toc(log = T)
    time_RF <- difftime(e_RF, s_RF, units = "secs")[[1]]
    
    pop_data_pred <- predict(fitRF, input.elements.no$pop)
    area_est_RF <- cbind.data.frame(pop_data_pred, idd = input.elements.no$pop[,group.id]) |>
      group_by(idd) |>
      dplyr     ::summarise(mw = mean(pop_data_pred))

    names(area_est_RF) <- c("domain", "Mean")
    area_est_RF$run <- no

    
    csv_path <- paste0(path, "/", sprintf("%03d", no), "_rf.csv")
    message(paste("csv_path:", csv_path))
    
    write.table(
      area_est_RF,
      file = csv_path,
      row.names = F,
      col.names = T
    )
    
    time_path <- paste0(path, "/", sprintf("%03d", no), "_time.txt")
    
    write(time_RF, file = time_path, ncolumns = 1)
    
    
  }
}

