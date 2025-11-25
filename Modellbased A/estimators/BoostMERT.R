for (i in r) {
  print(i)
  no <- i

  path <- paste0(mod, "/results/",
                 "MERT_l2/")
  dir.create(path, showWarnings = FALSE)
  
  # gradient process ----------------------------------------------------------
  if (file.exists(paste0(path, i, ".txt"))) {
    warning(paste("Run no.", i, "already completed, skipping"))
    next
  } else {
    write.table(x = i, paste0(path, i, ".txt"))
    input.elements.no <- input_elements[[no]]
    input.elements.no$pop[, group.id] <- as.factor(input.elements.no$pop[, group.id])
    input.elements.no$sample[, group.id] <- as.factor(input.elements.no$sample[, group.id])
    
    #-------------------------------------------------------------------------------
    #- MEM Boosting
    #-------------------------------------------------------------------------------
    model <- as.formula(paste(y.names, "~", paste(x.names, collapse = " + ")))
    random <- as.formula(paste("~ 1 |", group.id))
    conv_memboost <- 0.001
    maxIter_memboost <- 100
    verbose_memboost <- FALSE
    
    minIter_memboost <- 0
    cp <- 0.01
    n.trees <- 100
    interaction.depth <- 6
    shrinkage <- 0.1
    bag.fraction <- 0.5
    minsplit <- 20 
    
    tictoc::tic("MERT Boosting L2")
    s_BoostMERT_L2 <- Sys.time()
    fitBoostMERT_L2 <- mem_boost_mert(
      model, data = input.elements.no$sample,
      random = random,
      loss = "L2",
      n.trees = n.trees,
      interaction.depth = interaction.depth,
      shrinkage = shrinkage,
      bag.fraction = bag.fraction,
      minsplit = minsplit,
      maxIter_memboost = maxIter_memboost,
      minIter_memboost = minIter_memboost,
      verbose_memboost = verbose_memboost,
      cp = cp
    )
    e_BoostMERT_L2 <- Sys.time()
    time_BoostMERT_L2 <- difftime(e_BoostMERT_L2, s_BoostMERT_L2, units = "secs")[[1]]
    tictoc::toc()
    
    pop_data <- input.elements.no$pop
    smp_data <- input.elements.no$sample
    
    # Name the random effect IDs
    ids_re <- as.character(unique(smp_data[, group.id]))
    re_vec <- setNames(fitBoostMERT_L2$raneffs, ids_re)
    
    # Random effect per row in pop_data (NA -> 0 if domain not in sample)
    re_pop <- unname(re_vec[as.character(pop_data[, group.id])])
    
    # Prediction = f(x) + u_i
    fhat_Test1 <- predict.gtb(fitBoostMERT_L2$boosting_ensemble, pop_data)
    pop_data_pred <- as.numeric(fhat_Test1) + re_pop
    
    # Aggregate to area-level
    area_est_MERT_L2 <- data.frame(pop_data_pred, idd = pop_data[, group.id]) |>
      dplyr::group_by(idd) |>
      dplyr::summarise(mw = mean(pop_data_pred), .groups = "drop")
    
    names(area_est_MERT_L2) <- c("domain", "Mean")
    area_est_MERT_L2$run <- no
    
    csv_path <- paste0(path, "/", sprintf("%03d", no), "_MERT_l2.csv")
    message(paste("csv_path:", csv_path))
    
    write.table(
      area_est_MERT_L2,
      file = csv_path,
      row.names = FALSE,
      col.names = TRUE
    )
    time_path <- paste0(path, "/", sprintf("%03d", no), "_time.txt")
    
    write(time_BoostMERT_L2, file = time_path, ncolumns = 1)
  }
}
