for (i in r) {
  print(i)
  no <- i
  
  
  path <- paste0(mod, "/results/",
                 "xgboost_unit/")
  dir.create(path, showWarnings = F)
  # gradient process ----------------------------------------------------------
  if(file.exists(paste0(path, i, ".txt"))){
    warning(paste("Run no.", i, "has already been completed, skipping"))
    next
  } else{
    write.table(x = i, paste0(path, i, ".txt"))
    input.elements.no <- input_elements[[no]]
    
    
    data_matrix <- dummy_cols(
      input.elements.no$sample[, c(x.names, group.id)],
      select_columns = "idD",
      remove_first_dummy = TRUE,
      remove_selected_columns = TRUE
    ) |> as.matrix()

    labels <- input.elements.no$sample[, y.names]
    
    # Create the DMatrix
    x_data <- xgb.DMatrix(data = data_matrix, label = labels)
    
    gradient_params <- list(
      eta = 0.01,             # learning rate
      nrounds = 10000,        # number of trees
      max_depth = 3,          # maximal depth of trees
      min_child_weight = 3,   # minimal child weight
      lambda = 1,             # L2-regularisation
      alpha = 0,              # L1-regularisation
      gamma = 0.9
    )
    
    
    mod_gb <-
      xgb.cv(
        params = gradient_params[names(gradient_params) != "nrounds"],
        data = x_data,
        verbose = 0,
        nrounds = 10000,
        early_stopping_rounds = 50, print_every_n = 500,
        nthread = 7, nfold = 10, prediction = TRUE, stratified = T
      )
    
    xgb_fit <- xgboost(
      data = x_data,
      verbose = 0,
      params =  gradient_params[names(gradient_params) != "nrounds"],
      nrounds = mod_gb$best_iteration
    )

    pop_data_xgb <- dummy_cols(
      input.elements.no$pop[, c(x.names, group.id)],
      select_columns = "idD",
      remove_first_dummy = TRUE,
      remove_selected_columns = TRUE
    ) %>% as.matrix() |> xgb.DMatrix()
    
    # Combine the boosting estimates and the mixed-effects estimates
    # at unit-level for the population
    unit_preds <- predict(xgb_fit, pop_data_xgb)
    
    # aggregate for the area-level estimations
    unit_preds_ID <- cbind(input.elements.no$pop[, group.id, drop = F], unit_preds)
    
    f0 <- as.formula(paste0("unit_preds ", " ~ ", group.id))
    
    # aggregate estimates by their respective areas
    mean_preds <- aggregate(f0, data = unit_preds_ID,
                            FUN = mean)
    colnames(mean_preds) <- c("domain", "Mean")
    
    mean_preds$run <- no
    
    csv_path <- paste0(path, "/", sprintf("%03d", no), "_gradientEstimates.csv")
    message(paste("csv_path:", csv_path))

    write.table(
      mean_preds,
      file = csv_path,
      row.names = F,
      col.names = T
    )
    
    
    
}
}