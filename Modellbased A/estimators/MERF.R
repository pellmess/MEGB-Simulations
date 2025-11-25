merf_files <- list.files("merf_git/R", full.names = T)
lapply(merf_files, source)

for (i in r) {
  print(i)
  no <- i
  input.elements.no <- input_elements[[no]]
  path <- paste0(mod, "/results/",
                 "MERF","/")
  dir.create(path, showWarnings = F)
  # gradient process ----------------------------------------------------------
  if(file.exists(paste0(path, i, ".txt"))){
    warning(paste("Run no.", i, "has already been completed, skipping"))
    next
  } else{
    write.table(x = i, paste0(path, i, ".txt"))

    

    
  no <- i
  input.elements.no <- input_elements[[no]]
  
  
  # MERF Process ------------------------------------------------------------
  


    ts.MERF = Sys.time()
     MERF <- SAEforest_model(
      Y = input.elements.no$sample[, y.names],
      X = input.elements.no$sample[, x.names, drop = F],
      dName = group.id,
      smp_data = input.elements.no$sample,
      pop_data = input.elements.no$population,
      seed = 123, meanOnly = T
    )
    te.MERF = Sys.time()
    time.MERF <- te.MERF - ts.MERF
    MERF.results <- list()
    MERF.results[["MERF"]] <- MERF
    MERF.results[["time"]] <- time.MERF

    save(MERF.results,
         file = paste0(mod, "/results/MERF/MERF", sprintf("%03d", no), ".RData"))

    MERF.results.est <- as.data.frame(MERF.results$MERF$Indicators)
    names(MERF.results.est)[1] <- "domain"
    MERF.results.est$run <- no
    if (file.exists(paste0(mod, "/results/MERF/MERFEstimates.csv"))) {
      write.table(
        MERF.results.est,
        file = paste0(mod, "/results/MERF/MERFEstimates.csv"),
        append = T,
        row.names = F,
        col.names = F
      )
    } else{
      write.table(
        MERF.results.est,
        file = paste0(mod, "/results/MERF/MERFEstimates.csv"),
        append = F,
        row.names = F,
        col.names = T
      )
    }
  } 

  
  
}
