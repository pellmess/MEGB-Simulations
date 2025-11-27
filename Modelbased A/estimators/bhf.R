for (i in r) {
  no <- i
  input.elements.no <- input_elements[[no]]
  
  
  # BHF-model process -------------------------------------------------------
  
  #Function to obtain all information regarding the eblup.bhf for a simulation
  
  ts.bhf = Sys.time()
  bhf.results <- eblup.bhf.func(input.elements.no, formel)
  te.bhf = Sys.time()
  time.bhf <- te.bhf - ts.bhf
  bhf.results[["time"]] <- time.bhf

  bhf.results.est <- bhf.results$eblup
  names(bhf.results.est)[2] <- "Mean"
  # bhf.results.est$run <- 1
  bhf.results.est$run <- no
  if (file.exists(paste0(mod, "/results/bhf/bhfEstimates.csv"))) {
    write.table(
      bhf.results.est,
      file = paste0(mod, "/results/bhf/bhfEstimates.csv"),
      append = T,
      row.names = F,
      col.names = F
    )
  } else{
    write.table(
      bhf.results.est,
      file = paste0(mod, "/results/bhf/bhfEstimates.csv"),
      append = F,
      row.names = F,
      col.names = T
    )
  }
  
  
}