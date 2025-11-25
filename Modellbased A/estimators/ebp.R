for (i in r) {
  
  no <- i
  input.elements.no <- input_elements[[no]]
  
  
  # ebp-model process -------------------------------------------------------

  #Function to obtain all information regarding the eblup.bhf for a simulation

  ebp_start = Sys.time()
  ebp_results <- ebp_func(x = input.elements.no,
                          formul_ebp = formel,
                          domains =  group.id)

  ebp_end = Sys.time()
  time_ebp <- ebp_start - ebp_end
  ebp_results[["time"]] <- time_ebp


  ebp_results_est <- ebp_results$ind %>% dplyr::select(domain = Domain, Mean)
  stopifnot(all(dim(ebp_results_est) == c(50, 2)))
  # bhf.results.est$run <- 1
  ebp_results_est$run <- no
  if (file.exists(paste0(mod, "/results/ebp/ebpEstimates.csv"))) {
    write.table(
      ebp_results_est,
      file = paste0(mod, "/results/ebp/ebpEstimates.csv"),
      append = T,
      row.names = F,
      col.names = F
    )
  } else{
    write.table(
      ebp_results_est,
      file = paste0(mod, "/results/ebp/ebpEstimates.csv"),
      append = F,
      row.names = F,
      col.names = T
    )
  }

  
  
  
}
