
# Function to compute the ebp
ebp_func <- function(x, formul_ebp, domains) {
  emdi::ebp(
    fixed = formul_ebp,
    pop_data = x$population,
    smp_data =  x$sample,
    pop_domains = domains,
    smp_domains = domains
  )}


# Function to obtain the eblup estimates only
ebp_mean_func <- function(x) {
  x$ind[,1:2]
}

ebp_sim_results <- function(input.elements, formul, domains) {
  # Computing the ebp estimates
  formul <- as.formula(formul)
  
  ebp_sim <-
    lapply(input.elements, function(x) {
      ebp_func(x, as.formula(formul), domains)
    })
  
  # Calculating the ebp estimator for each simulated population dataset
  ebp_mean <- sapply(ebp_sim, ebp_mean_func)
  ebp_mean_matrix <- matrix(data =  unlist(ebp_mean[2, ]),
                                  nrow = Domains,
                                  byrow = F)

  results = list()
  results[["ebp"]] <- ebp_mean_matrix
  results
}
