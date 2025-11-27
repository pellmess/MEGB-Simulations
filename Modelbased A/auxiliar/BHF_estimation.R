if(!require("sae")) install.packages("sae")
library("sae")

# Function to compute the eblup.bhf 
eblup.bhf.func <- function(x,formul.bhf){
   
  sae::eblupBHF(formul.bhf,data =x$sample,
           dom =idD,meanxpop = x$covariates.mean[,-2],
           popnsize = x$covariates.mean[,c(1,2)] )

}

# Function to obtain the eblup estimates only
eblup.bhf.mean.func <- function(x){
  x$eblup
}

eblup.bhf.sim.results <- function(input.elements,formul){

  # Computing the eblup.bhf estimates
  formel.bhf <- as.formula(formul)
  
  eblup.bhf <-
    lapply(input.elements, function(x) {
      eblup.bhf.func(x, as.formula(formul))
    })
  
  # Calculating the eblup.bhf estimator for each simulated population dataset
  eblup.bhf.mean <- sapply(eblup.bhf, eblup.bhf.mean.func)
  eblup.bhf.mean.matrix <-
    matrix(data =  unlist(eblup.bhf.mean[2, ]),
           nrow = Domains,
           byrow = F)
  eblup.bhf.smpsize.matrix <-
    matrix(data =  unlist(eblup.bhf.mean[3, ]),
           nrow = Domains,
           byrow = F)
  results = list()
  results[["eblup.bhf"]] <- eblup.bhf.mean.matrix
  results[["eblup.bhf.smpsize.matrix"]] <- eblup.bhf.smpsize.matrix
  return(results)
}



