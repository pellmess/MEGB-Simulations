library("saeSim")
library("Pareto")
library("dplyr")
rm(list = ls())
# necessary for prediction in the BHF modell
covariates.mean.func <- function(x) {
  results <- list()
  results[["population"]] <- x$Pop
  results[["sample"]] <- x$Smp
  results[["covariates.mean"]] <- x$Pop %>%
    group_by(idD) %>%
    dplyr::summarise(N.i = n(), across(starts_with('x'), mean)) %>%
    as.data.frame()
  return(results)
}

simruns = 500
Domains = 50
set.seed(123)
sample_size <- round(runif(Domains, 5, 50), digits = 0)
pop_size = rep(1000, Domains)
sum(pop_size)
sum(sample_size)

calc_keepPop = function(dat)
{
  attr(dat, "pop") = dat
  dat
}
#==================================================================================================================================#
#-----------------------------------  Scenario 1: Linear Normal ---------------------------------------------------------------#
#==================================================================================================================================#


gen_x1 = function(dat, m = dat$muD, s = 3) {
  dat["x1"] = rnorm(nrow(dat), mean = m, sd = s)
  return(dat)
}

gen_x2 <- function(dat, m = dat$muD, s = 3) {
  dat["x2"] <-  rnorm(nrow(dat), mean = m, sd = s)
  # dat["x2"] <- rbinom(nrow(dat), 1, 0.8)
  return(dat)
}

gen_myE = function(dat, m = 0, s = 1000) {
  dat["e"] = rnorm(nrow(dat), mean = m, sd = s)
  return(dat)
}

setup = sim_base(data = base_id(nDomains = Domains, nUnits = pop_size)) %>%
  sim_gen(gen_generic(
    runif,
    min = -1,
    max = 1,
    groupVars = "idD",
    name = "muD"
  )) %>%
  sim_gen(gen_x1)       %>%
  sim_gen(gen_x2)       %>%
  as.data.frame %>%
  sim_gen(generator = gen_myE)         %>%
  sim_gen_v(mean = 0, sd = 500)         %>%
  sim_resp_eq(y =  5000 - 500 * x1 - 500 * x2 + v + e)


ts.sim = Sys.time()
# Simulating Population datasets
Pop = sim(setup, R = simruns)

#  Simulating sample datasets

# function for drawing the sample
sampler = function(DAT) {
  smp = as.data.frame(matrix(nrow = sum(sample_size) , ncol = ncol(DAT)))
  brd = append(0, cumsum(sample_size))
  for (i in 1:Domains) {
    # to generate sample with 0 obs. for some domains
    if (sample_size[i] != 0) {
      smp[((brd[i] + 1):brd[i + 1]),] = (DAT[DAT$idD == i, ])[sample(1:sum(DAT$idD == i),
                                                                     size = sample_size[i]), ]
    }
  }
  attr(smp, "pop") = DAT
  colnames(smp) = colnames(DAT)
  return(smp)
}

# creates the sample
Smp <- lapply(Pop, sampler)

# Domain Ids
inDom <- unique(Smp[[2]]$idD)

ch <- numeric()

# Adding additional features to population
for (p in 1:simruns) {
  Pop[[p]]$sam <- NA
  for (i in  1:length(Smp[[p]]$y)) {
    g <-  which((Smp[[p]]$y)[i] == Pop[[p]]$y)
    Pop[[p]]$sam[g] <- "in"
  }
  Pop[[p]]$sam <- ifelse(is.na(Pop[[p]]$sam), "out", "in")
  ch[p] <- table(Pop[[p]]$sam)[1]
}

dom.size.lrg.thn.zro <- sum(sample_size > 0)

for (p in 1:simruns) {
  Pop[[p]]$z1 <- NA
  Pop[[p]]$z1[is.element(Pop[[p]]$idD, inDom)] <-
    rnorm(dom.size.lrg.thn.zro * pop_size[1], 10, 1)
  Pop[[p]]$z1[is.na(Pop[[p]]$z1)] <-
    rnorm((Domains - dom.size.lrg.thn.zro) * pop_size[1], -10, 1)
  for (i in 1:length(Smp[[p]]$y)) {
    Smp[[p]]$z1[i] <- Pop[[p]]$z1[which((Smp[[p]]$y)[i] == Pop[[p]]$y)]
  }
}

sampling.frame <- list()
for (i in 1:simruns) {
  aux.list <- list(Pop = Pop[[i]], Smp = Smp[[i]])
  sampling.frame[[i]] <- aux.list
}

overall.pop.mean <- sapply(
  Pop,
  FUN = function(X) {
    mean(X$y)
  }
)
domain.pop.mean <-
  t(sapply(
    Pop,
    FUN = function(X) {
      tapply(X$y, X$idD, mean)
    }
  ))
domain.pop.var <-
  t(sapply(
    Pop,
    FUN = function(X) {
      tapply(X$y, X$idD, var)
    }
  ))

# for bhf model
input_elements <- lapply(sampling.frame, covariates.mean.func)

te.sim = Sys.time()
time.sim <- te.sim - ts.sim

# saves complete R workspace
save.image("mod1/model_objects.RData")
# saves specified R objects listed below
save(input_elements,
     setup,
     aux.list,
     Domains,
     simruns,
     sample_size,
     time.sim,
     file = "mod1/model_specf_objects.RData")

alt.matrix <-
  data.frame(
    mean = as.vector(t(domain.pop.mean)),
    domain = rep(1:Domains, times = simruns),
    run = rep(1:simruns, each = Domains)
  )
write.csv(alt.matrix, "mod1/results/pop.value.csv", row.names = F)

