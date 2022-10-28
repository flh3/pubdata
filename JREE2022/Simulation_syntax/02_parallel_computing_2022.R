### Parallel computing for simulation
### Francis L. Huang
### 2021.10.21
### 2022.06.03 #revisiting


### make sure all required packages are installed

library(nlme)
library(parallel)
library(lme4)
library(CR2) #make sure this is v0.1.1
library(performance)

### Change HD path / working directory
### Helps if all the syntax files are in the same directory
### Load simBM function first (can be done using the source function)
source("01_sim_syntax.R")

pars <- data.frame(reps = 2000, start.c = 1:40, end.c = 1:40, seed = 111)
ll <- split(pars, pars$start.c)

st <- Sys.time()
cl <- makeCluster(detectCores() - 1) #all cores less one
clusterExport(cl, c("simBM", 'glmer', 'lme', 'clustSE', 'MatSqrtInverse'), 
              envir = environment())
out <- parLapply(cl, ll, function(z) {
  simBM(z[,1], z[,2], z[,3], z[,4])
})
Sys.time() - st
stopCluster(cl)
sm2 <- Reduce(rbind, out)

save(sm2, file = 'sim_results_2000.rda')

### On a Dell Precision 5000, takes a little less than 4 hours.

