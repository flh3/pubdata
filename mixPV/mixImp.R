## Combining imputations in WeMix
## Francis L. Huang / huangf@missouri.edu
## 2024.03.05: testing out pooling

rm(list = ls())

library(dplyr)
library(WeMix)
library(mitml)
library(MLMusingR)

data(pisa2012)
dat <- select(pisa2012, pv1math, escs, 
              gender = st04q01, w1 = sc14q02, 
              w_fstuwt, schoolid)
dat$one <- 1
dat$schoolid <- as.numeric(dat$schoolid)
dat <- dplyr::arrange(dat, schoolid)
xx <- sjmisc::frq(dat$schoolid)
ns <- xx[[1]]$frq #ns

ns2 <- 1:length(ns)

tab <- data.frame(schoolid = xx[[1]]$val, dist = ceiling(seq_along(ns2)/3))

dat <- left_join(dat, tab, by = 'schoolid')
set.seed(123) #remove data at random
dat$escs[sample(1:nrow(dat), 200)] <- NA
dat$gender[sample(1:nrow(dat), 200)] <- NA
mice::md.pattern(dat, rotate.names = TRUE)

# now impute the data!
fml <- list(escs + gender ~ pv1math + (1|schoolid),
            w1 ~ 1)
# just impute quickly; for test purposes ONLY!!
imps <- jomoImpute(data = dat, formula = fml, 
    n.burn = 10, n.iter = 10, m = 5, seed = 123)
comp <- mitmlComplete(imps)

## convert from mids to just a list needed? Not sure

# res <- lapply(comp, \(x) mix(pv1math ~ escs + gender + w1 + (1|schoolid), data = x,
#                       weights = c('w_fstuwt', 'one')))
#FE

summary.mixPV2 <- function(res){
  tmp <- lapply(res, FUN = \(x) summary(x)$coef)
  cfs <- lapply(tmp, FUN = \(x) x[,1])
  ses <- lapply(tmp, FUN = \(x) x[,2])
  cfs1 <- do.call(rbind, cfs)
  
  if(ncol(cfs1) == 1) colnames(cfs1) <- "(Intercept)"
  
  ses1 <- do.call(rbind, ses)
  
  #RE
  re <- lapply(res, FUN = \(x) summary(x)$vars)
  renms <- rownames(re[[1]])
  re1 <- lapply(re, FUN = \(x) x[,1])
  reses <- lapply(re, FUN = \(x) x[,2])
  cfs2 <- do.call(rbind, re1)
  ses2 <- do.call(rbind, reses)
  colnames(cfs2) <- renms
  
  # combined FE and RE
  cfs2pool <- cbind(cfs1, cfs2)
  ses2pool <- cbind(ses1, ses2)
  
  ####
  
  allcfs <- colMeans(cfs2pool) #average of coefficients
  m <- nrow(cfs2pool) #how many imputations
  
  ## Rubin's rules
  ubar <- apply(ses2pool, 2, FUN = \(x) sum(x^2) / m)
  Bb2 <- apply(cfs2pool, 2, var)
  SE.2 <- sqrt(ubar + (1 + (1 / m)) * Bb2)
  res2 <- data.frame(estimate = allcfs, std.error = SE.2, t = allcfs / SE.2)
  # print(nm)
  print(res[[1]]$ngroups)
  print(res2)
}

##

mixImp <- function(fml, mids = NULL, silent = FALSE, ...){
  res <- list() #empty list
  
  if(any(class(mids) %in% 'list') == FALSE) stop("mids must be a list of multiply imputed datasets")
  
  cat("Analyzing", length(mids), "datasets:\n")
  res <- lapply(mids, function(z){
    cat(".") #progress
    mix(fml, ..., data = z)
  })
  cat("Done! \n")
  
  class(res) <- c("mixPV2", 'list')
  return(res)
}

#### run the model and get the results!

res2 <- mixImp(pv1math ~ escs + w1 + (1|schoolid),
       weights = c('w_fstuwt', 'one'), mids = comp)
summary(res2)

res3 <- mixImp(pv1math ~ escs + w1 + (1|dist/schoolid),
 weights = c('w_fstuwt', 'one', 'one'), mids = comp)
summary(res3)

##

r1 <- mix(pv1math ~ escs + w1 + (1|dist/schoolid),
          weights = c('w_fstuwt', 'one', 'one'),
          data = dat, cWeights = T)
summary(r1)

r2 <- mix(pv1math ~ escs + w1 + (1|dist/schoolid),
          weights = c('w_fstuwt', 'one', 'one'),
          data = dat, cWeights = F)
summary(r2)

###
