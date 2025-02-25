### r2 nakagawa
### Francis Huang
### for WeMix
### 2025.02.02

library(WeMix)
library(MLMusingR)

data(pisa2012) #load dataset
pisa2012$one <- 1 #constant weight

##  testing using lmer and mix
# pisa2012$nwt <- pisa2012$w_fstuwt / mean(pisa2012$w_fstuwt)
m1 <- mix(pv1math ~ st04q01 + escs + st29q03 + (1|schoolid),
          weights = c('one', 'one'), cWeights = T,
          data = pisa2012)
# table(predict(m1))
# table(pisa2012$st04q01)
m2 <- lmer(pv1math ~ st04q01 + escs + st29q03 + (1|schoolid),
           data = pisa2012,
           weights = one, REML = F)

####


r2_ns <- function(x, data){
  ## nakagawa and schielzeth r2
  ## only for two level RI models for now
  ## BETA :: BETA :: BETA
  ## 2025.02.02
  
  tmp <- all.vars(formula(x$call))
  n <- length(tmp)
  fml <- paste(tmp[2:(n-1)], collapse = ' + ')
  fml <- paste("~", fml)
  X <- model.matrix(as.formula(fml), data)
  
  ## if there are unused levels
  torem <- which(apply(X, 2, var) == 0)
  if(length(torem) > 1){
    torem <- torem[-1]
    X <- X[,-torem]
  }
  
  ## 
  sigmaf <- var(X %*% matrix(x$coef))
  # sigmaf2 <- var(predict(x))
  t00 = x$vars[1]
  sig2 = x$vars[2]
  
  r2m = (sigmaf) / (sigmaf + t00 + sig2)
  r2c = (sigmaf + t00) / (sigmaf + t00 + sig2)
  
  return(data.frame(r2m = r2m, r2c = r2c))
  
}

performance::r2_nakagawa(m2) #using lmer
r2_ns(m1, pisa2012) #using mix

