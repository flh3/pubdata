## tidying WeMixResults
## 2024.02.01
## Francis Huang / flhuang2000@yahoo.com

tidy.WeMixResults <- function(m1, ...){
  re.tmp <- m1$varDF
  cfs.res <- unname(c(m1$coef, re.tmp$vcov))
  ses.res <- unname(c(m1$SE, re.tmp$SEvcov))
  crit <- 1.96
  conf.low <- cfs.res - crit * ses.res 
  conf.high <- cfs.res + crit * ses.res 
  zstat <- cfs.res / ses.res
  pv <- 2 * pt(-abs(zstat), Inf)
  
  final <- data.frame(term = c(names(m1$coef), re.tmp$fullGroup),
                      estimate = cfs.res,
                      std.error = ses.res,
                      statistic = zstat,
                      dof = Inf,
                      conf.low = conf.low,
                      conf.high = conf.high,
                      p.value = round(pv, 4))
  row.names(final) <- NULL
  return(final)
}

glance.WeMixResults <- function(m1, ...){
  tmp <- m1$ngroups
  p <- length(tmp)
  ng <- tmp[p]
  lnl <- m1$lnl
  ns <- tmp[1]
  icc <- m1$ICC
  gof <- data.frame(Nobs = ns, lnl = lnl,
                    N.grps = ng, icc.x = icc)
  return(gof)
}



# library(WeMix)
# data(pisa2012, package = 'MLMusingR')
# 
# m0 <- mix(pv1math ~ (1|schoolid), data = pisa2012,
#           weights = c('w_fstuwt', 'w_fschwt'))
# m1 <- mix(pv1math ~ escs + (1|schoolid), data = pisa2012,
#           weights = c('w_fstuwt', 'w_fschwt'))
# m2 <- mix(pv1math ~ escs + st04q01 + (1|schoolid), data = pisa2012,
#           weights = c('w_fstuwt', 'w_fschwt'))
# 
# library(modelsummary)
# modelsummary(list(m0, m1, m2))
