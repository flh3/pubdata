## wrappers for WeMix to allow for plausible values
## 2023.05.10
## 2023.05.22
## 2023.11.16 added AICbar

mixPV <- function(fml, silent = FALSE, ...){
  res <- list() #empty list
  xx <- deparse1(fml[[2]])
  xx <- gsub(" ", "", xx)
  outc <- unlist(strsplit(xx, split = '\\+'))
  #nout <- length(outc) #number of outcomes
  pred <- fml[[3]]
  
  res <- lapply(outc, function(x){
    if (silent == FALSE) cat("Analyzing plausible value:", x, "\n")
    newf <- reformulate(deparse(pred), response = x)
    mix(newf, ...)
  }
  )
  
  class(res) <- c("mixPV", 'list')
  return(res)
}


tidy.mixPV <- function(out, ...){
  m <- length(out)
  
  re <- sapply(out, \(x) x$vars)
  ns <- nrow(re)
  rese <- sapply(out, \(x) x$varDF$SEvcov)
  cfs <- rbind(sapply(out, coef), re)
  ses <- rbind(sapply(out, \(x) x$SE), rese[1:ns,]) #this is the SE
 
  cfs.res <- rowMeans(cfs)
  if(names(cfs.res)[1] == '') names(cfs.res)[1] <- "(Intercept)"
  B <- apply(cfs, 1, var) #Vb
  ubar <- apply(ses, 1, \(x) mean(x^2)) #Vw
  adj <- (1 + (1/m))
  combvar <- ubar + adj * (B)
  ses.res <- sqrt(combvar)
  tstat <- cfs.res / ses.res
  dof <- (m - 1) * (1 + ((m * ubar) / ((m + 1) * B)))^2
  #from Graham
  RIV <- (B + (B / m)) / ubar #same
  term <- names(cfs.res)
  pv <- 2 * pt(-abs(tstat), dof)
  crit <- abs(qt(p = .025, df = dof)) #.05 / 2
  conf.low <- cfs.res - crit * ses.res 
  conf.high <- cfs.res + crit * ses.res 
  final <- data.frame(term = term,
                      estimate = cfs.res, 
                      std.error = ses.res, 
                      statistic = tstat,
                      dof = dof,
                      conf.low = conf.low,
                      conf.high = conf.high,
                      p.value = round(pv, 4),
                      RIV = RIV)
  row.names(final) <- NULL
  return(final)
}

glance.mixPV <- function(out, ...){
  m <- length(out)
  ns <- summary(out[[1]])$ngroups[1]
  ll <- sapply(out, \(x) x$lnl)
  p <- (nrow(out[[1]]$varDF) + length(coef(out[[1]])))
  AICbar <- mean(-2 * ll) + (2 * p)
  BICbar <- mean(-2 * ll) + (log(ns) * p)
  gof <- data.frame(Nobs = ns, N.pv = m, AICbar = AICbar,
    BICbar = BICbar)
  return(gof)
}

summary_all <- function(x){
  lapply(x, summary)
}

pool_pv <- function(Bs, SEs){
  cfs <- do.call(rbind, Bs)
  ses <- do.call(rbind, SEs)
  m <- nrow(cfs) #number of imputations
  cfs.res <- colMeans(cfs)
  ns <- ncol(cfs) #number of coefficients / estimates
  # does not inc number of covariances for rs models
  B <- apply(data.frame(cfs[,1:ns]), 2, var) #Vb
  ubar <- apply(data.frame(ses[,1:ns]), 2, \(x) mean(x^2)) #Vw
  adj <- (1 + (1/m))
  combvar <- ubar + adj * (B)
  ses.res <- sqrt(combvar)
  tstat <- cfs.res / ses.res
  dof <- (m - 1) * (1 + ((m * ubar) / ((m + 1) * B)))^2
  #from Graham
  RIV <- (B + (B / m)) / ubar #same
  term <- names(cfs.res)
  pv <- 2 * pt(-abs(tstat), dof)
  crit <- abs(qt(p = .025, df = dof)) #.05 / 2
  conf.low <- cfs.res - crit * ses.res 
  conf.high <- cfs.res + crit * ses.res 
  output <- cbind(
    estimate = cfs.res, 
    std.error = ses.res, 
    statistic = tstat,
    df = dof,
    conf.low = conf.low,
    conf.high = conf.high,
    p.value = round(pv, 4),
    RIV = RIV,
    "Pr(>t)" = round(pv, 4)
  )
  return(output)
  
}

summary.mixPV <- function(out){
  require(broom)
  m <- length(out)
  
  #random effects
  re <- lapply(out, \(x) x$vars)
  rese <- lapply(out, \(x) x$varDF$SEvcov)
  re.final <- pool_pv(re, rese)
  
  #fixed effects
  cfs <- lapply(out, coef)
  ses <- lapply(out, \(x) x$SE) #this is the SE
  ttable <- pool_pv(cfs, ses)

  res <- list(ttable = ttable[drop = F],
              reff = re.final,
              gof = glance(out)
              )
  class(res) <- c("mPV", "list")
  return(res)
}

print.mPV <- function(x, ...){
  cat("Results of multilevel analyses with", x$gof$N.pv, "plausible values.\n")
  cat("Number of observations:", x$gof$Nobs, '\n')
  cat("\nEstimates for random effects: \n")
  printCoefmat(x$reff[,-c(5:8)], digits = 3)
  cat("\nEstimates for fixed effects: \n")
  printCoefmat(x$ttable[,-c(5:8), drop = F], digits = 3)
  
}


lrtPV <- function(mf, mr){ #for mixPV
  nll <- sapply(mr, \(x) x$lnl)
  fll <- sapply(mf, \(x) x$lnl)
  a1 <- summary(mr[[1]])
  a2 <- summary(mf[[1]])
  k.r <- length(mr[[1]]$theta) + length(mr[[1]]$coef)
  k.f <- length(mf[[1]]$theta) + length(mf[[1]]$coef)
  k <- k.f - k.r
  lldif <- -2 * (nll - fll)
  m <- length(nll)
  dbar <- mean(lldif)
  r2 <- (1 + (1 / m)) * var(sqrt(lldif))
  Fs <- (dbar/k - ((m + 1)/(m - 1)) * r2) / (1 + r2) #from package
  v <- k^(-3 / m) * (m - 1) * (1 + r2^(-1))^2
  pv <- pf(Fs, k, v, lower.tail = FALSE)
  return(data.frame(F = Fs, df1 = k, df2 = v, r = r2, pv = round(pv, 4)))
}
