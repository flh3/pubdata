## Monte-Carlo simulation for JREE paper
## Using Robust Standard Errors for the Analysis of Binary Outcomes with a Few Clusters
## AUT: Francis Huang
## 2022.06.11 {Revised / tested}

simBM <- function(reps, start.c, end.c, seed){
  
coverage <- function(b, se, true, level = .95, dof = dof){ # Estimate,
  qtile <- level + (1 - level) / 2 # Compute the proper quantile
  lower.bound <- b - qt(qtile, df = dof) * se # Lower bound
  upper.bound <- b + qt(qtile, df = dof) * se # Upper bound
  true.in.ci <- ifelse(true >= lower.bound & true <= upper.bound, 1, 0)
  cp <- mean(true.in.ci, na.rm = T) # The coverage probability
  return(coverage.probability = cp * 100)
} # coverage rate


## modified pql function (from MASS package together with NLME) using REML
## Venables, W. N. & Ripley, B. D. (2002) Modern Applied
## Statistics with S. Fourth Edition. Springer, New York. ISBN
##  0-387-95457-0
##
## Pinheiro J, Bates D, DebRoy S, Sarkar D, R Core Team (2021).
## _nlme: Linear and Nonlinear Mixed Effects Models_. R package
## version 3.1-153, <URL: https://CRAN.R-project.org/package=nlme>.

pql <- function (fixed, random, family, data, correlation, weights, 
                 control, niter = 10, verbose = TRUE, ...) 
{
  if (!requireNamespace("nlme", quietly = TRUE)) 
    stop("package 'nlme' is essential")
  if (is.character(family)) 
    family <- get(family)
  if (is.function(family)) 
    family <- family()
  if (is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }
  m <- mcall <- Call <- match.call()
  nm <- names(m)[-1L]
  keep <- is.element(nm, c("weights", "data", "subset", 
                           "na.action"))
  for (i in nm[!keep]) m[[i]] <- NULL
  allvars <- if (is.list(random)) 
    allvars <- c(all.vars(fixed), names(random), unlist(lapply(random, 
                                                               function(x) all.vars(formula(x)))))
  else c(all.vars(fixed), all.vars(random))
  Terms <- if (missing(data)) 
    terms(fixed)
  else terms(fixed, data = data)
  off <- attr(Terms, "offset")
  if (length(off <- attr(Terms, "offset"))) 
    allvars <- c(allvars, as.character(attr(Terms, "variables"))[off + 
                                                                   1])
  if (!missing(correlation) && !is.null(attr(correlation, "formula"))) 
    allvars <- c(allvars, all.vars(attr(correlation, "formula")))
  Call$fixed <- eval(fixed)
  Call$random <- eval(random)
  m$formula <- as.formula(paste("~", paste(allvars, collapse = "+")))
  environment(m$formula) <- environment(fixed)
  m$drop.unused.levels <- TRUE
  m[[1L]] <- quote(stats::model.frame)
  mf <- eval.parent(m)
  off <- model.offset(mf)
  if (is.null(off)) 
    off <- 0
  wts <- model.weights(mf)
  if (is.null(wts)) 
    wts <- rep(1, nrow(mf))
  mf$wts <- wts
  fit0 <- glm(formula = fixed, family = family, data = mf, 
              weights = wts, ...)
  w <- fit0$prior.weights
  eta <- fit0$linear.predictors
  zz <- eta + fit0$residuals - off
  wz <- fit0$weights
  fam <- family
  nm <- names(mcall)[-1L]
  keep <- is.element(nm, c("fixed", "random", "data", 
                           "subset", "na.action", "control"))
  for (i in nm[!keep]) mcall[[i]] <- NULL
  fixed[[2L]] <- quote(zz)
  mcall[["fixed"]] <- fixed
  mcall[[1L]] <- quote(nlme::lme)
  mcall$random <- random
  mcall$method <- "REML" #just changed this
  if (!missing(correlation)) 
    mcall$correlation <- correlation
  mcall$weights <- quote(nlme::varFixed(~invwt))
  mf$zz <- zz
  mf$invwt <- 1/wz
  mcall$data <- mf
  for (i in seq_len(niter)) {
    if (verbose) 
      message(gettextf("iteration %d", i), domain = NA)
    fit <- eval(mcall)
    etaold <- eta
    eta <- fitted(fit) + off
    if (sum((eta - etaold)^2) < 1e-06 * sum(eta^2)) 
      break
    mu <- fam$linkinv(eta)
    mu.eta.val <- fam$mu.eta(eta)
    mf$zz <- eta + (fit0$y - mu)/mu.eta.val - off
    wz <- w * mu.eta.val^2/fam$variance(mu)
    mf$invwt <- 1/wz
    mcall$data <- mf
  }
  attributes(fit$logLik) <- NULL
  fit$call <- Call
  fit$family <- family
  fit$logLik <- as.numeric(NA)
  oldClass(fit) <- c("glmmPQL", oldClass(fit))
  fit
}

library(CR2)
library(MASS)
library(lme4)

NG <- c(10, 20, 30, 40, 50)
GS <- c(20, 100)
cv <- c(.25) #appendix c; cv = 0 || appendix d; cv = .50
inter <- c(1.1, 2.2)
iccsd <- c(.608, .917)
beta <- .3

cond <- expand.grid(NG = NG, GS = GS, cv = cv, inter = inter,
                    iccsd = iccsd, beta = beta)[start.c:end.c, ]
nconds <- nrow(cond)
results <- list()
reps <- reps #rev

cfs1 <- cfs2 <- cfs3 <- cfs4 <- cfs5 <- 
  ses1 <- ses2 <- ses3 <- ses4 <- ses5 <-
  pvs1 <-
  pvs2 <- pvs3 <- pvs4 <- pvs5 <- pvs6 <-
  ses.cr0 <- pvs.cr0.bm <- pvs.cr0.bw <- 
  dofbm <- matrix(NA, ncol = 5, nrow = reps)
l2icc <- t00 <- numeric(reps)


for (j in 1:nconds){
  
  NG <- cond[j, 'NG']
  GS <- cond[j, 'GS']
  cv <- cond[j, 'cv']
  inter <- cond[j, 'inter']
  iccsd <- cond[j, 'iccsd']
  beta <- cond[j, 'beta']
  sd1 <- GS * cv
 
set.seed(seed)
st <- Sys.time()

for (i in 1:reps){

## GENERATE DATA :::::::::::
  
  tr <- c(rep(0, NG / 2), rep(1, NG / 2))
  imbal <- round(rnorm(NG, mean = GS, sd = sd1), digits = 0)
  imbal[imbal < 5] <- 5
  total <- sum(imbal)

  tr <- rep(tr, imbal)
  w1 <- rep(rnorm(NG), imbal)
  e2 <- rep(rnorm(NG, 0, iccsd), imbal) #.608
  school <- rep(1:NG, imbal)

  t11 <- sqrt((iccsd^2) / 5)
 
  e2rs = 0 #no rs
  # e2rs <- rep(rnorm(NG, 0, t11), imbal) ### to include RS
  
  x1 <- rnorm(total)

  ystar <- inter + tr * beta + w1 * beta + x1 * .3 + w1 * x1 * beta + e2 + e2rs * x1
  y <- rbinom(total, 1, plogis(ystar)) #same as pp
  dat <<- data.frame(y, school, tr, w1, x1)

## ANALYZE DATA :::::::::::::
  
  nm <- tryCatch(glmmPQL(y ~ 1, random = ~1|school,
                data = dat, family = binomial, verbose = F),
                error = function(e){
                 print('nullm') 
                }
  )
  
  r1 <- tryCatch(glm(y ~ tr + w1 + x1 + w1 * x1, 
                     data = dat, family = binomial),
                 warning = function(w){
                   print('glm_conv')
                 }
  )

  ### ML
  r2 <- tryCatch(glmmPQL(y ~ tr + w1 + x1 + w1 * x1, 
                         random = ~1|school,
                         data = dat, family = binomial, verbose = F),
                 error = function(e){
                   print('glmm_conv')
                 }
  )

    
  ## REML
  r3 <- tryCatch(pql(y ~ tr + w1 + x1 + w1 * x1, data = dat, 
                     random = ~1|school,
            family = binomial, verbose = F),
    error = function(e){
    cat(i, 'pql_conv', '\n')
    }
  )
 
  if (length(class(r1)) == 2){ #glm lm
    tmp <- clustSE(r1, 'school')
    dofbm[i, ] <- tmp$dfBM
    cfs1[i,] <- coef(r1)
    ses1[i, ] <- tmp$CR2
    pvs1[i, ] <- tmp$CR2pv
    pvs4[i, ] <- tmp$CR2pv.n
    
    ses.cr0[i, ] <- tmp$CR0
    pvs.cr0.bm[i, ] <- tmp$CR0pv
    pvs.cr0.bw[i, ] <- tmp$CR0pv.n
   }
  
    if(length(class(r2)) == 2){ #glmmPQL
      
      tmp2 <- summary(r2)$tTable
      
      if (sum(tmp2[,2] < 20)){ #catch if SE is too big!!
        cfs2[i,] <- tmp2[,1]
        ses2[i, ] <- tmp2[,2]
        pvs2[i, ] <- tmp2[,5]
      }
    }
    
  ### GLMM
  
  if(length(class(nm)) == 2){
    l2icc[i] <- as.numeric(performance::icc(nm)[1])
  }
  
  if(length(class(r3)) == 2){
    t00[i] <- as.numeric(VarCorr(r3)[1,1]) #conditional variance
    tmp3 <- summary(r3)$tTable
    cfs3[i,] <- tmp3[,1]
    ses3[i, ] <- tmp3[,2]
    pvs3[i, ] <- tmp3[,5]
  }  
}

## END LOOP SIM
## BEGIN ANALYSIS
  
probl <- which(ses2[,2] > 5) #invalid SEs
ses2[probl,] <- cfs2[probl,] <- pvs2[probl, ] <- t00[probl] <- l2icc[probl] <- NA
ses3[probl,] <- cfs3[probl,] <- pvs3[probl, ] <- NA

probl2 <- which(ses3[,1] > 2 | ses3[,2] > 2)
cfs3[probl2, ] <- ses3[probl2, ] <- pvs3[probl2, ] <- 
  t00[probl2] <- l2icc[probl2] <- NA

nonc_glm <- sum(is.na(cfs1[,1]))
nonc_glmm <- sum(is.na(cfs3[,1])) #REML

setv1 <- apply(cfs1, 2, sd, na.rm = T) #GLM
setv2 <- apply(cfs2, 2, sd, na.rm = T) #ML
setv3 <- apply(cfs3, 2, sd, na.rm = T) #REML
# setv4 <- apply(cfs4, 2, sd, na.rm = T)
setv5 <- apply(cfs5, 2, sd, na.rm = T) #glmer

iccm <- mean(l2icc, na.rm =T)
t00[t00 > 5] <- NA #false convergence / problems
t00m <- mean(t00, na.rm =T) #not ICC, variance

sebias.cr2 <- (apply(ses1, 2, mean, na.rm = T) / setv1) - 1 #CR2
sebias.ml <- (apply(ses2, 2, mean, na.rm = T) / setv2) - 1 #ML
sebias.reml <- (apply(ses3, 2, mean, na.rm = T) / setv3) - 1 #REML
# (apply(ses4, 2, mean, na.rm = T) / setv4) - 1 #GEE
sebias.glmer <- (apply(ses5, 2, mean, na.rm = T) / setv5) - 1 #glmer
sebias.cr0 <- (apply(ses.cr0, 2, mean, na.rm = T) / setv1) - 1 #CR0

power.cr2bm <- colMeans(pvs1 <= .05, na.rm = T) #CR2 BM
power.ml <- colMeans(pvs2 <= .05, na.rm = T) #ML
power.reml <- colMeans(pvs3 <= .05, na.rm = T) #REML
power.cr2bw <- colMeans(pvs4 <= .05, na.rm = T) #BW
# colMeans(pvs5 <= .05, na.rm = T) #GEE
power.glmer <- colMeans(pvs6 <= .05, na.rm = T) #glmer
power.cr0bm <- colMeans(pvs.cr0.bm <= .05, na.rm = T) #CR0 BM
power.cr0bw <- colMeans(pvs.cr0.bw <= .05, na.rm = T) #CR0 BW

colMeans(cfs1, na.rm = T) #GLM
colMeans(cfs2, na.rm = T) #PQL ML
colMeans(cfs3, na.rm = T) #PQL REML
# colMeans(cfs4, na.rm = T) #GEE
colMeans(cfs5, na.rm = T) #glmer

bm <- colMeans(dofbm, na.rm = T)
dfn <- tmp$dfn

adj <- sqrt(3.29 / (3.29 + iccsd^2)) 
tv.cs <- c(inter, beta, .3, .3, .3) #glmm #1.1 and 2.2
tv.pa <- tv.cs * adj #conversion done here

cp.cr2.bw <- sapply(1:5, function(x) #CR2 w BW
  coverage(cfs1[,x], ses1[,x], true = tv.pa[x], dof = dfn[x])
)

cp.cr0.bw <- sapply(1:5, function(x) #CR0 w BW
  coverage(cfs1[,x], ses.cr0[,x], true = tv.pa[x], dof = dfn[x])
)

cp.ml <- sapply(1:5, function(x) #ML
  coverage(cfs2[,x], ses2[,x], true = tv.cs[x], dof = dfn[x])
)

cp.reml <- sapply(1:5, function(x) #ML
  coverage(cfs3[,x], ses3[,x], true = tv.cs[x], dof = dfn[x])
)

cp.cr2.bm <- sapply(1:5, function(x) #CR2 w BM
  coverage(cfs1[,x], ses1[,x], true = tv.pa[x], dof = bm[x])
)

cp.cr0.bm <- sapply(1:5, function(x) #CR0 w BM
  coverage(cfs1[,x], ses.cr0[,x], true = tv.pa[x], dof = bm[x])
)


tm <- as.numeric(difftime(Sys.time(), st, units = 'mins'))

results[[j]] <- list(
  NG, GS, cv, inter, beta, iccsd, iccm, nonc_glm, nonc_glmm, i, #9
  sebias.cr2, sebias.ml, sebias.reml, #sebias.glmer, #20
  sebias.cr0,
  power.cr2bm, power.cr2bw, power.ml, power.reml, #power.glmer, #25
  power.cr0bm, power.cr0bw,
  cp.cr2.bm, cp.cr2.bw, cp.ml, cp.reml, #cp.lapl #25
  cp.cr0.bm, cp.cr0.bw, tm, t00m
)

}

res2 <- data.frame(matrix(unlist(results), ncol = 92,
                          byrow = T))

names(res2) <- c(
  'NG', 'GS', 'cv', 'inter', 'beta', 'iccsd', 'iccm', 'nonc_glm', 
  'nonc_glmm', 'reps',
  paste0('sebias.cr2.b', 0:4), 
  paste0('sebias.ml.b', 0:4), 
  paste0('sebias.reml.b', 0:4), 
  #paste0('sebias.glmer.b', 0:4),
  paste0('sebias.cr0.b', 0:4),
  paste0('power.cr2bm.b', 0:4), paste0('power.cr2bw.b', 0:4), 
  paste0('power.ml.b', 0:4), 
  paste0('power.reml.b', 0:4), #paste0('power.glmer.b', 0:4),
  paste0('power.cr0bm.b', 0:4),
  paste0('power.cr0bw.b', 0:4),
  paste0('cp.cr2bm.b', 0:4), paste0('cp.cr2bw.b', 0:4), 
  paste0('cp.ml.b', 0:4), 
  paste0('cp.reml.b', 0:4), 
  #paste0('cp.glmer.b', 0:4)
  paste0('cp.cr0bm.b', 0:4), paste0('cp.cr0bw.b', 0:4),
  'time', 'T00'
  )

  return(res2)
}
