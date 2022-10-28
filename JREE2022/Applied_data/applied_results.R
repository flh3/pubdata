## Applied example
## Journal of Research on Educational Effectiveness
## Cluster randomized controlled trial with 18 schools
## AUT: Francis Huang / huangf@missouri.edu
## 2022.06.04 // 2022.07.01

## 1. load function
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

## END pql function

## 2. load package and the dataset
library(dplyr)
library(sjmisc)
library(sjmisc)
library(dplyr)
library(tableone)
library(survey)

dat.rct <- rio::import(file = 'https://github.com/flh3/CR2/raw/master/data/crct.rda')
names(dat.rct)

## 3. Get some descriptives
frq(dat.rct, race_Black, race_Hispanic, female, odr_pre)
dplyr::filter(dat.rct, !duplicated(usid)) %>% 
  with(., range(size))
gsize <- dat.rct %>% group_by(usid) %>% 
  summarise(n = n())
sd(gsize$n) / mean(gsize$n) #coefficient of variation

des <- svydesign(data = select(dat.rct, -starts_with('stype'), -size), ids = ~usid,
                 weights = ~1)
t1 <- svyCreateTableOne(data = des, 
 vars = c('odr_pre', 'female', 
 'race_Black', 'race_Hispanic'),
 strata = 'trt')
print(t1, smd = TRUE, digits = 3)

l2vars <- select(dat.rct, size, trt, stype, usid) %>%
  filter(!duplicated(usid))
l2vars$stype <- rio::factorize(l2vars$stype)
t2 <- CreateTableOne(data = l2vars, vars = c('size', 'stype'),
                     strata = 'trt')
print(t2, smd = TRUE)

## 4. Run the analysis
m2a <- glm(odr_post ~ odr_pre + female + race_Hispanic + 
             trt + size + stype_elem + stype_hs, data = dat.rct,
           family = binomial)
out <- CR2::clustSE(m2a, 'usid', digits = 8) #orig S2
out
# out$CR2

library(nlme)
ml0 <- pql(odr_post ~ 1,
           random = ~1|usid,
           data = dat.rct,
           family = binomial)
performance::icc(ml0)
summary(ml0)
.881^2 / (3.29 + .881^2) #manual ICC computation

## RI; not used
ml1 <- pql(odr_post ~ odr_pre + female + race_Hispanic  +
             trt + size + stype_elem + stype_hs,
           random = ~1|usid,
           data = dat.rct,
           family = binomial)
summary(ml1)

## RS: used
ml2 <- pql(odr_post ~ odr_pre + female + race_Hispanic +
             trt + size + stype_elem + stype_hs,
           random = ~odr_pre|usid, #RS
           data = dat.rct,
           family = binomial)
summary(ml2)

glmm <- summary(ml2)$tTable

# rio::export(glmm, "glmmoutput2.xlsx")
### ENDS HERE 
