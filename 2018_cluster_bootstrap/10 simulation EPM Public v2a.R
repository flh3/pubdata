### SIMULATION: MLM vs CBSE
### Author: Francis L. Huang
### 2016


### 0 clear workspace

rm(list = ls())

### 1 load packages

# library(lme4)
library(rms)
library(lmerTest) # new
library(pbkrtest) # new

### 2 generate functions

icc <- function(dat) {
  b <- VarCorr(dat)
  residual_var <- attr(b, "sc")^2
  intercept_var <- attr(b$group, "stddev")[1]^2
  return(icc = intercept_var / (intercept_var + residual_var))
} # computes the ICC



### 3. START simulation

simbs <- function(l2var, reps = 10) {
  # residual variance at level 2
  if (l2var == .05) {
    w2var <- 0.20
  } # .18 originally
  if (l2var == .10) {
    w2var <- 0.45
  } # .41
  if (l2var == .20) {
    w2var <- 0.79
  } # .71
  if (l2var == .30) {
    w2var <- 1.08
  } # .96

  # begin
  options(warn = 0, error = NULL) # checking for error trapping

  L2 <- c(10, 20, 30) # l2 number of clusters
  L1 <- c(10, 30, 50) # l1 group size

  set.seed(1234) # for reproducability

  for (a in 1:length(L2)) {
    collect <- matrix(0, reps, 29) # how much data to collect, 29 columns
    l2size <- L2[a]
    for (b in 1:length(L1)) {
      nG <- eachc <- L1[b] * 1.1 # n per group
      l1size <- l2size * nG
      t1 <- Sys.time() # start time
      # progress bar
      pb <- txtProgressBar(min = 0, max = reps, style = 3)

      for (i in 1:reps) {
        # simulate data
        x <- rnorm(l1size, 0, 1)
        resid <- rnorm(l1size, 0, 1) # level 1 residual variance
        w <- rnorm(l2size, 0, 1) # l2 coefficient
        w <- rep(w, each = eachc)
        wr <- rnorm(l2size, 0, 1)
        wr <- rep(wr, each = eachc) # l2 residual variance
        t <- rep(0:1, each = l1size / 2) # treatment
        # made this fixed, can also use rbinom

        y <- .3 * t + .3 * w + .8 * x + (w2var * wr) + (resid * 1.5) # icc05
        # y<-.3*t+0*w+.3*x+(w2var*wr)+(resid*1.5) #icc05
        group <- gl(l2size, eachc) # grouping variable
        strata <- gl(2, l1size / 2)

        ms <<- unique(as.integer(round(rnorm((l1size - (L1[b] * l2size)) * 2, l1size / 2, l1size / 6), 2)))
        # cat(l1size,L1[b],l2size,length(ms),"\n")
        ms2 <- ms[ms > 0 & ms < (l1size)]
        # print(ms2)
        # cat(l1size,L1[b],l2size,length(ms2))
        ms2 <- sample(ms2, l1size - (L1[b] * l2size))
        # print(length(ms2))

        dat <- data.frame(y, x, w, t, wr, resid, group, strata)
        dat <- dat[-ms2, ] # remove missing random indicators
        ssize <- nrow(dat)
        # get ICC
        m1 <- lmer(y ~ 1 + (1 | group), data = dat)
        ICC <- icc(m1)

        m2 <- lmer(y ~ w + x + t + (1 | group), data = dat)
        mlm1 <- summary(m2)
        vadj <- sqrt(diag(vcovAdj(m2)))
        l2.mlm.pe <- mlm1$coef[2, 1]
        l2.mlm.se <- vadj[2] # mlm1$coef[2,2]
        l1.mlm.pe <- mlm1$coef[3, 1]
        l1.mlm.se <- vadj[3] # mlm1$coef[3,2]
        l2.mlmT.pe <- mlm1$coef[4, 1]
        l2.mlmT.se <- vadj[4] # mlm1$coef[4,2]

        m3 <- lm(y ~ w + x + t, data = dat)
        ol <- summary(m3) # save once to speed up
        l2.ols.pe <- ol$coef[2, 1]
        l2.ols.se <- ol$coef[2, 2]
        l1.ols.pe <- ol$coef[3, 1]
        l1.ols.se <- ol$coef[3, 2]
        l2.olsT.pe <- ol$coef[4, 1]
        l2.olsT.se <- ol$coef[4, 2]

        ## THIS IS USING MAXIMUM LIKELIHOOD
        m4 <- lmer(y ~ w + x + (1 | group), data = dat, REML = F)
        mlm2 <- summary(m4)
        l2.mlmML.pe <- mlm2$coef[2, 1]
        l2.mlmML.se <- mlm2$coef[2, 2]
        l1.mlmML.pe <- mlm2$coef[3, 1]
        l1.mlmML.se <- mlm2$coef[3, 2]

        # bols<-ols(y~w+x,data=dat,x=T,y=T)
        bols <- ols(y ~ w + x + t, data = dat, x = T, y = T)
        c1 <- bootcov(bols, cluster = dat$group, group = dat$strata, B = 1000)
        # c1<-bootcov(bols,cluster=dat$group,B=1000)
        l2.boot.pe <- coef(c1)[2] # should be the same as OLS
        l2.boot.se <- sqrt(diag(vcov(c1)))[2]
        l1.boot.pe <- coef(c1)[3] # should be the same as OLS
        l1.boot.se <- sqrt(diag(vcov(c1)))[3]
        l2.bootT.pe <- coef(c1)[4] # should be the same as OLS
        l2.bootT.se <- sqrt(diag(vcov(c1)))[4]
        # l2.bootT.pe<-coef(c1)[3] #tmp
        # l2.bootT.se<-sqrt(diag(vcov(c1)))[3] #tmp

        successb <- c1$B

        # cat("iteration:",i,"\n")
        setTxtProgressBar(pb, i)

        collect[i, 1] <- l2.mlm.pe
        collect[i, 2] <- l2.mlm.se
        collect[i, 3] <- l1.mlm.pe
        collect[i, 4] <- l1.mlm.se
        collect[i, 5] <- l2.ols.pe
        collect[i, 6] <- l2.ols.se
        collect[i, 7] <- l1.ols.pe
        collect[i, 8] <- l1.ols.se
        collect[i, 9] <- l2.mlmML.pe
        collect[i, 10] <- l2.mlmML.se
        collect[i, 11] <- l1.mlmML.pe
        collect[i, 12] <- l1.mlmML.se
        collect[i, 13] <- l2.boot.pe
        collect[i, 14] <- l2.boot.se
        collect[i, 15] <- l1.boot.pe
        collect[i, 16] <- l1.boot.se
        collect[i, 17] <- ICC
        collect[i, 18] <- l2size
        collect[i, 19] <- l1size
        collect[i, 20] <- m2@optinfo$conv$opt # check for errors
        collect[i, 21] <- m4@optinfo$conv$opt
        collect[i, 22] <- l2.mlmT.pe
        collect[i, 23] <- l2.mlmT.se
        collect[i, 24] <- l2.olsT.pe
        collect[i, 25] <- l2.olsT.se
        collect[i, 26] <- l2.bootT.pe
        collect[i, 27] <- l2.bootT.se
        collect[i, 28] <- successb
        collect[i, 29] <- ssize
      }

      close(pb) # for the progress bar


      t2 <- Sys.time() # end of simulation
      t2 - t1 # time difference between start and end


      colMeans(collect[, c(1, 5, 13)])
      (true.mlm <- sd(collect[, 1]))
      (true.mlmML <- sd(collect[, 9]))
      (true.lm <- sd(collect[, 5]))
      (true.igm <- sd(collect[, 13]))
      (true.mlmT <- sd(collect[, 22]))
      (true.lmT <- sd(collect[, 24]))
      (true.bootT <- sd(collect[, 26]))

      bias.mlm <- (collect[, 2] - true.mlm) / true.mlm
      bias.ols <- (collect[, 6] - true.lm) / true.lm
      bias.mlmML <- (collect[, 10] - true.mlmML) / true.mlmML
      bias.boot <- (collect[, 14] - true.igm) / true.igm
      bias.mlmT <- (collect[, 23] - true.mlmT) / true.mlmT
      bias.olsT <- (collect[, 25] - true.lmT) / true.lmT
      bias.bootT <- (collect[, 27] - true.bootT) / true.bootT


      deff <- 1 + (collect[, 17] * ((l1size / l2size) - 1))
      bias.deff <- ((collect[, 6] * sqrt(deff)) - true.lm) / true.lm


      if (a == 1 & b == 1) {
        nn <- data.frame(
          n = l1size, G = l2size,
          ols = round(mean(bias.ols), 3),
          deff = round(mean(bias.deff), 3),
          mlmML = round(mean(bias.mlmML), 3),
          mlm = round(mean(bias.mlm), 3),
          boot = round(mean(bias.boot), 3),
          mlmT = round(mean(bias.mlmT), 3),
          olsT = round(mean(bias.olsT), 3),
          bootT = round(mean(bias.bootT), 3),
          icc = round(mean(collect[, 17]), 3), time = t2 - t1
        )
        results <- collect
      } else {
        nn <- rbind(nn, data.frame(
          n = l1size, G = l2size,
          ols = round(mean(bias.ols), 3),
          deff = round(mean(bias.deff), 3),
          mlmML = round(mean(bias.mlmML), 3),
          mlm = round(mean(bias.mlm), 3),
          boot = round(mean(bias.boot), 3),
          mlmT = round(mean(bias.mlmT), 3),
          olsT = round(mean(bias.olsT), 3),
          bootT = round(mean(bias.bootT), 3),
          icc = round(mean(collect[, 17]), 3), time = t2 - t1
        ))
        results <- rbind(results, collect)
      }
    }
  }

  colnames(results) <- c(
    "l2.mlm.pe",
    "l2.mlm.se",
    "l1.mlm.pe",
    "l1.mlm.se",
    "l2.ols.pe",
    "l2.ols.se",
    "l1.ols.pe",
    "l1.ols.se",
    "l2.mlmML.pe",
    "l2.mlmML.se",
    "l1.mlmML.pe",
    "l1.mlmML.se",
    "l2.boot.pe",
    "l2.boot.se",
    "l1.boot.pe",
    "l1.boot.se",
    "ICC",
    "l2size",
    "l1size",
    "errorREML",
    "errorML",
    "l2.mlmT.pe",
    "l2.mlmT.se",
    "l2.olsT.pe",
    "l2.olsT.se",
    "l2.bootT.pe",
    "l2.bootT.se",
    "Breps", "ssize"
  )

  Icc <- round(mean(results[, 17]), 2)
  ff <- paste0("collect_icc.", Icc, ".csv")
  fr <- paste0("summary_icc.", Icc, ".rdata")
  fz <- paste0("collect_icc.", Icc, ".rdata")
  write.csv(results, file = ff, row.names = F)
  save(nn, file = fr)
  save(results, file = fz)
}

# 4 Run simulation
simbs(l2var = .05, reps = 1000)
simbs(l2var = .10, reps = 1000)
simbs(l2var = .20, reps = 1000)
simbs(l2var = .30, reps = 1000)
