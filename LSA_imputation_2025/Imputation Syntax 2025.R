## Francis Huang 2025
## Huang, F., & Keller, B. (2025). Working with missing data in large-scale assessments. Large-scale Assessments in 
## Education. doi: 10.1186/s40536-025-00248-9
## Load the required packages
## Updated 2025.04.18 / using updated version of rblimp which can use factor variables

library(MLMusingR) # for mixPV and summary for pooling
library(WeMix) # for mix
library(tidyr) # to convert from wide to long
library(dplyr) # data management
library(rblimp) #for creating the imputations

## Read in the data
comb <- rio::import("https://github.com/flh3/pubdata/raw/refs/heads/main/miscdata/belgium.rds",
  trust = TRUE) #combined original dataset
wmiss <- rio::import("https://github.com/flh3/pubdata/raw/refs/heads/main/miscdata/belgiumwmiss.rds",
  trust = TRUE) #dataset with additional missing data

## Fit the models with complete and missing data

l1a <- mixPV(pv1math + pv2math + pv3math + pv4math +
  pv5math + pv6math + pv7math + pv8math +
  pv9math + pv10math ~ gender + escs + immig2 +
  stubeha + lackstaff + (1|cntschid),
  weights = c('w_fstuwt', 'w_schgrnrabwt'),
  data = comb, mc = TRUE)
l1b <- mixPV(pv1math + pv2math + pv3math + pv4math +
  pv5math + pv6math + pv7math + pv8math +
  pv9math + pv10math ~ gender + escs + immig2 +
  stubeha + lackstaff + (1|cntschid),
  weights = c('w_fstuwt', 'w_schgrnrabwt'),
  data = wmiss, mc = TRUE)

## Investigate missing data
MLMusingR::nmiss(wmiss) 

## Reshape the data
tall <- pivot_longer(wmiss, pv1math:pv10math, values_to = 'math')
m <- 20 #number of datasets to impute
ns <- nrow(wmiss) #count how many observations there are
tall$.pv <- as.numeric(gsub("[^0-9]", "", tall$name)) #extract the numeric value 
nopv <- length(table(as.character(tall$.pv))) #the number of plausible values

# NOTE: there is no need to do this step anymore with the updated version of rblimp
# tall_numeric <- mutate(tall,
#  across(everything(), as.numeric) #blimp will only work with numeric data
# )  

## Impute the data
mymodel <- rblimp(
  data = tall, #this is using a different dataset as a result of the updated rblimp
  nominal = 'gender immig2 lackstaff',
  # ordinal = '',
  clusterid = 'cntschid',
  fixed = 'gender w_fstuwt', 
  model = 'math ~ gender escs immig2 
           stubeha lackstaff w_fstuwt', 
  options = 'labels',
  seed = 1234,
  nimps = m
) |> by_group('.pv') 


## Check the imputations
mymodel |> sapply(\(x) tail(x@psr, 1) |> max(na.rm = TRUE))
# lapply(mymodel, psr)

## Extract the data
impdat <- lapply(mymodel, as.mitml) |>
  unlist(recursive = FALSE) |>
  do.call(rbind, args = _)

impdat$.implist <- rep(1:(m * nopv), each = ns)

## NO NEED FOR THIS STEP ANYMORE with the updated rblimp
## Return factor attributes to the data
# table(impdat$lackstaff)
# attributes(impdat$lackstaff) <- attributes(comb$lackstaff)
# table(impdat$lackstaff)
# attributes(impdat$gender) <- attributes(comb$gender)
# attributes(impdat$immig2) <- attributes(comb$immig2)

## Create a list for analysis
alldat <- split(impdat, impdat$.implist)

## Analyze the data

## Specify the model
modfit <- function(x) {
  mix(math ~ gender + escs + immig2 + stubeha +
  lackstaff + (1|cntschid),
  weights = c('w_fstuwt', 'w_schgrnrabwt'),
  data = x)
}

#### Method 1: Serial operation (slow)
# allres <- lapply(alldat, FUN = modfit)

#### Method 2: Parallel computation (much faster)
library(parallel)
library(doParallel)
cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl)
clusterExport(cl, list('modfit', 'mix'))
allres <- parLapply(cl, alldat, fun = modfit)
stopCluster(cl)

## Get ready to pool results, need to change the class
class(allres) <- 'mixPV'
summary(allres)

# library(modelsummary)
