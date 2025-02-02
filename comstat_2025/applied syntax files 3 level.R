# Analyzing results for applied example in
# Huang, F., & Zhang, B. (2025). Cluster-robust standard errors with 
# three-level data. Communications in Statistics - Theory and Methods. 
# doi: 10.1080/03610926.2025.2461609.

library(MLMusingR)
library(modelsummary) #to display output
library(CR2) #contains dataset
library(dplyr) #for data management
library(estimatr)
library(clubSandwich)
library(tidyr)
library(lmerTest)

data(crct)
head(crct)

crct$id <- 1:nrow(crct) #create id per student

n_distinct(crct$usid) #no of schools
#shaping from wide to tall/long
tall <- pivot_longer(crct, c(odr_post, odr_pre))
tall$time <- ifelse(tall$name == 'odr_post', 1, 0)

#fitting an unconditional model to get ICCs
unc <- lmer(value ~ (1|usid/id),
           data = tall)
summary(unc)
performance::icc(unc, by_group = TRUE)

############# ANALYZING RESULTS USING DIFFERENT APPROACHES ::::

m1a <- lm(value ~ stype  + female + race_Black + trt + time + trt * time,
                data = tall)
cr0bm <- coef_test(m1a, vcov = vcovCR(m1a, cluster = tall$usid, type = "CR0"))
cr2bm <- coef_test(m1a, vcov = vcovCR(m1a, cluster = tall$usid, type = "CR2"))
cr0g1 <- coef_test(m1a, vcov = vcovCR(m1a, cluster = tall$usid, type = "CR0"),
          test = 'naive-t')
cr2g1 <- coef_test(m1a, vcov = vcovCR(m1a, cluster = tall$usid, type = "CR2"),
          test = 'naive-t')


m2a <- lmer(value ~ stype  + female + race_Black + trt  + 
          time + trt * time + (1|usid) + (1|usid:id),
          data = tall) #RI
m2b <- lmer(value ~ stype  + female + race_Black + trt  + 
          time + trt * time + (time|usid) + (1|usid:id),
          data = tall) #RS
anova(m2a, m2b) #RS fits better

## need to display output using modelsummary
tidy.coef_test_clubSandwich <- function(xx, conf.int = FALSE, conf.level = 0.95, ...) {
  yy <- data.frame(xx)
  colnames(yy) <- c('term', 'estimate', 'std.error',
                    'statistic', 'xx', 'p.value')
  return(yy)
}

# output of models
modelsummary(list("CR2_BM" = cr2bm, 
                  "CR2_G1" = cr2g1, 
                  "CR0_BM" = cr0bm, 
                  "CR0_G1" = cr0g1,
                  # "MLM-RI" = m2a, #not shown in orig
                  "MLM-RS" = m2b), 
             gof_omit = ".*",
             coef_omit = 'SD|Cor',
             statistic = "({std.error})
                           {p.value}{stars}")

# Results will differ from original analysis as these
# are simulated data and do not include all the covariates
