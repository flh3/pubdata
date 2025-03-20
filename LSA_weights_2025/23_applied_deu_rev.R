## Applied example showing TIMSS Germany data
## AUT: XXXX
## 2024.06.25
## 2024.10.04

# read in the data


stu <- rio::import("c:/data/timss/deu/asgdeum7.sav")
sch <- rio::import("c:/data/timss/deu/acgdeum7.sav")

library(sjmisc)
library(dplyr)
library(WeMix)
library(jtools)
library(modelsummary)

names(stu) <- tolower(names(stu))
names(sch) <- tolower(names(sch))

# select variables

stu2 <- select(stu, books = asbg04,
               likemath = asbm02e, gender = itsex, totwgt,
               starts_with("asmmat"), idschool, idclass)
sch2 <- select(sch, idschool, econdis =  acbg03a, schwgt, 
               empsucc = acbgeas, resmath = acbgmrs)

# recode
sch2$econdis <- rio::factorize(sch2$econdis)
stu2$gender <- rio::factorize(stu2$gender)
stu2$books <- rio::factorize(stu2$books)
# stu2$likemath <- rio::factorize(stu2$likemath)

# merge
comb <- left_join(stu2, sch2, by = 'idschool')
n_distinct(comb$idschool)
### are weights informative?
mice::md.pattern(comb, rotate.names = T)
tmp <- select(comb, idschool, starts_with("asm"), gender, books,
        econdis, empsucc, totwgt, schwgt) %>%
  tidyr::drop_na()
n_distinct(tmp$idschool)
#mice::md.pattern(tmp, rotate.names = T)
range(tmp$empsucc)
weighted.mean(tmp$empsucc, tmp$schwgt)

psych::describe(select(tmp, starts_with("asm")))
wtest1 <- lm(asmmat01 ~ gender + books + econdis + empsucc,
         data = tmp) #just using one pv
wtest2 <- update(wtest1, asmmat02 ~ .)
wtest3 <- update(wtest1, asmmat03 ~ .)
wtest4 <- update(wtest1, asmmat04 ~ .)
wtest5 <- update(wtest1, asmmat05 ~ .)

weights_tests(wtest1, weights = 'totwgt', data = tmp)
weights_tests(wtest2, weights = 'totwgt', data = tmp)
weights_tests(wtest3, weights = 'totwgt', data = tmp)
weights_tests(wtest4, weights = 'totwgt', data = tmp)
weights_tests(wtest5, weights = 'totwgt', data = tmp)

### DD: <.001  <.001 =.02 .003 .001 /// PS: ns .03 ns ns .04

## NOTE: in this case, the DD test was stat significant but the 
## PS was not. In this case, it is not about using one test over the 
## other

# doing this manually DD

wtest6 <- lm(asmmat01 ~ (gender + books + econdis + empsucc) *
               totwgt, data = tmp)
anova(wtest, wtest6)
# compare F and p value: it is the same as above

# doing this manually PS
wtest7 <- lm(resid(wtest) ~ totwgt, data = tmp)
summary(wtest7)

wtest8 <- lm(I(resid(wtest)^2) ~ totwgt, data = tmp)
summary(wtest8)



####  creating some weights

comb$one <- 1 #weight of one
comb$nwt <- comb$totwgt / mean(comb$totwgt) #normalized weight
comb$pwt1 <- comb$totwgt / comb$schwgt #conditional weight; SAS, Mplus, etc will use this

## cluster weights (mPlus)

# Using base R
# comb$ns <- as.numeric(ave(comb$idschool, comb$idschool, FUN = length)) #how many in cluster (numerator)
# comb$swt <- ave(comb$pwt1, comb$idschool, FUN = sum) #sum of wij (denominator)
# comb$clustw <- with(comb, pwt1 * (ns / swt)) #wij x adjustment


comb <- comb %>% 
  group_by(idschool) %>%
  mutate(ns = n(),
         swt = sum(pwt1),
         clustw = pwt1 * (ns / swt)
        )

## ecluster (effective sample size)
# comb$ess <- sum(comb$pwt1)^2 / sum(comb$pwt1^2)
comb <- comb %>% group_by(idschool) %>%
  mutate(ess = sum(pwt1)^2 / sum(pwt1^2), #effective sample size
         eclust = pwt1 * (ess / sum(pwt1)))


range(comb$clustw) ##both 1
range(comb$eclust)

## with this example of timss, there is no variation within school
comb %>% group_by(idschool) %>%
  summarise(vr = var(pwt1)) %>% filter(vr > 0)


### Descriptives

library(tableone)
CreateTableOne(data = select(tmp, gender, books, asmmat01))

CreateTableOne(data = select(tmp, idschool, econdis, empsucc) %>%
                 filter(!duplicated(idschool)), vars = c('econdis', 'empsucc'))
range(tmp$empsucc)
### NOTE: both are equal to 1-- no variation

#### fit the models

# source("https://raw.githubusercontent.com/flh3/pubdata/main/mixPV/wemix_modelsummary.R")
source("https://raw.githubusercontent.com/flh3/pubdata/refs/heads/main/mixPV/mixPVv2.R")

### unconditional

l2l1 <- mixPV(asmmat01 + asmmat02 + asmmat03 +
                asmmat04 + asmmat05 ~ (1|idschool),
            weights = c('totwgt', 'schwgt'),
            data = comb, mc = TRUE)

nowgt <- mixPV(asmmat01 + asmmat02 + asmmat03 +
               asmmat04 + asmmat05 ~  (1|idschool),
             weights = c('one', 'one'),
             data = comb, cWeights = TRUE, mc = TRUE)

l2 <- mixPV(asmmat01 + asmmat02 + asmmat03 +
            asmmat04 + asmmat05 ~ (1|idschool),
          weights = c('one', 'schwgt'),
          data = comb, cWeights = TRUE, mc = TRUE)

l1 <- mixPV(asmmat01 + asmmat02 + asmmat03 +
            asmmat04 + asmmat05 ~ (1|idschool),
          weights = c('nwt', 'one'),
          data = comb, cWeights = TRUE, mc = TRUE)


## unconditional results
modelsummary(list("two" = l2l1, "L2" = l2, 
                  "L1" = l1, "no" = nowgt), 
             stars = TRUE)

### conditional

l2l1 <- mixPV(asmmat01 + asmmat02 + asmmat03 +
              asmmat04 + asmmat05 ~ gender + books + econdis + empsucc + (1|idschool),
          weights = c('totwgt', 'schwgt'),
          data = comb, mc = TRUE)

nowgt <- mixPV(asmmat01 + asmmat02 + asmmat03 +
               asmmat04 + asmmat05 ~ gender + books + econdis  + empsucc + (1|idschool),
          weights = c('one', 'one'),
          data = comb, cWeights = TRUE, mc = TRUE)

l2 <- mixPV(asmmat01 + asmmat02 + asmmat03 +
            asmmat04 + asmmat05 ~ gender + books + econdis  + empsucc + (1|idschool),
          weights = c('one', 'schwgt'),
          data = comb, cWeights = TRUE, mc = TRUE)

l1 <- mixPV(asmmat01 + asmmat02 + asmmat03 +
            asmmat04 + asmmat05 ~ gender + books + econdis + empsucc + (1|idschool),
          weights = c('nwt', 'one'),
          data = comb, cWeights = TRUE, mc = TRUE)

## to output in modelsummary
# source("https://raw.githubusercontent.com/flh3/pubdata/main/mixPV/wemix_modelsummary.R")

library(modelsummary)
modelsummary(list("two" = l2l1, "L2" = l2, 
                  "L1" = l1, "no" = nowgt), 
             stars = TRUE)

# modelsummary(list("two" = l2l1, "L2" = l2, 
#                   "L1" = l1, "no" = nowgt), 
#              stars = TRUE, out = "appliedres2.docx")


