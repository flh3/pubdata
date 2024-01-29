File to accompany:

https://scholarworks.umass.edu/pare/vol23/iss1/2/

# Using Instrumental Variable Estimation to Evaluate Randomized Experiments with Imperfect Compliance

DOI
https://doi.org/10.7275/k0p6-yj16

Abstract:

Among econometricians, instrumental variable (IV) estimation is a commonly used technique to estimate the causal effect of a particular variable on a specified outcome. However, among applied researchers in the social sciences, IV estimation may not be well understood. Although there are several IV estimation primers from different fields, most manuscripts are not readily accessible by researchers who may only be familiar with regression-based techniques. The manuscript provides a conceptual framework of why and how IV works in the context of evaluating treatment effects using randomized evaluations. I discuss the issue of imperfect treatment compliance, explain the logic of IV estimation, provide a sample dataset, and syntax for conducting IV analysis using R. A goal of the current manuscript is to demystify the use of IV estimation and make evaluation studies that use this technique more readily understood by researchers.

    dat <- rio::import("https://raw.githubusercontent.com/flh3/pubdata/main/IV/ivexample.csv")
    head(dat)

    ##   assign takeup y
    ## 1      0      0 0
    ## 2      0      0 0
    ## 3      0      0 0
    ## 4      0      0 0
    ## 5      0      0 0
    ## 6      0      0 0

    dim(dat)

    ## [1] 200   3

    xtabs(~assign + takeup, data = dat)

    ##       takeup
    ## assign  0  1
    ##      0 91  9
    ##      1 22 78

    library(estimatr)
    m1 <- lm_robust(y ~ assign, data = dat)
    m2 <- iv_robust(y ~ takeup | assign, data = dat)

    ## what if we block those always takers?
    dat$takeup2 <- ifelse(dat$assign == 0 & dat$takeup == 1, 0, dat$takeup)
    xtabs(~assign + takeup2, data = dat)

    ##       takeup2
    ## assign   0   1
    ##      0 100   0
    ##      1  22  78

    m3 <- iv_robust(y ~ takeup2 | assign, data = dat)
    summary(m3)

    ## 
    ## Call:
    ## iv_robust(formula = y ~ takeup2 | assign, data = dat)
    ## 
    ## Standard error type:  HC2 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value  Pr(>|t|) CI Lower CI Upper  DF
    ## (Intercept)     0.90     0.2887   3.118 2.094e-03   0.3307    1.469 198
    ## takeup2         8.91     0.3977  22.406 3.371e-56   8.1260    9.694 198
    ## 
    ## Multiple R-squared:  0.8126 ,    Adjusted R-squared:  0.8117 
    ## F-statistic:   502 on 1 and 198 DF,  p-value: < 2.2e-16

    library(modelsummary)
    modelsummary(list( "ITT" = m1, "LATE" = m2, "TOT" = m3))

<table class="table" style="width: auto !important; margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:left;">
</th>
<th style="text-align:center;">
ITT
</th>
<th style="text-align:center;">
LATE
</th>
<th style="text-align:center;">
TOT
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
(Intercept)
</td>
<td style="text-align:center;">
0.900
</td>
<td style="text-align:center;">
-0.007
</td>
<td style="text-align:center;">
0.900
</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:center;">
(0.289)
</td>
<td style="text-align:center;">
(0.031)
</td>
<td style="text-align:center;">
(0.289)
</td>
</tr>
<tr>
<td style="text-align:left;">
assign
</td>
<td style="text-align:center;">
6.950
</td>
<td style="text-align:center;">
</td>
<td style="text-align:center;">
</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:center;">
(0.519)
</td>
<td style="text-align:center;">
</td>
<td style="text-align:center;">
</td>
</tr>
<tr>
<td style="text-align:left;">
takeup
</td>
<td style="text-align:center;">
</td>
<td style="text-align:center;">
10.072
</td>
<td style="text-align:center;">
</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:center;">
</td>
<td style="text-align:center;">
(0.153)
</td>
<td style="text-align:center;">
</td>
</tr>
<tr>
<td style="text-align:left;">
takeup2
</td>
<td style="text-align:center;">
</td>
<td style="text-align:center;">
</td>
<td style="text-align:center;">
8.910
</td>
</tr>
<tr>
<td style="text-align:left;box-shadow: 0px 1px">
</td>
<td style="text-align:center;box-shadow: 0px 1px">
</td>
<td style="text-align:center;box-shadow: 0px 1px">
</td>
<td style="text-align:center;box-shadow: 0px 1px">
(0.398)
</td>
</tr>
<tr>
<td style="text-align:left;">
Num.Obs.
</td>
<td style="text-align:center;">
200
</td>
<td style="text-align:center;">
200
</td>
<td style="text-align:center;">
200
</td>
</tr>
<tr>
<td style="text-align:left;">
R2
</td>
<td style="text-align:center;">
0.475
</td>
<td style="text-align:center;">
0.978
</td>
<td style="text-align:center;">
0.813
</td>
</tr>
<tr>
<td style="text-align:left;">
R2 Adj.
</td>
<td style="text-align:center;">
0.472
</td>
<td style="text-align:center;">
0.978
</td>
<td style="text-align:center;">
0.812
</td>
</tr>
<tr>
<td style="text-align:left;">
Std.Errors
</td>
<td style="text-align:center;">
HC2
</td>
<td style="text-align:center;">
HC2
</td>
<td style="text-align:center;">
HC2
</td>
</tr>
<tr>
<td style="text-align:left;">
p.value.weakinst
</td>
<td style="text-align:center;">
</td>
<td style="text-align:center;">
</td>
<td style="text-align:center;">
</td>
</tr>
<tr>
<td style="text-align:left;">
statistic.endogeneity
</td>
<td style="text-align:center;">
</td>
<td style="text-align:center;">
</td>
<td style="text-align:center;">
</td>
</tr>
<tr>
<td style="text-align:left;">
p.value.endogeneity
</td>
<td style="text-align:center;">
</td>
<td style="text-align:center;">
</td>
<td style="text-align:center;">
</td>
</tr>
<tr>
<td style="text-align:left;">
statistic.weakinst
</td>
<td style="text-align:center;">
</td>
<td style="text-align:center;">
</td>
<td style="text-align:center;">
</td>
</tr>
<tr>
<td style="text-align:left;">
p.value.overid
</td>
<td style="text-align:center;">
</td>
<td style="text-align:center;">
</td>
<td style="text-align:center;">
</td>
</tr>
<tr>
<td style="text-align:left;">
statistic.overid
</td>
<td style="text-align:center;">
</td>
<td style="text-align:center;">
</td>
<td style="text-align:center;">
</td>
</tr>
</tbody>
</table>
