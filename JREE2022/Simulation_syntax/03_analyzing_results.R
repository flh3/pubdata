### Analyzing results
### simBM binary outcomes
### AUT: Francis Huang
### 2022.06.12

rm(list=ls())

load('sim_results_2000.rda') #processed sim data
dir.create("2022") #just create a subdirectory for results, if not existing
res2 <- sm2

# For tables to be nicely formatted, this has to be installed
# Not from CRAN.
# devtools::install_github(repo="haozhu233/kableExtra", ref="a6af5c0")
# rio::export(res2, file = 'res.xlsx', overwrite = T)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggh4x)
library(kableExtra)
library(sjmisc) #for %nin%

res2$icc <- ifelse(res2$iccsd == .608, "ICC = .10", "ICC = .20")
sum(res2$nonc_glm)
table(res2$icc)
table(res2$iccsd)

## SE BIAS
sum(res2$time)
sum(res2$nonc_glm) #no issues
sum(res2$nonc_glmm)

conv <- dplyr::select(res2, 1:4, icc, nonc_glmm, reps)
conv$nonc <- round((conv$nonc_glmm / conv$reps) * 100, 2)
con2 <- dplyr::select(conv, NG, GS, icc, inter, nonc)

## Table 1 output
(t1 <- reshape2::dcast(con2, GS + NG ~ icc + inter))
#reshape2::dcast(con2, GS + NG ~ inter + icc)
table1 <- t1 %>%
  kableExtra::kable(col.names = c('GS', 'NG', rep(c('PR = 75%', 'PR = 90%'), 2))) %>%
  collapse_rows(c(1,2), valign = "top") %>%
  kable_classic(full_width = F, html_font = "Times") %>%
  add_header_above(c(" " = 2, "ICC = .10" = 2, "ICC = .20" = 2), line = T) 

save_kable(table1, file = '2022/table1.html')
# rio::export(t1, file = '2022/table1.xlsx', overwrite = T)

## Figure 1: Standard error bias
bias <- dplyr::select(res2, 1:4, icc, starts_with('sebias'))
bias1.1 <- pivot_longer(bias, names_to = 'Cond', values_to = 'bias', starts_with('sebias'))
bias1.1$xx <- gsub('sebias.', '', bias1.1$Cond)
bias <- bias1.1 %>% separate(xx, into = c('type', 'coef'))
bias <- mutate(bias,
      type = case_when(
        type == 'cr0' ~ 'CR0',
        type == 'cr2' ~ 'CR2',
        type == 'reml' ~ 'GLMM',
        type == 'ml' ~ 'ML'
      ),
      PR = ifelse(inter == 1.1, 'PR = 75%', 'PR = 90%'),
      size = ifelse(GS == 20, "GS = 20", "GS = 100")
      )
table(bias$size)
bias$size <- relevel(factor(bias$size), ref = 'GS = 20')
sjmisc::frq(bias$type)
#no ml
sebias <- dplyr::filter(bias, type != 'ML') %>%
ggplot(aes(x = NG, y = bias, 
                 group = type)) +
  geom_rect(aes(xmin=10, xmax=50, ymin=-.1, 
                        ymax=.1), 
            fill="grey90", alpha=0.5) +
  geom_hline(yintercept = 0, lwd = .25) +
  geom_line(lty = 'dotted', lwd = .5) +
  geom_point(shape = 21, size = 1.2, alpha = .75, aes(fill = type)) + theme_bw() + 
  facet_nested(size + icc ~ PR + coef) +
  theme(panel.spacing.y = unit(0, "lines")) +
  theme(panel.spacing.x = unit(.25, "lines")) +
  scale_fill_manual(values = c('white','black','grey')) +
  theme(legend.position= 'bottom', legend.direction = 'horizontal') +
  theme(axis.text = element_text(size = 7),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(strip.text = element_text(colour = 'black'),
        strip.background = element_rect(
          color="black", fill="white")) +
  labs(x = 'Number of groups', y = 'Standard error bias',
       fill = '') +
  theme(axis.text = element_text(size = 6)) 
sebias
ggsave('2022/fig1.jpg', sebias, width = 8, height = 5,
       dpi = 300)

### describing SEs ...

subset(res2, NG != 10) %>% dplyr::select(., starts_with('sebias.cr2')) %>%
apply(., 2, range) %>% round(3)
subset(res2, sebias.cr2.b2 < -.10 & NG == 20)
subset(res2, NG != 10) %>% select(., starts_with('sebias.reml')) %>%
  apply(., 2, range)
tmp <- subset(res2, NG != 10) %>% 
  select(., starts_with('sebias.cr2.b2'), 1:4)
tmp
tmp <- subset(res2, NG != 10) %>% 
  select(., starts_with('sebias.reml.b3'), 1:4)

#### Figure 2: Coverage probabilities

cov <- dplyr::select(res2, 1:4, icc, starts_with('cp.'))
tmp2 <- pivot_longer(cov, names_to = 'Cond', values_to = 'coverage', starts_with('cp.'))
tmp2$xx <- gsub('cp.', '', tmp2$Cond)
tmp2$xx <- gsub(' ', '', tmp2$xx) #xtra space?
table(tmp2$xx) 
tmp2 <- separate(tmp2, xx, into = c('type', 'coef'))
sjmisc::frq(tmp2, type)

cov2 <- mutate(tmp2,
               type = case_when(
                 type == 'cr0bm' ~ 'CR0.BM',
                 type == 'cr2bm' ~ 'CR2.BM',
                 type == 'cr0bw' ~ 'CR0.BW',
                 type == 'cr2bw' ~ 'CR2.BW',
                 type == 'reml' ~ 'GLMM',
                 type == 'ml' ~ 'ML'
               ),
               PR = ifelse(inter == 1.1, 'PR = 75%', 'PR = 90%'),
               size = ifelse(GS == 20, "GS = 20", "GS = 100")
)
table(cov2$size)
cov2$size <- relevel(factor(cov2$size), ref = 'GS = 20')
table(cov2$type)  

#####

cpimage <- dplyr::filter(cov2, type != 'ML') %>%
  ggplot(aes(x = NG, y = coverage, 
             shape = type, fill = type)) +
  geom_rect(aes(xmin=10, xmax=50, ymin=92.5, 
                ymax=97.5), 
            fill="grey90", alpha=0.5) +
  geom_hline(yintercept = c(95), lwd = .5) +
  geom_line(lty = 'dotted', lwd = .5) +
  geom_point(size = 1.2) +
  theme_bw() + 
  facet_nested(size + icc ~ PR + coef) + 
  theme(panel.spacing.y = unit(0, "lines")) +
  theme(panel.spacing.x = unit(.25, "lines")) +
  scale_shape_manual(name = '',
                     values = c(21, 24, 21, 24, 25)) +
  scale_fill_manual(name = '',
                    values = c('white','white', 'black', 'black', 'green')) +
  theme(legend.position= 'bottom', legend.direction = 'horizontal') +
  labs(x = "Number of groups", y = 'Coverage', fill = '') +
  theme(strip.background = element_rect(
    color="black", fill="white"),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(strip.text = element_text(colour = 'black')) +
  theme(axis.text = element_text(size = 6)) 
cpimage
ggsave("2022/fig2.jpg", plot = cpimage, 
       width = 8,
       height = 5,
       dpi = 300)

range(res2$cp.cr2bm.b1)
range(res2$cp.cr2bm.b1)
range(res2$cp.cr2bm.b2)
range(res2$cp.cr2bm.b3)
range(res2$cp.cr2bm.b4)

dplyr::select(res2, starts_with("cp.cr2bm.")) %>% range()

#### Tabl 2: Type 1 Error (1 - cp)

cov2$t1 <- round(100 - cov2$coverage, 2)
t1table <- dplyr::filter(cov2, type != 'ML') %>%
  reshape2::dcast(coef + GS + NG ~ icc + PR + type, value.var = 't1')

t1table[,4:22] <- round(t1table[,4:22], 1)
table2 <- t1table[,1:13] %>% #only when ICC = .10
  kableExtra::kable(col.names = c('Coef', 'GS', 'NG', 
  rep(c('CR0.BM', 'CR0.BW', 'CR2.BM', 'CR2.BW', 'GLMM'), 2))) %>%
  collapse_rows(c(1, 2, 3), valign = "top") %>%
  kable_classic(full_width = F, html_font = "Times") %>%
  add_header_above(c(" " = 3, "PR = 75%" = 5, "PR = 90%" = 5), line = T) 

save_kable(table2, file = '2022/table2.html')
# rio::export(t1table, file = '2022/table2.xlsx')

### Figure 3: Power

pwr <- dplyr::select(res2, 1:4, icc, starts_with("power"),
                     -ends_with("b0"), -cv)
pwr <- pivot_longer(pwr, names_to = 'est', values_to = 'power', starts_with('power'))
table(pwr$est)
pwr$xx <- gsub('power.', '', pwr$est)
table(pwr$xx) 
pwr <- separate(pwr, xx, into = c('type', 'coef'))
sjmisc::frq(pwr, type)

pwr <- mutate(pwr,
              type = case_when(
                type == 'cr0bm' ~ 'CR0.BM',
                type == 'cr2bm' ~ 'CR2.BM',
                type == 'cr0bw' ~ 'CR0.BW',
                type == 'cr2bw' ~ 'CR2.BW',
                type == 'reml' ~ 'GLMM',
                type == 'ml' ~ 'ML'
              ),
              PR = ifelse(inter == 1.1, 'PR = 75%', 'PR = 90%'),
              size = ifelse(GS == 20, "GS = 20", "GS = 100")
)

pwr$size <- relevel(factor(pwr$size), ref = 'GS = 20')

fig3 <- filter(pwr, type %nin% c('ML', 'CR0.BW', 'CR0.BM', 'CR2.BW')) %>%
  ggplot(aes(x = NG, y = power, 
             shape = type)) +
  geom_line(lty = 'dotted', lwd = .5) +
  geom_point(size = 1.2) +
  theme_bw() +
  facet_nested(size + icc ~ PR + coef) + 
  theme(panel.spacing.y = unit(0, "lines")) +
  theme(panel.spacing.x = unit(.25, "lines")) +
  scale_shape_manual(name = '',
                     values = c(21, 24, 21, 24, 25)) +
  theme(legend.position= 'bottom', legend.direction = 'horizontal') +
  labs(x = "Number of groups", y = 'Power', fill = '') +
  theme(strip.background = element_rect(
    color="black", fill="white"),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(strip.text = element_text(colour = 'black')) +
  theme(axis.text = element_text(size = 6)) +
  scale_y_continuous(limits = c(0, 1), breaks = c(.25, .5, .75, 1))
fig3

ggsave("2022/fig3.jpg", plot = fig3, 
       width = 7,
       height = 5,
       dpi = 300)

filter(res2, icc == 'ICC = .10', GS == 100, NG == 10) %>%
  dplyr::select(1:5, power.cr2bm.b4, power.reml.b4)

filter(res2, icc == 'ICC = .10', GS == 100, NG == 20) %>%
  dplyr::select(1:5, power.cr2bm.b4, power.reml.b4)

filter(res2, icc == 'ICC = .10', GS == 100, NG == 30) %>%
  select(1:5, power.cr2bm.b1, power.reml.b1)

### Power table

xx <- sm2
xx$icc <- ifelse(xx$iccsd == .608, "ICC = .10", "ICC = .20")
pw <- dplyr::select(xx, 1:4, icc, power.cr2bm.b0, power.cr2bm.b1, power.cr2bm.b2, power.cr2bm.b3, power.cr2bm.b4, 
                    power.reml.b0, power.reml.b1, power.reml.b2, power.reml.b3, power.reml.b4, 
                     -cv) %>% relocate(GS)
tp <- pivot_longer(pw, starts_with('power'))
tp$est <- ifelse(stringr::str_detect( tp$name, 'reml'), "GLMM", "CR2.BM")
tp$cf[stringr::str_detect( tp$name, 'b1')] <- 'b1'
tp$cf[stringr::str_detect( tp$name, 'b2')] <- 'b2'
tp$cf[stringr::str_detect( tp$name, 'b3')] <- 'b3'
tp$cf[stringr::str_detect( tp$name, 'b4')] <- 'b4'
#tp$cf[stringr::str_detect( tp$name, 'b0')] <- 'b0' #for the intercept
tp$value <- round(tp$value * 100)
tp <- dplyr::filter(tp, !is.na(cf)) %>%
  select(-name)

tp2 <- pivot_wider(tp, values_from = value,
            names_from = c(inter, icc, est), 
            id_cols = c(cf, GS, NG))
tp3 <- arrange(tp2, cf, GS, NG)

pcomp <- tp3 %>% kableExtra::kable(col.names = c("Coef", 'GS', 'NG', rep(c('CR2.BM', 'GLMM'), 4))) %>%
  collapse_rows(c(1,2), valign = "top") %>%
  kable_classic(full_width = F, html_font = "Times") %>% 
  add_header_above(c(" " = 3, rep(c("PR = 75%" = 2, "PR = 90%" = 2), 2))) %>%
  add_header_above(c(" " = 3, "ICC = .10" = 4, "ICC = .20" = 4), line = TRUE) 
pcomp

save_kable(pcomp, file = '2022/table3.html')

