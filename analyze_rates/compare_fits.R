library(tidyverse)

model.fits <- read.csv("model_fits.csv")

model.fits %>% 
    group_by(dataset) %>%
    arrange(dataset, AICc) %>%
    mutate(AICc_rank = rank(AICc)) %>% 
    mutate(deltaAICc = min(AICc) - AICc) -> model.aicc
model.aicc

model.aicc %>% 
    filter(AICc_rank == 1) %>% 
    group_by(model, type) %>%
    tally() %>%
    arrange(type, desc(n))
# A tibble: 9 x 3
# Groups:   model [4]
#      model   type     n
#     <fctr> <fctr> <int>
# 1 gcpREV-G chloro    56
# 2    JTT-G chloro     2
# 3  mtVer-G chloro     1
# 4    WAG-G chloro     1
# 5    WAG-G enzyme    95
# 6 gcpREV-G enzyme     3
# 7    JTT-G enzyme     2
# 8  mtVer-G   mito     9
# 9 gcpREV-G   mito     4


model.aicc %>% 
    filter(AICc_rank == 2) %>% 
    ungroup() %>%
    summarize(mean(deltaAICc), median(deltaAICc), min(abs(deltaAICc)), max(abs(deltaAICc)))
#   `mean(deltaAICc)` `median(deltaAICc)` `min(abs(deltaAICc))`
#               <dbl>               <dbl>                 <dbl>
#1         -1145.358           -735.0449              2.446922
# # ... with 1 more variables: `max(abs(deltaAICc))` <dbl>

model.aicc %>%
    filter(AICc_rank == 2) %>% 
    ungroup() %>%
    filter(abs(deltaAICc) <= 10) %>%
    select(dataset, type, model, deltaAICc) %>%
    arrange(dataset, type)
#   dataset   type     model deltaAICc
#    <fctr> <fctr>    <fctr>     <dbl>      top AICc model
# 1    petN chloro gcpREV-No -2.446922  ---> gcpREV-G 
# 2    psbF chloro gcpREV-No -5.801242  ---> gcpREV-G
# 3   rps12 chloro  gcpREV-G -6.366939  ---> JTT-G
# 

################### Fit difference between gamma and no rv across all models ####################


model.fits %>% 
    separate(model, c("model", "rv"), "-") %>% 
    select(-AICc, -GammaShape) %>% 
    group_by(dataset, type, model) %>%  
    spread(rv, logl) %>% 
    mutate(LRT = 2*(G - No), 
           p.raw = (1 - pchisq(LRT, 1)),
           p.correct = ifelse((p.raw * 6 >= 1), 1, p.raw*6)) -> rv.lrt  ### 6 tests per dataset


## Only 4 datasets, all choloroplast, show any evidence that Gamma is not an improvement
## petN:  4 models
## psbF:  3 models
## psbL:  1 model
## rpl36: 1 model
######################
rv.lrt %>% 
    filter(p.correct >= 0.01) %>% 
    select(dataset, model, type, p.correct)
# A tibble: 5 x 4
# Groups:   dataset, type, model [5]
#  dataset  model   type  p.correct
#   <fctr>  <chr> <fctr>      <dbl>
#  1    petN gcpREV chloro 0.19425957
#  2    petN    JTT chloro 0.07801858
#  3    petN     LG chloro 0.06825117
#  4    petN    WAG chloro 0.13412816
#  5    psbF gcpREV chloro 0.02978061
#  6    psbF    JTT chloro 0.02781492
#  7    psbF  mtVer chloro 0.01987518
#  8    psbF    WAG chloro 0.02432178
#  9    psbL     LG chloro 0.01401419
# 10   rpl36    WAG chloro 1.00000000















