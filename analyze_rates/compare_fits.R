library(tidyverse)
datadir <- "summarized_data/"
model.fits <- read.csv(paste0(datadir, "model_fits.csv"))

model.fits %>% 
    group_by(dataset, type) %>%
    mutate(AICc_rank = rank(AICc)) %>% 
    mutate(deltaAICc = min(AICc) - AICc) -> model.aicc


model.aicc %>% 
    filter(AICc_rank == 1) %>% 
    group_by(model, type) %>%
    tally() %>%
    arrange(type, desc(n))
#  1  gcpREV-G chloro    40
#  2     JTT-G chloro    14
#  3   mtVer-G chloro     3
#  4 gcpREV-No chloro     1
#  5      LG-G chloro     1
#  6      LG-G enzyme    62
#  7     WAG-G enzyme    36
#  8     JTT-G enzyme     2
#  9   mtVer-G   mito    11
# 10  gcpREV-G   mito     1
# 11     JTT-G   mito     1
# 12     JTT-G  virus    53
# 13   mtVer-G  virus    19
# 14  gcpREV-G  virus    15
# 15    JTT-No  virus    13
# 16  mtVer-No  virus     9
# 17 gcpREV-No  virus     6
# 18      LG-G  virus     1
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

rv.lrt %>% filter(p.correct >= 0.01) %>% group_by(type) %>% tally()
#1 chloro     9
#2  virus   271

# 1 chloro gcpREV    petN    37
# 2 chloro    JTT    petN    37
# 3 chloro     LG    petN    37
# 4 chloro  mtVer    petN    37
# 5 chloro    WAG    petN    37
# 6 chloro gcpREV    psbF    50
# 7 chloro   JC69    psbL    68
# 8 chloro    JTT    psbL    68
# 9 chloro     LG    psbL    68
# 
rv.lrt %>% filter(p.correct >= 0.01)  %>% filter(type == "chloro") %>% ungroup() %>% select(dataset) %>% unique()
# 3 chloroplast datasets
rv.lrt %>% filter(p.correct >= 0.01)  %>% filter(type == "virus") %>% ungroup() %>% select(dataset) %>% unique()
# 55 virus datasets

###############################################################################################












