###
### SJS
### This script processes the rate inferences to obtain correlations
### Generates the following files:
###   1) rate_inferences_all.csv, all rate inferences. **No normalization or filtering applied**
###   2) correlations_gamma_noRV.csv, per-gene correlations (both Pearson on log-transformed data and Spearman) for each model with and without gamma, i.e. WAG norv vs WAG gamma, etc.
###   3) correlations_between_models.csv, per-gene correlations (both Pearson on log-transformed data and Spearman) between models, specifically norv inferences
###

require(tidyverse)
require(purrr)
require(broom)
mito <- c("ATP6", "ATP8", "COX1", "COX2", "COX3", "CYTB", "ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6") 


###### Read in organelle rates #######
fpath = "../rate_inference/organelle_data-inference/"
files <- dir(path = fpath, pattern = "*csv")
data_frame(filename = files) %>%
  mutate(file_contents = map(filename,
           ~ read_csv(file.path(fpath, .)))) %>%
  unnest() %>%
  separate(filename, c("dataset", "byebye1", "model", "byebye2", "byebye3"), sep="\\.") %>%
  dplyr::select(dataset, model, site, MLE, Lower, Upper) %>%
  mutate(type = ifelse(dataset %in% mito, "mito", "chloro")) -> raw.organelle.rates

###### Read in enzyme rates and merge with organelle to create single tibble of rates #######
fpath = "../rate_inference/enzyme_data-inference/"
files <- dir(path = fpath, pattern = "*csv")
data_frame(filename = files) %>%
  mutate(file_contents = map(filename,
           ~ read_csv(file.path(fpath, .)))) %>%
  unnest() %>%
  separate(filename, c("dataset", "byebye1", "model", "byebye2", "byebye3"), sep="\\.") %>%
  dplyr::select(dataset, model, site, MLE, Lower, Upper) %>% 
  mutate(type = "enzyme") %>%
  rbind(raw.organelle.rates) -> raw.rates

##### SAVE RATES TO rate_inferences_all.csv
write_csv(raw.rates, "rate_inferences_all.csv")


# 
# ##### For fun, confirm that sites with rate = 0 are consistently 0 for all methods
# > raw.rates %>% 
# +     group_by(dataset, site) %>% 
# +     tally(MLE == 0 & Lower == 0) %>% 
# +     filter(n!=12, n!=0)
# A tibble: 0 x 3
# Groups:   dataset [0]
# ... with 3 variables: dataset <chr>, site <int>, n <int>


### Spread the rates, adding 1e-8 all around to allow median normalization sans NA's and for downstream logging
raw.rates %>% 
    group_by(dataset, model) %>%
    select(-Lower,-Upper) %>%
    mutate(MLE = MLE + 1e-8) %>%
    mutate(MLE = MLE/median(MLE)) %>%
    spread(model, MLE) -> rates.mednorm

    
##### Rank correlations for noRV - Gamma
rates.mednorm %>%
    group_by(dataset, type) %>%
    do(LG = cor(.$`LG-G`, .$`LG-No`, method = "spearman"),
       WAG = cor(.$`WAG-G`, .$`WAG-No`, method = "spearman"),
       JTT = cor(.$`JTT-G`, .$`JTT-No`, method = "spearman"),
       mtMet = cor(.$`mtMet-G`, .$`mtMet-No`, method = "spearman"),
       cpREV = cor(.$`cpREV-G`, .$`cpREV-No`, method = "spearman"),
       JC69 = cor(.$`JC69-G`, .$`JC69-No`, method = "spearman")) %>%
    unnest(LG, WAG, JTT, mtMet, cpREV, JC69) %>%
    gather(model, r, LG:JC69) %>%
    mutate(correlation.type = "Spearman") -> rank.corr.g.norv


##### Pearson on transformed data correlations for noRV - Gamma, and merge with spearman correlations to create a single data frame
rates.mednorm %>%
  group_by(dataset, type) %>%
  do(LG = cor(log(.$`LG-G`), log(.$`LG-No`)),
     WAG = cor(log(.$`WAG-G`), log(.$`WAG-No`)),
     JTT = cor(log(.$`JTT-G`), log(.$`JTT-No`)),
     mtMet = cor(log(.$`mtMet-G`), log(.$`mtMet-No`)),
     cpREV = cor(log(.$`cpREV-G`), log(.$`cpREV-No`)),
     JC69 = cor(log(.$`JC69-G`), log(.$`JC69-No`))) %>%
    unnest(LG, WAG, JTT, mtMet, cpREV, JC69) %>%
    gather(model, r, LG:JC69) %>%
    mutate(correlation.type = "Pearson") %>%
  rbind(rank.corr.g.norv) -> corr.gamma.norv


##### SAVE NoRV-Gamma correlations to correlations_gamma_noRV.csv
write_csv(corr.gamma.norv, "correlations_gamma_noRV.csv")




##### Rank correlations between each pair of models. Note that I'm sure there's a better way to do this, but I don't know what that way is.
rates.mednorm %>%
    group_by(dataset, type) %>%
    do(`LG-WAG`      = cor(.$`LG-No`, .$`WAG-No`, method = "spearman"),
       `LG-JTT`      = cor(.$`LG-No`, .$`JTT-No`, method = "spearman"),
       `LG-mtMet`    = cor(.$`LG-No`, .$`mtMet-No`, method = "spearman"),
       `LG-cpREV`    = cor(.$`LG-No`, .$`cpREV-No`, method = "spearman"),
       `LG-JC69`     = cor(.$`LG-No`, .$`JC69-No`, method = "spearman"),
       `WAG-JTT`     = cor(.$`WAG-No`, .$`JTT-No`, method = "spearman"),
       `WAG-mtMet`   = cor(.$`WAG-No`, .$`mtMet-No`, method = "spearman"),
       `WAG-cpREV`   = cor(.$`WAG-No`, .$`cpREV-No`, method = "spearman"),
       `WAG-JC69`    = cor(.$`WAG-No`, .$`JC69-No`, method = "spearman"),
       `JTT-mtMet`   = cor(.$`JTT-No`, .$`mtMet-No`, method = "spearman"),
       `JTT-cpREV`   = cor(.$`JTT-No`, .$`cpREV-No`, method = "spearman"),
       `JTT-JC69`    = cor(.$`JTT-No`, .$`JC69-No`, method = "spearman"),
       `mtMet-cpREV` = cor(.$`mtMet-No`, .$`cpREV-No`, method = "spearman"),
       `mtMet-JC69`  = cor(.$`mtMet-No`, .$`JC69-No`, method = "spearman"),
       `cpREV-JC69`  = cor(.$`cpREV-No`, .$`JC69-No`, method = "spearman")) %>%
    unnest(`LG-WAG`, `LG-JTT`, `LG-mtMet`, `LG-cpREV`, `LG-JC69`, `WAG-JTT`, `WAG-mtMet`, `WAG-cpREV`, `WAG-JC69`, `JTT-mtMet`, `JTT-cpREV`, `JTT-JC69`, `mtMet-cpREV`, `mtMet-JC69`, `cpREV-JC69`) %>% 
    gather(comparison, r, `LG-WAG`:`cpREV-JC69`) %>%
    mutate(correlation.type = "Spearman") -> corr.models.rank



##### Pearson on log-tranformed data correlations between each pair of models, and merge with spearman to create single data frame
rates.mednorm %>%
    group_by(dataset, type) %>%
    do(`LG-WAG`      = cor(log(.$`LG-No`), log(.$`WAG-No`)),
       `LG-JTT`      = cor(log(.$`LG-No`), log(.$`JTT-No`)),
       `LG-mtMet`    = cor(log(.$`LG-No`), log(.$`mtMet-No`)),
       `LG-cpREV`    = cor(log(.$`LG-No`), log(.$`cpREV-No`)),
       `LG-JC69`     = cor(log(.$`LG-No`), log(.$`JC69-No`)),
       `WAG-JTT`     = cor(log(.$`WAG-No`), log(.$`JTT-No`)),
       `WAG-mtMet`   = cor(log(.$`WAG-No`), log(.$`mtMet-No`)),
       `WAG-cpREV`   = cor(log(.$`WAG-No`), log(.$`cpREV-No`)),
       `WAG-JC69`    = cor(log(.$`WAG-No`), log(.$`JC69-No`)),
       `JTT-mtMet`   = cor(log(.$`JTT-No`), log(.$`mtMet-No`)),
       `JTT-cpREV`   = cor(log(.$`JTT-No`), log(.$`cpREV-No`)),
       `JTT-JC69`    = cor(log(.$`JTT-No`), log(.$`JC69-No`)),
       `mtMet-cpREV` = cor(log(.$`mtMet-No`), log(.$`cpREV-No`)),
       `mtMet-JC69`  = cor(log(.$`mtMet-No`), log(.$`JC69-No`)),
       `cpREV-JC69`  = cor(log(.$`cpREV-No`), log(.$`JC69-No`))) %>%
    unnest(`LG-WAG`, `LG-JTT`, `LG-mtMet`, `LG-cpREV`, `LG-JC69`, `WAG-JTT`, `WAG-mtMet`, `WAG-cpREV`, `WAG-JC69`, `JTT-mtMet`, `JTT-cpREV`, `JTT-JC69`, `mtMet-cpREV`, `mtMet-JC69`, `cpREV-JC69`) %>% 
    gather(comparison, r, `LG-WAG`:`cpREV-JC69`) %>%
    mutate(correlation.type = "Pearson") %>%
  rbind( corr.models.rank ) -> corr.models


##### SAVE between model correlations to correlations_between_models.csv
write_csv(corr.models, "correlations_between_models.csv")




#####################################################################
########################## Plotting #################################

theme_set(theme_classic())
corrtype.levels <- c("Pearson", "Spearman")
comparison.levels <- c("LG-WAG", "LG-JTT", "LG-mtMet", "LG-cpREV", "WAG-JTT", "WAG-mtMet", "WAG-cpREV", "JTT-mtMet", "JTT-cpREV", "mtMet-cpREV", "LG-JC69", "WAG-JC69", "JTT-JC69", "mtMet-JC69", "cpREV-JC69")
model.levels <- c("LG", "WAG", "JTT", "mtMet", "cpREV", "JC69")
type.levels <- c("enzyme", "chloro", "mito")
type.labels <- c("Enzyme", "Chloroplast", "Mitochondria")



########## Correlations between Gamma, No RV ##########
corr.gamma.norv$model <- factor(corr.gamma.norv$model, levels=model.levels)
corr.gamma.norv$type <- factor(corr.gamma.norv$type, levels=type.levels, labels=type.labels)
corr.gamma.norv$correlation.type <- factor(corr.gamma.norv$correlation.type, levels=corrtype.levels)
corr.gamma.norv %>% 
    ggplot(aes(x = model, y = r, fill = type))+
    geom_point(pch=21, position = position_jitterdodge()) + 
    facet_wrap(~correlation.type, nrow=2) + theme(strip.background = element_rect(fill = "grey90")) + 
    ylab("Correlation") + xlab("Model") + scale_fill_discrete(name = "Dataset") -> rv.jitter.plot
### Outliers are rpl36 (WAG, JTT, JC), and psbL (JC)
ggsave("gamma_norv_jitter.pdf", rv.jitter.plot)




########## Jitter plot of correlations between models ##########
corr.models$comparison <- factor(corr.models$comparison, levels=comparison.levels)
corr.models$type <- factor(corr.models$type, levels=type.levels, labels=type.labels)
corr.models$correlation.type <- factor(corr.models$correlation.type, levels=corrtype.levels, labels=corrtype.labels)
corr.models %>% 
    ggplot(aes(x = comparison, y = r, fill = type))+
    geom_point(pch=21, position = position_jitterdodge()) + 
    facet_wrap(~correlation.type, nrow=2) + theme(axis.text.x = element_text(angle=10, face = "bold", size=7), strip.background = element_rect(fill = "grey90")) + 
    ylab("Correlation") + xlab("Compared models") + scale_fill_discrete(name = "Dataset") -> between.jitter.plot
#    guides(fill = guide_legend(override.aes = list(size = 2))) + 
ggsave("between_models_jitter.pdf", between.jitter.plot)
