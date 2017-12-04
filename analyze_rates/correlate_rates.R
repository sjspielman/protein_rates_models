###
### SJS
### This script processes the rate inferences to obtain correlations
### Generates the following files:
###   1) rate_inferences_all.csv, all rate inferences. **No normalization or filtering applied**
###   2) cleaned_rate_inferences_all.csv, rate inferences where unreliable rates are removed (for all methods), specifically <Lower == 0, Upper == 10000, MLE==0)> or any Upper==10000. Generally unbounded rates are not super.
###   3) correlations_gamma_noRV.csv, per-gene correlations (both Pearson on log-transformed data and Spearman) for each model with and without gamma, i.e. WAG norv vs WAG gamma, etc.
###   4) correlations_between_models.csv, per-gene correlations (both Pearson on log-transformed data and Spearman) between models, specifically norv inferences
###

require(tidyverse)
require(purrr)
require(broom)
mito <- c("ATP6", "ATP8", "COX1", "COX2", "COX3", "CYTB", "ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6") 


###### Read in organelle rates #######
fpath = "organelle_data-inference/"
files <- dir(path = fpath, pattern = "*csv")
data_frame(filename = files) %>%
  mutate(file_contents = map(filename,
           ~ read_csv(file.path(fpath, .)))) %>%
  unnest() %>%
  separate(filename, c("dataset", "byebye1", "model", "byebye2", "byebye3"), sep="\\.") %>%
  dplyr::select(dataset, model, site, MLE, Lower, Upper) %>%
  mutate(type = ifelse(dataset %in% mito, "mito", "chloro")) -> raw.organelle.rates

###### Read in enzyme rates and merge with organelle to create single tibble of rates #######
fpath = "enzyme_data-inference/"
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

##### Confirm that sites with rate = 0 are consistently 0 for all methods
raw.rates %>% 
    group_by(dataset, site) %>% 
    tally((MLE == 0 & Lower == 0) | Upper == 10000) %>% 
    filter(n!=12, n!=0)
# A tibble: 0 x 3
# Groups:   dataset [0]
# ... with 3 variables: dataset <chr>, site <int>, n <int>


#### Remove all MLE = 0 and change all MLE>=1e2 to 1e2, as these rates are necessarily unstable. Normalize raw rates by median (per gene), and spread ####
raw.rates %>%
  mutate(keep = ifelse(, TRUE,
  
  filter(MLE > 0 & Lower > 0) %>%
  filter(Upper != 10000 & Lower 
  mutate(MLE = ifelse(MLE >= 1e2, 1e2, MLE)) %>%
  ungroup() %>%
  group_by(dataset, model) %>%
  mutate(MLE = MLE/median(MLE)) %>%
  spread(model, MLE) -> rates.nonzero.mednorm
  

##### Rank correlations for noRV - Gamma
rates.nonzero.mednorm %>%
    group_by(dataset) %>%
    do(corLG = cor(.$`LG-G`, .$`LG-No`, method = "spearman"),
       corWAG = cor(.$`WAG-G`, .$`WAG-No`, method = "spearman"),
       corJTT = cor(.$`JTT-G`, .$`JTT-No`, method = "spearman"),
       cormtMAM = cor(.$`mtMAM-G`, .$`mtMAM-No`, method = "spearman"),
       corcpREV = cor(.$`cpREV-G`, .$`cpREV-No`, method = "spearman"),
       corJC69 = cor(.$`JC69-G`, .$`JC69-No`, method = "spearman")) %>%
    unnest(corLG, corWAG, corJTT, cormtMAM, corcpREV, corJC69) %>%
    mutate(correlation = "spearman") -> rank.corr.g.norv


##### Pearson on transformed data correlations for noRV - Gamma, and merge with spearman correlations to create a single data frame
rates.nonzero.mednorm %>%
  group_by(dataset) %>%
  do(corLG = cor(log(.$`LG-G`), log(.$`LG-No`)),
    corWAG = cor(log(.$`WAG-G`), log(.$`WAG-No`)),
    corJTT = cor(log(.$`JTT-G`), log(.$`JTT-No`)),
    cormtMAM = cor(log(.$`mtMAM-G`), log(.$`mtMAM-No`)),
    corcpREV = cor(log(.$`cpREV-G`), log(.$`cpREV-No`)),
    corJC69 = cor(log(.$`JC69-G`), log(.$`JC69-No`))) %>%
  unnest(corLG, corWAG, corJTT, cormtMAM, corcpREV, corJC69) %>%
  mutate(correlation = "pearson.on.log") %>%
  rbind(rank.corr.g.norv) -> corr.gamma.norv


##### SAVE NoRV-Gamma correlations to correlations_gamma_noRV.csv
write_csv(corr.gamma.norv, "correlations_gamma_noRV.csv")




##### Rank correlations between each pair of models. Note that I'm sure there's a better way to do this, but I don't know what that way is.
rates.nonzero.mednorm %>%
    group_by(dataset, type) %>%
    do(LGWAG = cor(.$`LG-No`, .$`WAG-No`, method = "spearman"),
       LGJTT = cor(.$`LG-No`, .$`JTT-No`, method = "spearman"),
       LGmtMAM = cor(.$`LG-No`, .$`mtMAM-No`, method = "spearman"),
       LGcpREV = cor(.$`LG-No`, .$`cpREV-No`, method = "spearman"),
       LGJC = cor(.$`LG-No`, .$`JC69-No`, method = "spearman"),
       WAGJTT = cor(.$`WAG-No`, .$`JTT-No`, method = "spearman"),
       WAGmtMAM = cor(.$`WAG-No`, .$`mtMAM-No`, method = "spearman"),
       WAGcpREV = cor(.$`WAG-No`, .$`cpREV-No`, method = "spearman"),
       WAGJC = cor(.$`WAG-No`, .$`JC69-No`, method = "spearman"),
       JTTmtMAM = cor(.$`JTT-No`, .$`mtMAM-No`, method = "spearman"),
       JTTcpREV = cor(.$`JTT-No`, .$`cpREV-No`, method = "spearman"),
       JTTJC = cor(.$`JTT-No`, .$`JC69-No`, method = "spearman"),
       mtMAMcpREV = cor(.$`mtMAM-No`, .$`cpREV-No`, method = "spearman"),
       mtMAMJC = cor(.$`mtMAM-No`, .$`JC69-No`, method = "spearman"),
       cpREVJC = cor(.$`cpREV-No`, .$`JC69-No`, method = "spearman")) %>%
    unnest(LGWAG, LGJTT, LGmtMAM, LGcpREV, LGJC, WAGJTT, WAGmtMAM, WAGcpREV, WAGJC, JTTmtMAM, JTTcpREV, JTTJC, mtMAMcpREV, mtMAMJC, cpREVJC) %>% 
    gather(comparison, r, LGWAG:cpREVJC) %>%
    mutate(correlation_type = "spearman") -> corr.models.rank



##### Pearson on log-tranformed data correlations between each pair of models, and merge with spearman to create single data frame
rates.nonzero.mednorm %>%
    group_by(dataset, type) %>%
    do(LGWAG = cor(log(.$`LG-No`), log(.$`WAG-No`)),
       LGJTT = cor(log(.$`LG-No`), log(.$`JTT-No`)),
       LGmtMAM = cor(log(.$`LG-No`), log(.$`mtMAM-No`)),
       LGcpREV = cor(log(.$`LG-No`), log(.$`cpREV-No`)),
       LGJC = cor(log(.$`LG-No`), log(.$`JC69-No`)),
       WAGJTT = cor(log(.$`WAG-No`), log(.$`JTT-No`)),
       WAGmtMAM = cor(log(.$`WAG-No`), log(.$`mtMAM-No`)),
       WAGcpREV = cor(log(.$`WAG-No`), log(.$`cpREV-No`)),
       WAGJC = cor(log(.$`WAG-No`), log(.$`JC69-No`)),
       JTTmtMAM = cor(log(.$`JTT-No`), log(.$`mtMAM-No`)),
       JTTcpREV = cor(log(.$`JTT-No`), log(.$`cpREV-No`)),
       JTTJC = cor(log(.$`JTT-No`), log(.$`JC69-No`)),
       mtMAMcpREV = cor(log(.$`mtMAM-No`), log(.$`cpREV-No`)),
       mtMAMJC = cor(log(.$`mtMAM-No`), log(.$`JC69-No`)),
       cpREVJC = cor(log(.$`cpREV-No`), log(.$`JC69-No`))) %>%
  unnest(LGWAG, LGJTT, LGmtMAM, LGcpREV, LGJC, WAGJTT, WAGmtMAM, WAGcpREV, WAGJC, JTTmtMAM, JTTcpREV, JTTJC, mtMAMcpREV, mtMAMJC, cpREVJC) %>% 
    gather(comparison, r, LGWAG:cpREVJC) %>%
    mutate(correlation_type = "pearson.on.log") %>%
  rbind( corr.models.rank ) -> corr.models


##### SAVE between model correlations to correlations_between_models.csv
write_csv(corr.gamma.norv, "correlations_between_models.csv")


