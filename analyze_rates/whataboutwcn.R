wcn   <- read_csv("enzyme_wcnSC.csv")
rates <- read_csv("rate_inferences_all.csv")

rates %>% 
    filter(type == "enzyme") %>% 
    right_join(wcn) %>% 
    filter(MLE > 0) %>%
    select(-Lower, -Upper, -type) %>% 
    ungroup() %>% 
    group_by(dataset, model) %>%
    summarize(rho = cor(MLE, wcnSC, method = "spearman"))  -> rate.wcn.correlations


rate.wcn.correlations %>% 
    group_by(dataset) %>% 
    mutate(rhorank = rank(rho)) %>%
    filter(rhorank == 1) %>%
    select(dataset, model) %>%
    ungroup() -> rhoranks
# A tibble: 6 x 2
   model     n
   <chr> <int>
1 gcpREV    20
2   JC69    25
3    JTT     5
4     LG     9
5  mtVer    33
6    WAG     6

model.aicc %>% 
    filter(AICc_rank == 1, type == "enzyme") %>% 
    select(dataset, model) %>%
    mutate(aiccmodel = str_split(model, "-")) %>% 
    unnest() %>%
    filter(aiccmodel != "G") %>%
    select(-model) %>%
    left_join(rhoranks) %>% filter(aiccmodel == model)
# A tibble: 6 x 3
# Groups:   dataset [6]
#   dataset aiccmodel model
#     <chr>     <chr> <chr>
# 1  1A50_B       WAG   WAG
# 2  1CTN_A       WAG   WAG
# 3  1QI9_A       WAG   WAG
# 4  1UQT_A       WAG   WAG
# 5  1YCF_A       WAG   WAG
# 6  206L_A       WAG   WAG
########### 6/100 datasets show highest correlation for the best fitting model.

rate.wcn.correlations %>% 
    group_by(dataset) %>%
    summarize(sd.rho = sd(rho), range.rho = max(rho) - min(rho)) %>% filter(range.rho >= 0.05)
    summarize(max(sd.rho), max(range.rho))
# # A tibble: 1 x 2
#   `max(sd.rho)` `max(range.rho)`
#           <dbl>            <dbl>
# 1     0.04119974        0.1106538
