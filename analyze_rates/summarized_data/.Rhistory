library(tidyverse)
info <- read_csv("alignment_information.csv")
info$dataset
rates <- read_csv("rate_inferences_enzymes.csv")

rates %>% 
    filter(rv == "No", !(model %in% c("permWAG", "RAND")) ) %>% 
    select(dataset, site, MLE, model) %>% 
    left_join(info) -> here 
    
here %>% 
    group_by(dataset, model) %>%
    summarize(rho = cor(MLE, true_entropy, method = "spearman")) %>%
    ungroup() %>%
    group_by(model) %>%
    summarize(minrho = min(rho), maxrho = max(rho), meanrho = mean(rho), sdrho = sd(rho))
    
    ggplot(aes(x = rho)) + geom_histogram() + facet_grid(~model)
    
here %>%
    filter(dataset %in% c("13PK_A", "1A16_A", "1A4L_A", "1A50_B", "1A8H_A"), MLE>1e-5) %>%
    mutate(MLE = ifelse(MLE >= 100, 100, MLE)) %>%
    group_by(dataset, model) %>%
    ggplot(aes(x = MLE, y = true_entropy)) + geom_point() + 
    facet_grid(model~dataset) + scale_x_log10() 
    

   
    
    summarize(rho = cor(MLE, true_entropy, method = "spearman"))
