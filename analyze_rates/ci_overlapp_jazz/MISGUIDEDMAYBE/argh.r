library(tidyverse)
options(tibble.width = Inf)

enzyme <- read_csv("summarized_data/rate_inferences_enzymes.csv")
virus <- read_csv("summarized_data/rate_inferences_virus.csv")
org <- read_csv("summarized_data/rate_inferences_organelles.csv")

models <- c("WAG", "LG", "JTT", "JC69", "gcpREV", "mtVer")
upper.unbounded <- 1000

enzyme %>%
    rbind(virus) %>%
    rbind(org) %>% 
    filter(rv == "No") %>%
    filter(model %in% models) %>%
    select(dataset, model, site, MLE, Lower, Upper) -> rates

#### Convert MLE CI to *rank* CI by finding, for each bound, the rank of the closest MLE (per dataset per model)
#### If there are error ties, select min rank for Lower, max rank for Upper
rates %>% 
    ungroup() %>%
    group_by(dataset, model) %>% 
    mutate(MLErank = rank(MLE)) -> rates.rank

cirank <- tibble(dataset = character(), model = character(), site = integer(), LowerRank = integer(), UpperRank = integer())
for (d in all.datasets){
    print(d)
    rates.rank %>% filter(dataset == d) -> rates.rank.dataset
    rates.rank.dataset$dataset <- as.character(rates.rank.dataset$dataset)
        for (i in unique(rates.rank.dataset$site)){
            rates.rank.dataset %>% 
                filter(site == i) %>%
                ungroup() %>% 
                select(dataset, model, Lower, Upper) -> x

            rates.rank.dataset %>% 
                ungroup() %>% 
                select(-Lower, -Upper) %>%
                group_by(model) %>% 
                left_join(x) -> hinei
            
            hinei %>% 
                mutate(errlow = abs(MLE -Lower)/MLE) %>% 
                filter(errlow==min(errlow)) %>%
                filter(MLErank == min(MLErank)) %>%
                mutate(LowerRank = MLErank) %>%
                select(dataset, model, LowerRank)-> x2

             hinei %>%         
                mutate(errhigh = abs(MLE - Upper)/MLE) %>%
                filter(MLErank == max(MLErank)) %>%
                mutate(UpperRank = MLErank) %>%
                select(dataset, model, UpperRank) %>% 
                left_join(x2, by = c("dataset", "model")) %>%
                mutate(site = i) %>% 
                unique() -> temp
            cirank <- bind_rows(cirank, temp)
        }
}
write_csv(cirank, "CIrank_odpaam.csv")
