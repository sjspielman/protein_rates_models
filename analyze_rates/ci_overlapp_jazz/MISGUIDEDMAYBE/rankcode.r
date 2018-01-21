
cirank <- tibble(dataset = character(), model = character(), site = integer(), LowerRank = integer(), UpperRank = integer())
all.datasets <- unique(rates.rank$dataset)
for (d in unique(dataset)){
    print(d)
    rates.rank %>% filter(dataset == d) -> rates.rank.dataset
    rates.rank.dataset$dataset <- as.character(rates.rank.dataset$dataset)
    minrank <- min(rates.rank.dataset$MLErank)
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
                
            if (sum(x$Lower) == 0){
                x %>% select(dataset, model) %>% mutate(LowerRank = minrank) -> x2
            }
            else {                   
                hinei %>% 
                    mutate(errlow = abs(MLE -Lower)/MLE) %>% 
                    filter(errlow==min(errlow)) %>%
                    filter(MLErank == min(MLErank)) %>%
                    mutate(LowerRank = MLErank) %>%
                    select(dataset, model, LowerRank)-> x2
            }
            
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


lowerrank <- tibble(dataset = character(), model = character(), site = integer(), LowerRank = integer())

for (d in all.datasets)
{
    print(d)
    rates.rank %>% filter(dataset == d) -> rates.rank.dataset
    rates.rank.dataset$dataset <- as.character(rates.rank.dataset$dataset)
    minrank <- min(rates.rank.dataset$MLErank)
    for (i in unique(rates.rank.dataset$site)){
        rates.rank.dataset %>% 
            filter(site == i) %>%
            select(dataset, model, Lower)-> x
        if (sum(x$Lower) == 0){
            print(paste0("   ",i))
            x %>% mutate(LowerRank = minrank, site=i) %>%select(dataset, model, site, LowerRank) -> temp
            lowerrank <- bind_rows(lowerrank, temp)
        }
    }
}



upperrank <- tibble(dataset = character(), model = character(), site = integer(), UpperRank = integer())

for (d in all.datasets)
{
    print(d)
    rates.rank %>% filter(dataset == d) -> rates.rank.dataset
    rates.rank.dataset$dataset <- as.character(rates.rank.dataset$dataset)
    maxrank <- max(rates.rank.dataset$MLErank)
    for (i in unique(rates.rank.dataset$site)){
        rates.rank.dataset %>% 
            filter(site == i) %>%
            select(dataset, model, Lower)-> x
        if (sum(x$Lower) == 0){
            print(paste0("   ",i))
            x %>% mutate(LowerRank = minrank, site=i) %>%select(dataset, model, site, LowerRank) -> temp
            lowerrank <- bind_rows(lowerrank, temp)
        }
    }
}




write_csv(cirank, "CIrank_odpaam.csv")