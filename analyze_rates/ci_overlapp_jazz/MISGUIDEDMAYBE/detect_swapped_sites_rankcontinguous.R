##################################################################
## Kindly note that this script takes absolutely forever to run ##
##################################################################
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
#### We rank MLEs by ensuring that they are continguous; ties receive SAME RANK and as values increase there are no gaps. Makes standard deviation of ranks MUCH more reliable by removing artifacts of ties.
rates %>% 
    ungroup() %>%
    group_by(dataset, model) %>%
    select(-site, -Lower, -Upper) %>%
    unique() %>%
    arrange(MLE) %>%
    mutate(MLErank = 1:n(), medianrank = median(MLErank)) %>% 
    left_join(rates) %>%
    ungroup() -> rates.rank
    
cirank <- 


cirank <- read_csv("CIrank_FINAL.csv", col_types = list(col_character(), col_character(), col_integer(), col_double(), col_double())) 

rates.ranks2 <- left_join(rates.rank2, cirank)
rates.ranks2 %>% filter(dataset == "13PK_A", model =="WAG") %>% arrange(MLE) %>% head(120) %>% print.data.frame()
#### Detect all sites where all CIs do not fully contain all MLEs ####
sites.outside.ci <- data.frame(dataset = character(), site = numeric())
for (ds in unique(rates$dataset)){
    print(ds)
    rates.ranks %>% filter(dataset == ds) -> thisdata
    for (x in unique(thisdata$site)) {
        thisdata %>% filter(site == x) -> thissite
        total <- 0
        for (m in thissite$model){
            raterank <- thissite$MLErank[thissite$model == m]
            total <- total + sum(thissite$LowerRank > raterank) + sum(thissite$LowerRank < raterank)
        }
        if (total > 0 ){
            temp <- data.frame(dataset = ds, site = x)
            sites.outside.ci  <- rbind(sites.outside.ci, temp)
        }
    }
}




##### Detect sites for which ALL conditions stand:
######### 1) All rates are **upper bounded**. This excludes runaway sites with unreliable MLEs.
######### 2) At least one rate is swapped above/below the median rank. 
######### 3) The swapped rate is meaningful: It is not unbounded and it is more than 1 standard deviation away from the gene's median rank
### We inner_join at the end because ranks, median, sd will be incorrect otherwise. Must consider all sites for these measures.
rates %>% 
    ungroup() %>%
    group_by(dataset, model) %>% 
    mutate(MLErank = rank(MLE),
           medianrank = median(MLErank),
           rankside = sign(MLErank- medianrank),
           sdrank = sd(MLErank),
           withinsd = ifelse( abs(medianrank - MLErank) <= sdrank, TRUE, FALSE),
           unbounded = Upper >= upper.unbounded) -> rates.sd.rank

rates.sd.rank %>%
    inner_join(sites.outside.ci) %>% 
    ungroup() %>%
    filter(unbounded==FALSE) %>% 
    group_by(dataset, site) %>%
    mutate(n = n()) %>% 
    mutate( moderankside = ifelse(sum(rankside == 1) > sum(rankside == -1), 1, -1 ),
            all.unbounded = Upper >= upper.unbounded,
            goodsite = ifelse(rankside == moderankside, TRUE, FALSE), ## true when rank matches mode rank
            goodsite = ifelse(goodsite | withinsd, TRUE, FALSE)) %>% ## true when rank matches mode rank OR we are within 1 sd (so rank switch within 1 sd is ok)
    summarize(all.unbounded = sum(all.unbounded),
              goodsite = sum(goodsite),
              goodsite2 = abs(sum(rankside)),
              n =mean(n)) %>%
    filter(all.unbounded != n) %>%
    filter(goodsite != n) %>%
    filter(goodsite2 != n) %>%
    select(dataset, site) -> susceptible.sites

susceptible.sites %>% 
    left_join(rates.sd.rank) %>% 
    write_csv("swapped_sites.csv")






















########################################################################################################################################################




#### Detect all sites where all CIs do not fully contain all MLEs ####
sites.outside.ci <- data.frame(dataset = character(), site = numeric())
for (ds in unique(rates$dataset)){
    print(ds)
    rates %>% filter(dataset == ds) -> thisdata
    for (x in unique(thisdata$site)) {
        thisdata %>% filter(site == x) -> thissite
        total <- 0
        for (m in thissite$model){
            rate <- thissite$MLE[thissite$model == m]
            total <- total + sum(thissite$Lower > rate) + sum(thissite$Upper < rate)
        }
        if (total > 0 ){
            temp <- data.frame(dataset = ds, site = x)
            sites.outside.ci  <- rbind(sites.outside.ci, temp)
        }
    }
}

##### Detect sites for which ALL conditions stand:
######### 1) All rates are **upper bounded**. This excludes runaway sites with unreliable MLEs.
######### 2) At least one rate is swapped above/below the median rank. 
######### 3) The swapped rate is meaningful: It is not unbounded and it is more than 1 standard deviation away from the gene's median rank
### We inner_join at the end because ranks, median, sd will be incorrect otherwise. Must consider all sites for these measures.
rates %>% 
    ungroup() %>%
    group_by(dataset, model) %>% 
    mutate(MLErank = rank(MLE),
           medianrank = median(MLErank),
           rankside = sign(MLErank- medianrank),
           sdrank = sd(MLErank),
           withinsd = ifelse( abs(medianrank - MLErank) <= sdrank, TRUE, FALSE),
           unbounded = Upper >= upper.unbounded) -> rates.sd.rank

rates.sd.rank %>%
    inner_join(sites.outside.ci) %>% 
    ungroup() %>%
    filter(unbounded==FALSE) %>% 
    group_by(dataset, site) %>%
    mutate(n = n()) %>% 
    mutate( moderankside = ifelse(sum(rankside == 1) > sum(rankside == -1), 1, -1 ),
            all.unbounded = Upper >= upper.unbounded,
            goodsite = ifelse(rankside == moderankside, TRUE, FALSE), ## true when rank matches mode rank
            goodsite = ifelse(goodsite | withinsd, TRUE, FALSE)) %>% ## true when rank matches mode rank OR we are within 1 sd (so rank switch within 1 sd is ok)
    summarize(all.unbounded = sum(all.unbounded),
              goodsite = sum(goodsite),
              goodsite2 = abs(sum(rankside)),
              n =mean(n)) %>%
    filter(all.unbounded != n) %>%
    filter(goodsite != n) %>%
    filter(goodsite2 != n) %>%
    select(dataset, site) -> susceptible.sites

susceptible.sites %>% 
    left_join(rates.sd.rank) %>% 
    write_csv("swapped_sites.csv")






