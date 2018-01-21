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

#### Detect all sites where all CIs do not fully contain all MLEs ####
sites.outside.ci <- data.frame(dataset = character(), site = numeric(), model = character(), count.outside = numeric())
for (ds in unique(rates$dataset)){
    print(ds)
    rates %>% filter(dataset == ds) -> thisdata
    for (x in unique(thisdata$site)) 
    {
        thisdata %>% filter(site == x) -> thissite
        for (m in thissite$model)
        {
            total <- 0
            thissite %>% filter(model == m) -> mod
            anti_join(thissite, mod, by = c("dataset", "model", "site", "MLE", "Lower", "Upper")) -> notmod
            modmle <- mod$MLE
            total <- sum(modmle < notmod$Lower | modmle > notmod$Upper)
            sites.outside.ci <- bind_rows(sites.outside.ci, tibble(dataset = ds, site = x, model = m, count.outside=total))
        }
    }
}

rates %>% 
    filter(dataset == "1EG7_A", site == 203) %>% 
    arrange(MLE) %>% 
    mutate(y = 1:n()) %>% 
    ggplot() + 
    geom_segment(aes(x = Lower, xend = Upper, y = y, yend = y, color=model), size=1.25) + 
    geom_point(aes(x = MLE, y = y))+
    xlab("95% CI") + ylab("") + 
    theme_classic() + theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())


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
    write_csv("swapped_sites_CIrank.csv")






