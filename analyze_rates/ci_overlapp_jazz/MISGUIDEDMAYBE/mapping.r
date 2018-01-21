library(tidyverse)

get.outside.count <- function(dat){
    sitestuff <- tibble(model = character(), count.outside = numeric())
    for (m in dat$model)
    {
        total <- 0
        dat %>% filter(model == m) -> mod
        anti_join(dat, mod, by = c("model", "MLE", "Lower", "Upper")) -> notmod
        modmle <- mod$MLE
        total <- sum(modmle < notmod$Lower | modmle > notmod$Upper)
        sitestuff <- bind_rows(sitestuff, tibble(model = m, count.outside=total))
    }
    sitestuff %>% nest()
}



get.models.outside <- function(dat){
    sitestuff <- tibble(main.model = character(), diff.models = list())
    models.outside <- list()
    for (m in dat$model)
    {
        total <- 0
        dat %>% filter(model == m) -> mod
        anti_join(dat, mod, by = c("model", "MLE", "Lower", "Upper")) -> notmod
        modmle <- mod$MLE
        notmod %>% filter(mod$MLE < Lower | mod$MLE > Upper) %>% select(model) -> x
        sitestuff <- bind_rows(sitestuff, tibble(main.model = m, diff.models = list(x$model)))
    }
    sitestuff #%>% nest()
}



enzyme <- read_csv("summarized_data/rate_inferences_enzymes.csv")
virus <- read_csv("summarized_data/rate_inferences_virus.csv")
org <- read_csv("summarized_data/rate_inferences_organelles.csv")
bind_rows(enzyme, virus, org) %>%
    filter(rv == "No") %>%
    filter(model %in% c("WAG", "LG", "JTT", "JC69", "gcpREV", "mtVer")) %>%
    select(-LogL_global, -LogL_local, -rv) -> rates 


rates %>%
    group_by(dataset, site) %>% 
    nest() %>% 
    mutate(stuff = map(data, get.outside.count)) %>%
    select(dataset, site, stuff) %>% 
    unnest() %>% 
    unnest() %>%
    write_csv("outside_ci.csv")
    

boring <- read_csv("outplease/sites_consistentCI.csv")
anti_join(rates, boring) %>% 
    arrange(dataset, site)  -> rates2
    
rates2 %>%
    group_by(dataset, site) %>% 
    nest() %>% 
    mutate(stuff = map(data, get.models.outside)) %>%
    select(dataset, site, stuff) %>% 
    unnest() %>%  #dataset  site main.model diff.models
    unnest() %>%
    write_csv("whichmodelsoutsideCI.csv")


    