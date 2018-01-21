########## where actual means NOT RANKS.
library(tidyverse)
library(cowplot)
theme_set(theme_classic())
options(tibble.width = Inf)


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


enzyme <- read_csv("../summarized_data/rate_inferences_enzymes.csv")
virus <- read_csv("../summarized_data/rate_inferences_virus.csv")
org <- read_csv("../summarized_data/rate_inferences_organelles.csv")

models <- c("WAG", "LG", "JTT", "JC69", "gcpREV", "mtVer")
upper.unbounded <- 1000
upper.MLE <- 50

bind_rows(enzyme, virus, org) %>%
    filter(rv == "No") %>%
    filter(model %in% models) %>%
    select(dataset, model, site, MLE, Lower, Upper, type) %>%
    mutate(MLE = ifelse(MLE >= upper.MLE, upper.MLE, MLE)) %>%
    group_by(dataset, model) %>% 
    mutate(medianMLE = median(MLE)) %>%
    ungroup() -> rates

### Total number of sites per dataset ###
rates %>% 
    select(dataset, type, site) %>% 
    group_by(dataset) %>% 
    filter(site==max(site)) %>% 
    unique() %>% 
    rename(sitecount = site) -> dataset.sites

dataset.sites %>%
    ungroup() %>%
    group_by(type) %>%
    summarize(sitecount = sum(sitecount)) -> type.sites
#        
### GENERALLY PRECOMPUTED FOR RUNTIME SANITY, as in mental health. ####
# ## Any sites for which the CIs are NOT FULLY overlapping.
# rates %>%
#     group_by(dataset, site) %>% 
#     nest() %>% 
#     mutate(stuff = map(data, get.models.outside)) %>%
#     select(dataset, site, stuff) %>% 
#     unnest() %>%  #dataset  site main.model diff.models
#     unnest() -> outside.ci
# write_csv(outside.ci, "sites_outside_CI.csv")

outside.ci <- read_csv("sites_outside_ci.csv") 


#######################################################################################################
######### Find percentage of sites per alignment whose rates are in statistical agreement. ###########
######### 1. Any sites where all model CIs fully encompass all model MLEs. No evidence these rates differ from one another.
######### 2. Any sites where all model CIs fully encompass all gene median MLEs. No evidence these rates differ from median gene rate.

### Fully overlapping CI
rates %>% 
    anti_join(outside.ci) -> sites.overlap.ci

sites.overlap.ci %>%
    filter(model == "WAG") %>% ## arbitrary, just to get 1 row per site
    select(-model, -Lower, -Upper, -MLE) %>%
    group_by(dataset, type) %>%
    tally() %>%
    left_join(dataset.sites) -> overlap.ci.count


## Does CI for all models show full overlap with median? If so, we are consistent.
outside.ci %>% 
    select(dataset, site) %>%
    inner_join(rates) %>%
    unique() %>%
    mutate(differs.from.median = (Lower < medianMLE & Upper < medianMLE) | (Lower > medianMLE & Upper > medianMLE)) %>%
    group_by(dataset,type, site) %>%
    summarize(alldiff = sum(differs.from.median)) %>%  ## if sum==0, all sites consistent with median rate
    filter(alldiff==0) %>%
    select(-alldiff) -> sites.overlap.median

sites.overlap.median %>%
    ungroup() %>%
    group_by(dataset, type) %>%
    tally()-> overlap.median.count
    
### Bring together sites with overlapping CI across MLE or median 
overlap.median.count %>%  
    rename(overlap.median = n) %>% 
    left_join(overlap.ci.count) %>% 
    rename(ci.agree = n) %>% mutate(total.agreement = overlap.median + ci.agree) %>%
    mutate(percent.agreement = total.agreement / sitecount) -> fully.agreeing.count

sites.overlap.ci %>% 
    select(dataset, type, site) %>%
    bind_rows(sites.overlap.median) %>%
    unique() -> fully.agreeing.sites

#######################################################################################################
    
    

###########################################################################################################################################
######### Find percentage of sites per alignment with at least one disagreeing model, where disagreement SWAPS SIDES OF MEDIAN. ###########

outside.ci %>% 
    select(dataset, site) %>%
    inner_join(rates) %>%
    unique() %>%
    anti_join(fully.agreeing.sites) %>%
    mutate(side = sign(MLE - medianMLE)) -> disagreeing.rates
    
#### Which sites fall on different sides of the median across models? ###
disagreeing.rates %>%
    ungroup() %>%
    group_by(dataset, site) %>%
    summarize(sameside = abs(sum(side)) == length(models),
              modeside = sign(sum(side))) -> disagreeing.rates.sides
   

### non-ambiguous mode side
disagreeing.rates.sides %>%
    filter(modeside !=0) %>%
    filter(sameside == FALSE) %>%   ### retain where median side inconsistencies within site
    select(-sameside) %>%
    inner_join(disagreeing.rates) %>%
    filter(side != modeside) -> swapped.sites1 ### site/model where side has switched from mode
 ## 3,984 sites.
 ## Example does not inspire interesting-ness:
#  1  13PK_A gcpREV   140 0.9996021 0.7729046 1.268852 enzyme 0.9948173
# 2  13PK_A   JC69   140 1.2948740 1.0028944 1.640258 enzyme 0.9894930
# 3  13PK_A    JTT   140 0.9776193 0.7563245 1.240093 enzyme 0.9784735
# 4  13PK_A     LG   140 0.9967421 0.7710557 1.264270 enzyme 0.9727344
# 5  13PK_A  mtVer   140 0.8904826 0.6877131 1.132035 enzyme 1.0170045
# 6  13PK_A    WAG   140 1.1261013 0.8715962 1.427609 enzyme 0.9586352
# swapped.sites1 %>% mutate(mederror = abs(MLE-medianMLE)/medianMLE) %>% filter(mederror <=5) %>% ggplot(aes(x = mederror)) +geom_density()

### ambiguous mode side
disagreeing.rates.sides %>%
    filter(modeside ==0) %>%
    filter(sameside == FALSE) %>%   ### retain where median side inconsistencies within site
    select(-sameside) %>%
    inner_join(disagreeing.rates) -> swapped.sites2 ### site/model where side has switched from mode

bind_rows(swapped.sites1, swapped.sites2) -> all.swapped

all.swapped %>% 
    ungroup() %>% 
    select(dataset, site) %>% 
    unique() %>% 
    group_by(dataset) %>% 
    tally() -> swapped.count


disagreeing.rates.sides %>%
    filter(sameside == TRUE) %>%   ### retain where median side inconsistencies within site
    select(-sameside) %>%
    inner_join(disagreeing.rates) %>%
    filter(model == "WAG") %>% ## arbitrary, just to get 1 row per site
    group_by(dataset, type) %>% tally() -> disagree.sameside.count



######### THE TRULY DIFFERENT SITES: NEVER A MEDIAN OVERLAP ############
disagreeing.rates.sides %>%
    filter(sameside == FALSE) %>%   ### retain where median side inconsistencies within site
    select(-sameside) %>%
    inner_join(disagreeing.rates) %>%
    mutate(overlapmedian = medianMLE > Lower & medianMLE < Upper) %>%
    ungroup() %>%
    group_by(dataset, site) %>%
    summarize(medover = sum(overlapmedian)) %>% 
    filter(medover==0) -> entirely.different
# > entirely.different %>% inner_join(dataset.sites) %>% ungroup() %>% group_by(type) %>% tally()
# Joining, by = "dataset"
# # A tibble: 2 x 2
#     type     n
#    <chr> <int>
# 1 enzyme   111
# 2   mito     3    
    
# 
# 
# ### total sites of swapped-ness. ###
# all.swapped %>% 
#     ggplot(aes(x = model, fill = model)) + 
#     geom_bar() +
#     theme(legend.position = "none", axis.text.x = element_text(angle=30)) +
#     xlab("Model") + ylab("Count") +
#     scale_x_discrete(limits=c("JC69", "mtVer", "gcpREV", "JTT", "WAG", "LG")) +
#     facet_grid(~factor(type, levels=c("enzyme", "chloro", "mito"))) ->a
# all.swapped.mostest %>% 
#     ggplot(aes(x = model, fill = model)) + 
#     geom_bar() + 
#     theme(legend.position = "none", axis.text.x = element_text(angle=30)) +
#     xlab("Model") + ylab("Count") +
#     scale_x_discrete(limits=c("JC69", "mtVer", "gcpREV", "JTT", "WAG", "LG")) +
#         facet_grid(~factor(type, levels=c("enzyme", "chloro", "mito")))  -> b
# plot_grid(a,b,nrow=2,labels="auto")
# 

##################################### FINARRY #######################################
overlap.ci.count %>%
    rename(ci.overlap = n) %>%
    left_join(overlap.median.count) %>%
    rename(median.overlap = n) %>%
    left_join(disagree.sameside.count) %>%
    rename(nooverlap.sameside = n) %>%
    left_join(swapped.count) %>%
    rename(nooverlap.diffside = n) %>%
    rename(total.sites = sitecount) %>%
    gather(flavor, count, nooverlap.diffside, nooverlap.sameside, ci.overlap, median.overlap) %>%
    replace_na(list(count = 0)) -> total.counts

total.counts %>% 
    mutate(percent = count/total.sites)  %>%
    ungroup() %>%
    group_by(type, flavor) %>% 
    arrange(dataset) %>%
    summarize(meanp = mean(percent) * 100) -> ready
        
ready$type   <- factor(ready$type, levels=c("enzyme", "chloro", "mito", "virus"), labels=c("Enzyme", "Chloroplast", "Mitochondria", "Virus"))
ready$flavor <- factor(ready$flavor, levels=c("ci.overlap", "median.overlap", "nooverlap.sameside", "nooverlap.diffside"), labels=c("Agree, CI", "Agree, Median", "Disagree, all same side", "Disagree, >=1 different side"))


dodge <- position_dodge(width=0.9)
ggplot(ready, aes(x = type, y = meanp, fill=flavor)) +
    geom_bar(stat="identity", position=dodge) +
    geom_text(aes(x = type, y = 2 + meanp, label = round(meanp,2)), size=2.5, position=dodge) + 
    xlab("Dataset") + ylab("Average percentage of sites, per dataset") + 
    theme(legend.title = element_blank())

#######################################################################################

## CI plots with example dataset for types of agreement and disagreement
ds = "1A50_B"

rates %>%
    filter(dataset == ds) %>%
    filter(site %in% c(54, 111, 349, 393)) -> example.sites
    
example.sites$model <- factor(example.sites$model, levels=rev(c("LG", "WAG", "JTT", "mtVer", "gcpREV", "JC69")))
mle.size <- 2.25
median.size <- 2
ci.size <- 1.5
example.sites %>%
    filter(site==54) %>%
    ggplot() + 
    geom_segment(aes(x = Lower, xend = Upper, y = model, yend=model, color = model), size=ci.size) + 
    geom_point(aes(x = MLE, y = model), size=mle.size) + 
    theme(legend.position = "none") + 
    xlab("Rate and 95% CI") + ylab("Model") -> ci
    
example.sites %>%
    filter(site==111) %>%
    ggplot() + 
    geom_segment(aes(x = Lower, xend = Upper, y = model, yend=model, color = model), size=ci.size) + 
    geom_point(aes(x = MLE, y = model), size=mle.size) +  
    geom_point(aes(x = medianMLE, y = model), shape=23, size=median.size, fill="red") +  
    theme(legend.position = "none") + 
    xlab("Rate and 95% CI") + ylab("Model") -> med
    

example.sites %>%
    filter(site==349) %>%
    ggplot() + 
    geom_segment(aes(x = Lower, xend = Upper, y = model, yend=model, color = model), size=ci.size) +  
    geom_point(aes(x = MLE, y = model), size=mle.size) +  
    geom_point(aes(x = medianMLE, y = model), shape=23, size=median.size, fill="red") +  
    theme(legend.position = "none") + 
    xlab("Rate and 95% CI") + ylab("Model") -> diff.sameside ## good since some median overlap

example.sites %>%
    filter(site==393) %>%
    ggplot() + 
    geom_segment(aes(x = Lower, xend = Upper, y = model, yend=model, color = model), size=ci.size) + 
    geom_point(aes(x = MLE, y = model), size=mle.size) +  
    geom_point(aes(x = medianMLE, y = model), shape=23, size=median.size, fill="red") +  
    theme(legend.position = "none") + 
    xlab("Rate and 95% CI") + ylab("Model") -> diff.diffside 


plot_grid(ci, med, diff.sameside, diff.diffside, nrow=2,labels="auto") -> x
save_plot("grid_ci_plot.pdf", x, base_width=7, base_height=4)

example.sites %>%
    ggplot() + 
    geom_segment(aes(x = Lower, xend = Upper, y = model, yend=model, color = model), size=ci.size) + 
    geom_point(aes(x = MLE, y = model), size=mle.size) +  
    geom_point(aes(x = medianMLE, y = model), shape=23, size=median.size, fill="red") +  
    theme(legend.position = "none") + 
    facet_wrap(~site, scales="free_x")+
    xlab("Rate and 95% CI") + ylab("Model") -> facet.ci.plots
save_plot("facet_ci_plot.pdf", facet.ci.plots, base_width=7, base_height=5)










