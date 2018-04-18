#source("load.R") ### shared plotting bits

get.models.outside <- function(dat){
    sitestuff <- tibble(main.model = character(), diff.models = list())
    models.outside <- list()
    for (m in unique(dat$model))
    {
        total <- 0
        dat %>% filter(model == m) -> mod
        anti_join(dat, mod, by = c("model", "MLE", "Lower", "Upper")) -> notmod
        modmle <- mod$MLE
        notmod %>% filter(mod$MLE < Lower | mod$MLE > Upper) %>% select(model) -> x
        if (nrow(x) == 0){next}
        sitestuff <- bind_rows(sitestuff, tibble(main.model = m, diff.models = list(as.character(x$model))))
    }
    sitestuff 
}

upper.unbounded <- 1000
upper.MLE <- 50

rates %>%  
    filter(rv == "No") %>%
    select(dataset, model, site, MLE, Lower, Upper, type) %>%
    mutate(MLE = ifelse(MLE >= upper.MLE, upper.MLE, MLE)) %>%
    group_by(dataset, model) %>% 
    mutate(medianMLE = median(MLE)) %>%
    ungroup() -> rates2

### Total number of sites per dataset ###
rates2 %>% 
    select(dataset, type, site) %>% 
    group_by(dataset) %>% 
    filter(site==max(site)) %>% 
    unique() %>% 
    rename(sitecount = site) -> dataset.sites

dataset.sites %>%
    ungroup() %>%
    group_by(type) %>%
    summarize(sitecount = sum(sitecount)) -> type.sites


     
### GENERALLY PRECOMPUTED FOR RUNTIME SANITY, where sanity in this case is defined as mental health. ####
### Code should really be run parallelized, but this way, you can take a whole like 8 hours off and feel like you're working. ###
# ## Any sites for which the CIs are NOT FULLY overlapping.
# rates %>%
#     group_by(dataset, site) %>% 
#     nest() %>% 
#     mutate(stuff = map(data, get.models.outside)) %>%
#     select(dataset, site, stuff) %>%
#     unnest() %>%  #dataset  site main.model diff.models
#     unnest() -> outside.ci.neworg
# write_csv(outside.ci, paste0(datadir,"sites_outside_ci.csv"))

outside.ci <- read_csv(paste0(datadir, "sites_outside_ci.csv"))


#######################################################################################################
######### Find percentage of sites per alignment whose rates are in statistical agreement. ###########
######### 1. Any sites where all model CIs fully encompass all model MLEs. No evidence these rates differ from one another.
######### 2. Any sites where all model CIs fully encompass all gene median MLEs. No evidence these rates differ from median gene rate.

### Fully overlapping CI
rates2 %>% 
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
    inner_join(rates2) %>%
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
    inner_join(rates2) %>%
    unique() %>%
    anti_join(fully.agreeing.sites) %>%
    mutate(side = sign(MLE - medianMLE)) -> disagreeing.rates
    
#### Which sites fall on different sides of the median across models? ###
disagreeing.rates %>%
    ungroup() %>%
    group_by(dataset, site) %>%
    summarize(sameside = abs(sum(side)) == length(model.levels),
              modeside = sign(sum(side))) -> disagreeing.rates.sides
   

### non-ambiguous mode side
disagreeing.rates.sides %>%
    filter(modeside !=0) %>%
    filter(sameside == FALSE) %>%   ### retain where median side inconsistencies within site
    select(-sameside) %>%
    inner_join(disagreeing.rates) %>%
    filter(side != modeside) -> swapped.sites1 ### site/model where side has switched from mode
 
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
    spread(flavor, count) %>%
    mutate(agreement = ci.overlap + median.overlap) %>%
    select(-ci.overlap, -median.overlap) %>%
    gather(flavor, count, agreement, nooverlap.diffside, nooverlap.sameside) -> total.counts.agreemerged

total.counts.agreemerged %>% 
    mutate(percent = count/total.sites)  %>%
    ungroup() %>%
    group_by(type, flavor) %>% 
    arrange(dataset) -> percentages

percentages %>% summarize(meanp = mean(percent) * 100) -> ready
ready$flavor <- factor(ready$flavor, levels=c("agreement", "nooverlap.sameside", "nooverlap.diffside"), labels=c("Fully or mostly concordant sites", "Inconsistent sites", "Discordant or \ndiametrically opposed sites"))

# dodge <- position_dodge(width=0.9)
# ggplot(ready, aes(x = type, y = meanp, fill=flavor)) +
#     geom_bar(stat="identity", position=dodge) +
#     geom_text(aes(x = type, y = 2 + meanp, label = round(meanp,2)), size=2.5, position=dodge) + 
#     xlab("Dataset") + ylab("Average percentage of sites, per alignment") + 
#     theme(legend.title = element_blank()) -> percentage.flavors
#save_plot(paste0(figdir,"percentage_types_agreement.pdf"),percentage.flavors, base_width=7) 

dodge <- position_dodge(width=0.9)
ggplot(ready, aes(x = "", y = meanp, fill=flavor)) +
    facet_grid(~type) + 
    geom_bar(stat="identity", position=dodge) +
    geom_text(aes(x = "", y = 5 + meanp, label = round(meanp,2)), size=2.5,position=dodge) + 
    ylab("Percentage of sites") + xlab("") + 
    panel_border() + 
    theme(legend.title = element_blank(), axis.ticks.x = element_blank(), legend.text = element_text(size=8.5)) -> percentage.bars

agreement.legend <- get_legend(percentage.bars)

perc.bl <- read_csv(paste0(datadir,"dataset_properties.csv")) %>% 
    select(dataset, treelength) %>% 
    left_join(percentages, bl, by = "dataset")
perc.bl$flavor <- factor(perc.bl$flavor, levels=c("agreement", "nooverlap.sameside", "nooverlap.diffside"))

ggplot(perc.bl, aes(x = treelength, y = percent*100, color = flavor)) +
    geom_point(alpha = 0.4, size=1) + 
    facet_grid(~type, scales="free") + 
    geom_smooth(method="lm", se=FALSE) + 
    panel_border() + 
    xlab("Tree length") + ylab("Percentage of sites") + 
    theme(legend.position = "none") -> perc.bl.scatter

plot_grid(percentage.bars  + theme(legend.position="none"), perc.bl.scatter, nrow=2, labels="auto") -> g1
plot_grid(g1, agreement.legend, nrow=1, rel_heights=c(1, 0.1), rel_widths=c(1,0.3)) -> percent.grid
save_plot(paste0(figdir, "percentage_agreement_grid.pdf"), percent.grid, base_width=9)

# 
# > fit <- lm(percent ~ treelength*flavor, data=perc.bl)
# > summary(fit)
# Call:
# lm(formula = percent ~ treelength * flavor, data = perc.bl)
# 
# Residuals:
#      Min       1Q   Median       3Q      Max 
# -0.28792 -0.01963 -0.00357  0.02350  0.35223 
# 
# Coefficients:
#                                       Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                          9.377e-01  3.386e-03  276.97   <2e-16 ***
# treelength                          -3.124e-03  7.125e-05  -43.85   <2e-16 ***
# flavornooverlap.sameside            -8.788e-01  4.788e-03 -183.55   <2e-16 ***
# flavornooverlap.diffside            -9.344e-01  4.788e-03 -195.14   <2e-16 ***
# treelength:flavornooverlap.sameside  5.672e-03  1.008e-04   56.29   <2e-16 ***
# treelength:flavornooverlap.diffside  3.701e-03  1.008e-04   36.73   <2e-16 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.0581 on 1251 degrees of freedom
# Multiple R-squared:  0.9778,	Adjusted R-squared:  0.9777 
# F-statistic: 1.1e+04 on 5 and 1251 DF,  p-value: < 2.2e-16



#######################################################################################


###### Which models tend to outlie? ############################
##### ONLY for sites with a clear direction (no ambiguous) #####
disagreeing.rates.sides %>% 
    filter(sameside == FALSE, modeside != 0) %>%
    inner_join(all.swapped) -> outlying.models
    
rates2 %>% 
    mutate(inmedian = (medianMLE >= Lower &medianMLE <= Upper)) %>% 
    group_by(dataset, site) %>% 
    summarize(totalin = sum(inmedian)) %>% 
    filter(totalin == 0)  -> no.median.overlap.atall

## entirely jc69, where JC69 is higher and rest are lower. this is 100% consistent. damn. !!!!!!!!  
outlying.models %>% 
    inner_join(no.median.overlap.atall) %>%
    select(dataset, site, type) %>%
    ungroup() %>%
   inner_join(rates2) -> outlying.sites



### Can be explained by ILV:
process.ilv <- function(df, subdf, m)
{
    df %>%
        ungroup() %>%
        inner_join(subdf) %>%
        filter(subtype !="total", model == m) %>% 
        mutate(ilv = ifelse(subtype %in% c("IL","IV", "LV"), TRUE, FALSE)) %>% 
        select(-subtype) %>% 
        unique() %>%
        group_by(dataset, site, ilv) %>%
        mutate(sumcounts = sum(count)) %>% 
        select(-model,-count) %>% 
        unique() %>% 
        ungroup() %>% 
        spread(ilv, sumcounts) %>% 
        replace_na(list(`FALSE`=0, `TRUE`=0)) %>%
        rename(ilv = `TRUE`, notilv  =`FALSE`) %>%
        mutate(prop.ilv = ilv / (notilv + ilv)) -> ilv.spread
        
    ilv.spread
}

enzyme <- read_csv("summarized_data/enzyme_substitution_counts.csv")
org <- read_csv("summarized_data/organelle_substitution_counts.csv")
gpcr <- read_csv("summarized_data/gpcr_substitution_counts.csv")



full.subcounts <- bind_rows(enzyme, org, gpcr)
full.subcounts$type <- factor(full.subcounts$type, levels = type.levels, labels = type.labels)



######### Do for two models, and place all but one in SI ########
for (m in c("LG", "JC69")){

  
  ######### SAME SIDE ###########
  disagreeing.rates.sides %>%
      filter(sameside == TRUE) %>%   ### retain where median side inconsistencies within site
      select(-sameside) %>%
      inner_join(disagreeing.rates) %>%
      ungroup() %>%
      select(dataset, site, type) %>%
      unique() %>%
      process.ilv(full.subcounts, m) -> ilv.spread.sameside
  
      
  ######## DIFFERENT SIDES ############  
  all.swapped %>%
    anti_join(outlying.sites) %>%
      ungroup() %>%
      select(dataset, site, type) %>%
      unique() %>%
      process.ilv(full.subcounts, m) -> ilv.spread.swapped
  
  
  
  fully.agreeing.sites %>%
      process.ilv(full.subcounts, m) ->ilv.spread.agree
  
          
          
  outlying.sites %>% 
      select(dataset, site, type) %>%
      process.ilv(full.subcounts, m) ->ilv.spread.outlying
  
  
  
  bind_rows("agree" = ilv.spread.agree, 
            "sameside" = ilv.spread.sameside, 
            "diffside" = ilv.spread.swapped, 
            "outlying" = ilv.spread.outlying, 
            .id = "sitetype") %>%
            ungroup() %>%
            mutate(totalsubs = ilv + notilv) %>%
            group_by(sitetype) %>% 
            arrange(prop.ilv) %>%
            mutate(n = 1:n())-> subs
            
  subs$sitetype <- factor(subs$sitetype, levels=c("agree", "sameside", "diffside", "outlying"), 
                           labels=c("Fully or mostly\nconcordant sites", 
                                    "Inconsistent sites", 
                                    "Discordant sites",
                                    "Diametrically\nopposed sites"))
  
  subs %>% 
    group_by(sitetype) %>% 
    summarize(x1 = median(prop.ilv), x2 = median(totalsubs)) -> subs.means
  
  
  subs %>% 
    filter(totalsubs > 0) %>%
    ggplot(aes(x = prop.ilv, color = sitetype, fill = sitetype)) + 
      geom_density(alpha=0.5) + 
    geom_vline(data = subs.means, aes(xintercept = x1), color="grey40") +
    scale_x_continuous(limits=c(0,1)) + 
      scale_fill_manual(values = c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF")) + 
      scale_color_manual(values=c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF"))+        
      facet_wrap(~sitetype, nrow=1,scales="free") + 
      panel_border() + 
      theme(legend.position = "none", strip.text = element_text(size=9)) + 
      ylab("Density") + xlab("Proportion {I,L,V} substitutions per site")  -> ilv.prop.plot
  
  subs %>% 
    filter(prop.ilv >= 0.5) %>%
    ggplot(aes(x = totalsubs, color = sitetype, fill = sitetype)) + 
      geom_density(alpha=0.5) + 
      scale_fill_manual(values = c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF")) + 
      scale_color_manual(values=c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF"))+        
      facet_wrap(~sitetype, nrow=1,scales="free_y") +
      scale_x_continuous(limits=c(0,125)) + 
      geom_vline(data = subs.means, aes(xintercept = x2), color="grey40") +
      panel_border() + 
      theme(legend.position = "none", strip.text = element_text(size=9)) + 
      ylab("Density") + xlab("Total substitution counts") -> subcount.plot
  
  together <- plot_grid(ilv.prop.plot, subcount.plot, nrow=2, labels = "auto", scale=0.98)
  if (m== "JC69"){
    outname = "SI_ilv_plot_JC.pdf"
  }
  else{
    outname = "ilv_plot_LG.pdf"
  }
    save_plot(paste0(figdir, outname), together, base_width=8, base_height=4.25)
}

# > subs%>% group_by(sitetype) %>% summarize(min(prop.ilv))
# # A tibble: 4 x 2
# sitetype                         `min(prop.ilv)`
# <fct>                                      <dbl>
#   1 Fully or mostly concordant sites           0    
# 2 Inconsistent sites                         0    
# 3 Discordant sites                           0    
# 4 Diametrically opposed sites                0.560]





################################################################################################################

## CI plots with example dataset for types of agreement and disagreement
ds = repr[1] ### "1AJ8_A"
rates2 %>%
    filter(dataset == ds) %>%
    filter(site %in% c(467, 196, 51, 158, 87)) -> example.sites


theme_set(theme_classic() + theme(legend.position = "none", 
                                    strip.background = element_rect(fill = "grey90"),
                                    strip.text = element_text(size=9)))
mle.size <- 1.3
median.size <- 1.1
ci.size <- 1
example.sites %>% mutate(site2 = case_when(
                                       site == 467 ~ "Fully concordant site \n(467)",  
                                       site == 196 ~ "Mostly concordant site \n(196)", 
                                       site == 51 ~ "Inconsistent site \n(51)",
                                       site == 158 ~ "Discordant site \n(158)",
                                       site == 87 ~ "Diametrically opposed site \n(87)")) -> example.sites2
                                       
example.sites2$site2 <- factor(example.sites2$site2, levels=c("Fully concordant site \n(467)", "Mostly concordant site \n(196)",  "Inconsistent site \n(51)", "Discordant site \n(158)", "Diametrically opposed site \n(87)"))
example.sites2 %>% mutate(medianMLE = ifelse(site == 467, NA, medianMLE)) -> example.sites2 
example.sites2$model <- factor(example.sites2$model, levels=rev(levels(example.sites2$model)))

ggplot(example.sites2) +
    geom_segment(aes(x = Lower, xend = Upper, y = model, yend=model, color = model), size=ci.size) + 
    geom_point(aes(x = MLE, y = model), size=mle.size) + 
    geom_point(aes(x = medianMLE, y = model), shape=23, size=median.size, fill="red") +  
    facet_wrap(~site2, scales="free_x", nrow=1) +
    xlab("MLE and 95% CI") + ylab("Model") +
    panel_border() -> mleci
save_plot(paste0(figdir,"grid_ci_plot.pdf"), mleci, base_width=9.5, base_height=2.15)

# 
# example.sites %>%
#     filter(site==467) %>%
#     ggplot() + 
#     geom_segment(aes(x = Lower, xend = Upper, y = model, yend=model, color = model), size=ci.size) + 
#     geom_point(aes(x = MLE, y = model), size=mle.size) + 
#     xlab("MLE and 95% CI") + ylab("Model") + ggtitle("Fully concordant site (467)") -> ci
#     
# example.sites %>%
#     filter(site==196) %>%
#     ggplot() + 
#     geom_segment(aes(x = Lower, xend = Upper, y = model, yend=model, color = model), size=ci.size) + 
#     geom_point(aes(x = MLE, y = model), size=mle.size) +  
#     geom_point(aes(x = medianMLE, y = model), shape=23, size=median.size, fill="red") +  
#     xlab("MLE and 95% CI") + ylab("Model") + ggtitle("Generally concordant site (196)") -> med
#     
# 
# example.sites %>%
#     filter(site==51) %>%
#     ggplot() + 
#     geom_segment(aes(x = Lower, xend = Upper, y = model, yend=model, color = model), size=ci.size) +  
#     geom_point(aes(x = MLE, y = model), size=mle.size) +  
#     geom_point(aes(x = medianMLE, y = model), shape=23, size=median.size, fill="red") +  
#     xlab("MLE and 95% CI") + ylab("Model") + ggtitle("Inconsistent site (51)")-> diff.sameside 
# 
# example.sites %>%
#     filter(site==158) %>%
#     ggplot() + 
#     geom_segment(aes(x = Lower, xend = Upper, y = model, yend=model, color = model), size=ci.size) + 
#     geom_point(aes(x = MLE, y = model), size=mle.size) +  
#     geom_point(aes(x = medianMLE, y = model), shape=23, size=median.size, fill="red") +  
#     xlab("MLE and 95% CI") + ylab("Model") + ggtitle("Discordant site (158)")-> diff.diffside1 
# 
# example.sites %>%
#     filter(site==87) %>%
#     ggplot() + 
#     geom_segment(aes(x = Lower, xend = Upper, y = model, yend=model, color = model), size=ci.size) + 
#     geom_point(aes(x = MLE, y = model), size=mle.size) +  
#     geom_point(aes(x = medianMLE, y = model), shape=23, size=median.size, fill="red") +  
#     xlab("MLE and 95% CI") + ylab("Model") + ggtitle("Diametrically opposite site (87)") -> diff.diffside2
# plot_grid(ci, med, diff.sameside, diff.diffside1, diff.diffside2, nrow=1,labels="auto", label_size=11) -> x
# save_plot(paste0(figdir,"grid_ci_plot.pdf"), x, base_width=10, base_height=1.8)
# 
# 
