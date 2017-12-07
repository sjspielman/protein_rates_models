###### Crux analysis file for rates
###### SJS
##### docu forthcoming.....

require(cowplot)
require(tidyverse)
require(broom)
mito <- c("ATP6", "ATP8", "COX1", "COX2", "COX3", "CYTB", "ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6") 


raw.rates <- read_csv("rate_inferences_all.csv")


### Plotting things ###
theme_set(theme_classic() + theme(strip.background = element_rect(fill = "grey90")))
comparison.levels <- c("LG-WAG", "LG-JTT", "LG-mtVer", "LG-cpREV", "WAG-JTT", "WAG-mtVer", "WAG-cpREV", "JTT-mtVer", "JTT-cpREV", "mtVer-cpREV", "LG-JC69", "WAG-JC69", "JTT-JC69", "mtVer-JC69", "cpREV-JC69")
model.levels <- c("LG", "WAG", "JTT", "mtVer", "gcpREV", "JC69")
type.levels <- c("enzyme", "chloro", "mito")
type.labels <- c("Enzyme", "Chloroplast", "Mitochondria")






### Obtain a dataframe of rates only, performing the following modifications, right-censoring >=1e3 --> 1e3 ###
raw.rates %>% 
    group_by(dataset, model, rv) %>%
    select(-Lower,-Upper) %>%
    mutate(MLE  = ifelse(MLE >= 1e3, 1e3, MLE)) -> rates



##########################################################################################
############################# Representative scatterplots ################################
##########################################################################################

repr <- c("1A50_B", "matK", "CYTB")
rates$model <- factor(rates$model, levels = model.levels)
rates$type <- factor(rates$type, levels = type.levels, labels = type.labels)

rates %>% 
    filter(dataset %in% repr) %>%
    group_by(model, type) %>%
    filter(MLE >= 1e-8) %>%
    spread(rv, MLE) %>%
    ggplot(aes(x = No, y = G)) + 
        geom_point(size=1) + 
        geom_abline(color = "blue") + 
        facet_grid(type~model) + 
        panel_border() + 
        scale_x_log10() + scale_y_log10() + 
        xlab("MLE, no rate variation") + ylab("MLE, Gamma") -> scatter.grid.rv

ggsave("scatterplot_grid_ratevariation.pdf", scatter.grid.rv, width=9, height=3.5)


rates %>%
    ungroup() %>%
    filter(rv == "G", dataset %in% repr, MLE >= 1e-8) %>%
    select(-rv) %>%
    group_by(type) %>%
    spread(model, MLE) %>%
    gather(model, MLE, WAG:JC69) -> lg.vs.all
lg.vs.all$model <- factor(lg.vs.all$model, levels = model.levels[c(-1)])
#lg.vs.all$type <- factor(lg.vs.all$type, levels = type.levels, labels = type.labels)


lg.vs.all %>%
    ggplot(aes(x = LG, y = MLE)) + 
        geom_point(size=1) + 
        geom_abline(color = "blue") + 
        facet_grid(type~model) + 
        panel_border() + 
        scale_x_log10() + scale_y_log10() + 
        xlab("LG MLE") + ylab("Model MLE") -> scatter.grid.lg.vs.all
ggsave("scatterplot_grid_LG_vs_all.pdf", scatter.grid.lg.vs.all, width=9, height=3.5)

#########################################################################################


rates %>%
    ungroup() %>% 
    group_by(dataset, type, model) %>%
    spread(rv, MLE) %>%
    summarize( rho = cor(G, No, method = "spearman")) -> rv.correlations


####### PLOT ######
#rv.correlations$model <- factor(rv.correlations$model, levels=model.levels)
#rv.correlations$type <- factor(rv.correlations$type, levels=type.levels, labels=type.labels)
rv.correlations %>% 
    ggplot(aes(x = model, y = rho, fill = type))+
    geom_point(pch=21, position = position_jitterdodge())  + 
    ylab("Rank Correlation") + xlab("Model") + scale_fill_discrete(name = "Dataset") -> rv.jitter.plot
### Outliers are rpl36 (WAG, JTT, JC), and psbL (JC)
ggsave("jitter_ratevariation.pdf", rv.jitter.plot, width = 6, height=2)


 
##### SAVE NoRV-Gamma correlations to correlations_rate_variation.csv
write_csv(rv.correlations, "correlations_rate_variation.csv")



##########################################################################################
##########################################################################################
##########################################################################################



run.correlation.models <- function(dat, model1, model2) {
    if (model1 == model2) {
        rhos <- data.frame(dataset = dat$dataset,
                           type = dat$type,
                           rho = 1,
                           "model1" = model1, 
                           "model2" = model2)
    }
    else{
        dat %>% 
            select(dataset, type, model1, model2) %>%
            ungroup() -> subdat
        names(subdat) <- c("dataset", "type", "m1", "m2")
        subdat %>% 
            group_by(type, dataset) %>%
            summarize( rho = cor(m1, m2, method = "spearman")) %>%
            mutate("model1" = model1, "model2" = model2) -> rhos
    }
    as.data.frame(rhos)
}

##### Rank correlations between each pair of models, for gamma only
rates %>%
    ungroup() %>%
    filter(rv == "G") %>%
    select(-rv) %>%
    group_by(dataset, type) %>% 
    spread(model, MLE) -> spread.gamma

corr.between.models <- data.frame("dataset" = character(), "type" = character(), "rho" = numeric(), "model1" = character(), "model2" = character())
for (m1 in model.levels){
    for (m2 in model.levels) {
            thisrho <- run.correlation.models(spread.gamma, m1, m2)
            corr.between.models <- rbind( corr.between.models, thisrho) 
    }
}

write_csv(corr.between.models, "correlation_between_models.csv")

##########################################################################################
############################ Upper triangle heatmap #####################################
extract.corr.type <- function(dat, mytype){
    meanrho %>% 
        ungroup() %>%
        filter(type == mytype) %>%
        select(-type) %>%
        spread(model2, x) %>% 
        select(-model1) %>% as.matrix() -> w
    rownames(w) <- model.levels
    w[upper.tri(w, diag=FALSE)] <- NA
    as.tibble(w) %>% mutate(type = mytype) -> w
    w
}

#### mean ####
corr.between.models %>%
    ungroup() %>%
    group_by(type, model1, model2) %>%
    summarize(x = mean(rho)) -> meanrho
ecorr <- extract.corr.type(meanrho, "Enzyme")
ccorr <- extract.corr.type(meanrho, "Chloroplast")
mcorr <- extract.corr.type(meanrho, "Mitochondria")
ecorr %>% 
    rbind(ccorr) %>% 
    rbind(mcorr) %>%
    mutate(model2 = rep(model.levels, 3)) %>%
    gather(model1, x, LG:JC69) %>%
    rename(meanrho = x) -> all.meanrho

### standard deviation ###
corr.between.models %>%
    ungroup() %>%
    group_by(type, model1, model2) %>%
    summarize(x = sd(rho)) -> sdrho
ecorr <- extract.corr.type(sdrho, "Enzyme")
ccorr <- extract.corr.type(sdrho, "Chloroplast")
mcorr <- extract.corr.type(sdrho, "Mitochondria")
ecorr %>% 
    rbind(ccorr) %>% 
    rbind(mcorr) %>%
    mutate(model2 = rep(model.levels, 3)) %>%
    gather(model1, x, LG:JC69) %>%
    rename(sdrho = x) -> all.sdrho
    
### join up mean and sd and format the label ###
left_join(all.meanrho, all.sdrho) %>%
    mutate(label = paste0( round(meanrho, 2), "\n(",  round(sdrho, 3), ")") ) %>%
    mutate(label = ifelse(label == "NA\n(NA)", "", label)) -> all.corr
    
    
### factor order
all.corr$model1 <- factor(all.corr$model1, levels=model.levels)
all.corr$model2 <- factor(all.corr$model2, levels=model.levels)
all.corr$type <- factor(all.corr$type, levels=c("Enzyme", "Chloroplast", "Mitochondria"))


### PLOT!
all.corr %>%
    ggplot(aes(x = model1, y = model2, fill = meanrho)) + 
    geom_tile(color="white") + geom_text(aes(label = label), size=2.25) +
    scale_fill_gradient(name = "Rank correlation", low = "red", high = "yellow", na.value = "grey90") + 
    facet_grid(~type) + 
    xlab("Model 1") + ylab("Model 2") +
    theme(axis.text.x = element_text(angle=15)) -> heatmap.rho.between.models
ggsave("between_heatmap.pdf", heatmap.rho.between.models, width=10, height=3)    



##########################################################################################

###################### OLD PLOTTING BELOW IN EVENT OF RESUSCITATION #######################
    
# Code commented out here will make a symmetric heatmap.
# 
# ### heatmap
# corr.between.models$model1 <- factor(corr.between.models$model1, levels=model.levels)
# corr.between.models$model2 <- factor(corr.between.models$model2, levels=model.levels)
# corr.between.models$type <- factor(corr.between.models$type, levels=type.levels, labels=type.labels)
# corr.between.models %>%
#     ungroup() %>%
#     group_by(type, model1, model2) %>%
#     summarize(meanrho = mean(rho)) %>%
#     ggplot(aes(x = model1, y = model2, fill = meanrho)) + 
#     geom_tile(color="white") + geom_text(aes(label = round(meanrho, 2)), size=2.5) +
#     scale_fill_gradient(name = "Mean correlation", low = "red", high = "yellow") + facet_grid(~type) + 
#     xlab("Model 1") + ylab("Model 2")
#     
#     
#     
#     
#     
    
# 
# Code commented out below is an older jitter plot for correlations between models. Now using heatmap in a valiant attempt to diversify figures
# ########## Jitter plot of correlations between models ##########
# corr.models$comparison <- factor(corr.models$comparison, levels=comparison.levels)
# corr.models$type <- factor(corr.models$type, levels=type.levels, labels=type.labels)
# corr.models$correlation.type <- factor(corr.models$correlation.type, levels=corrtype.levels, labels=corrtype.labels)
# corr.models %>% 
#     ggplot(aes(x = comparison, y = r, fill = type))+
#     geom_point(pch=21, position = position_jitterdodge()) + 
#     facet_wrap(~correlation.type, nrow=2) + theme(axis.text.x = element_text(angle=10, face = "bold", size=7), strip.background = element_rect(fill = "grey90")) + 
#     ylab("Correlation") + xlab("Compared models") + scale_fill_discrete(name = "Dataset") -> between.jitter.plot
# #    guides(fill = guide_legend(override.aes = list(size = 2))) + 
# ggsave("between_models_jitter.pdf", between.jitter.plot, height = 3)
