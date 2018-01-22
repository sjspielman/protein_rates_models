###### Crux analysis file for rates
###### SJS
##### docu forthcoming.....
source("load.R") ### shared plotting bits

############################################################################################
################################## Functions ###############################################
############################################################################################

############################################################################################
############################################################################################




############################################################################################
########################### Read in rates and factor as needed  ############################
############################################################################################

### Obtain a dataframe of rates only, performing the following modifications, right-censoring >=1e3 --> 1e3 ###
raw.rates1 <- read_csv(paste0(datadir, "rate_inferences_enzymes.csv"))
raw.rates2 <- read_csv(paste0(datadir,"rate_inferences_organelles.csv"))
raw.rates3 <- read_csv(paste0(datadir,"rate_inferences_virus.csv"))
bind_rows(raw.rates1, raw.rates2, raw.rates3) %>% 
    group_by(dataset, model, rv) %>%
    mutate(MLE  = ifelse(MLE >= 1e3, 1e3, MLE)) %>%
    ungroup() -> rates
    
    
    
### Excluding RAND and permWAG ###
rates %>% filter(model %in% model.levels) -> rates.real ####### exclude RAND and permWAG
rates.real$model <- factor(rates.real$model, levels = model.levels)
rates.real$type <- factor(rates.real$type, levels = type.levels, labels = type.labels)


### With LG, WAG, RAND, permWAG ###
rates %>% filter(model %in% model.levels2) -> rates.fake
rates.fake$model <- factor(rates.fake$model, levels = model.levels2)
rates.fake$type <- factor(rates.fake$type, levels = type.levels, labels = type.labels)

############################################################################################
############################################################################################




############################################################################################
###################### Create, save correlations between models ############################
############################################################################################


#### No RV - Gamma. 
rates.real %>%
    select(-Lower, -Upper, -LogL_global, -LogL_local) %>% 
    group_by(dataset, type, model) %>%
    spread(rv, MLE) %>%
    summarize( rho = cor(G, No, method = "spearman")) -> rv.correlations
write_csv(rv.correlations, paste0(datadir, "correlations_between_rv.csv"))


## Spread rates
rates.real %>%
    filter(rv == "No") %>%
    select(-rv, -Lower, -Upper, -LogL_global, -LogL_local) %>%
    group_by(dataset, type) %>% 
    spread(model, MLE) %>%
    ungroup() -> spread.rates

## Correlate pairwise
corr.between.models <- data.frame("dataset" = character(), "type" = character(), "rho" = numeric(), "model1" = character(), "model2" = character())
for (m1 in model.levels){
    for (m2 in model.levels) {
    
        if (m1 == m2) {
            rhos <- data.frame(dataset = spread.rates$dataset,
                           type = spread.rates$type,
                           rho = 1,
                           "model1" = m1, 
                           "model2" = m2)
        } else{
            spread.rates %>% 
                select(dataset, type, m1, m2) %>%
                ungroup() -> subdat
            names(subdat) <- c("dataset", "type", "m1", "m2")
            subdat %>% 
                group_by(type, dataset) %>%
                summarize( rho = cor(m1, m2, method = "spearman")) %>%
                mutate("model1" = m1, "model2" = m2) %>% 
                as.data.frame() -> rhos
        }
        corr.between.models <- rbind(corr.between.models, rhos) 
    }
}
corr.between.models <- as.tibble(corr.between.models)
write_csv(corr.between.models, paste0(datadir, "correlations_between_models.csv"))
############################################################################################
############################################################################################




##########################################################################################
############################# Representative scatterplots ################################
##########################################################################################

rates.real %>% 
    filter(dataset %in% repr, MLE >= 1e-8) %>%
    group_by(model, type) %>%
    select(-Lower, -Upper, -LogL_global, -LogL_local) %>% 
    spread(rv, MLE) %>%
    ggplot(aes(x = No, y = G)) + 
        geom_point(size=1) + 
        geom_abline(color = "blue") + 
        facet_grid(type~model) + 
        panel_border() + 
        scale_x_log10(labels = c(0.1, 1., 10., 100.)) + 
        scale_y_log10(labels = c(0.1, 1., 10., 100.)) + 
        xlab("MLE, no rate variation") + ylab("MLE, Gamma") -> scatter.grid.rv
ggsave(paste0(figdir, "scatterplot_grid_ratevariation.pdf"), scatter.grid.rv, width=9, height=4)


####### "regular" ######
rates.real %>%
    filter(rv == "No", dataset %in% repr, MLE >= 1e-8) %>%
    select(-rv, -Lower, -Upper, -LogL_global, -LogL_local) %>% 
    group_by(type) %>%
    spread(model, MLE) %>%
    gather(model, MLE, WAG:JC69) -> lg.vs.all
lg.vs.all$model <- factor(lg.vs.all$model, levels = model.levels[c(-1)])


lg.vs.all %>%
    ggplot(aes(x = LG, y = MLE)) + 
        geom_point(size=1) + 
        geom_abline(color = "blue") + 
        facet_grid(type~model) + 
        panel_border() + 
        scale_x_log10(labels = c(0.1, 1., 10., 100.)) + 
        scale_y_log10(labels = c(0.1, 1., 10., 100.)) + 
        xlab("MLE, LG") + ylab("MLE, compared model") -> scatter.grid.lg.vs.all
ggsave(paste0(figdir, "scatterplot_grid_LG_vs_all.pdf"), scatter.grid.lg.vs.all, width=9, height=4)
############################################################################################
############################################################################################



##########################################################################################
###################### Jitter plot: Correlations between no rv/gamma #####################
##########################################################################################

rv.correlations %>% 
    ggplot(aes(x = model, y = rho, fill = type))+
    geom_point(pch=21, position = position_jitterdodge())  + 
    ylab("Rank Correlation") + xlab("Model") + scale_fill_discrete(name = "Dataset")  -> rv.jitter.plot
ggsave(paste0(figdir, "jitter_ratevariation.pdf"), rv.jitter.plot, width = 6, height=3)

############################################################################################
############################################################################################





##########################################################################################
###################### Heatmap:  Correlations between models #############################
##########################################################################################

### Used to create a correlation matrix from dataframe of correlations. Hands down, most efficient solution. ;) ###
extract.corr.type <- function(dat, mytype, rnames){
    dat %>% 
        ungroup() %>%
        filter(type == mytype) %>%
        select(-type) %>%
        spread(model2, x) %>% 
        select(-model1) %>% as.matrix() -> w
    rownames(w) <- rnames
    w[upper.tri(w, diag=FALSE)] <- NA
    as.tibble(w) %>% mutate(type = mytype) -> w
    w
}


##### Heatmap ######

#### mean ####
corr.between.models %>%
    ungroup() %>%
    group_by(type, model1, model2) %>%
    summarize(x = mean(rho)) -> meanrho
ecorr <- extract.corr.type(meanrho, "Enzyme", model.levels)
ccorr <- extract.corr.type(meanrho, "Chloroplast", model.levels)
mcorr <- extract.corr.type(meanrho, "Mitochondria", model.levels)
vcorr <- extract.corr.type(meanrho, "Virus", model.levels)
ecorr %>% 
    rbind(ccorr) %>% 
    rbind(mcorr) %>% 
    rbind(vcorr) %>%
    mutate(model2 = rep(model.levels,4)) %>%
    gather(model1, x, LG:JC69) %>%
    rename(meanrho = x) -> all.meanrho

### standard deviation ###
corr.between.models %>%
    ungroup() %>%
    group_by(type, model1, model2) %>%
    summarize(x = sd(rho)) -> sdrho
ecorr <- extract.corr.type(sdrho, "Enzyme", model.levels)
ccorr <- extract.corr.type(sdrho, "Chloroplast", model.levels)
mcorr <- extract.corr.type(sdrho, "Mitochondria", model.levels)
vcorr <- extract.corr.type(sdrho, "Virus", model.levels)
ecorr %>% 
    rbind(ccorr) %>% 
    rbind(mcorr) %>%
    rbind(vcorr) %>%
    mutate(model2 = rep(model.levels, 4)) %>%
    gather(model1, x, LG:JC69) %>%
    rename(sdrho = x) -> all.sdrho
    
### join up mean and sd and format the label ###
left_join(all.meanrho, all.sdrho) %>%
    mutate(label = paste0( round(meanrho, 2), "\n(",  round(sdrho, 3), ")") ) %>%
    mutate(label = ifelse(label == "NA\n(NA)", "", label)) -> all.corr
    
### factor order
all.corr$model1 <- factor(all.corr$model1, levels=model.levels)
all.corr$model2 <- factor(all.corr$model2, levels=model.levels)
all.corr$type <- factor(all.corr$type, levels=c("Enzyme", "Chloroplast", "Mitochondria", "Virus"))


### Plot, save the heatmap
all.corr %>%
    ggplot(aes(x = model1, y = model2, fill = meanrho)) + 
    geom_tile(color="white") + 
    geom_text(aes(label = label), size=2) + 
    scale_fill_gradient(name = "Rank correlation", low = "red", high = "yellow", na.value = "grey90") + 
    xlab("Model 1") + ylab("Model 2") +
    facet_grid(~type) + 
    theme(axis.text.x = element_text(angle=15, size=8)) -> heatmap.rho.between.models
ggsave(paste0(figdir, "heatmap_correlations.pdf"), heatmap.rho.between.models, width=11, height=3)    


############################################################################################
############################################################################################


##########################################################################################################
########## in conclusion, rates decoupled from matrix. models contribute nothing? #########
########## entropy against rates. high entropy portend high rates as you would expect but is there significant deviation??? 

#you move the fit by sadding rates but by same amount per model
#model fit stems entirely from branch klength optimization and not rates themselves
######## PLOT LOGL ########
rates.real %>%  
    group_by(dataset, type, model, rv) %>% 
    summarize(global = sum(LogL_global), local = sum(LogL_local)) %>% 
    mutate(logl_diff = global-local) -> data.logl.plot
data.logl.plot$rv <- factor(data.logl.plot$rv, levels=c("G", "No"), labels=c("Gamma", "No rate variation"))
data.logl.plot %>% 
    ungroup() %>%
    select(-dataset) %>%
    rename(Dataset = type) %>%
    ggplot(aes(x = model, y = logl_diff, fill = Dataset)) + 
        geom_boxplot(outlier.size = 0.5, size=0.25) +
        facet_wrap(~rv, nrow=2, scales="free_y") +  
        xlab("Model") + ylab("LogL Global - LogL Local") +
        background_grid() -> logl_global_local
save_plot(paste0(figdir, "logl_global_local.pdf"), logl_global_local, base_width=7, base_height=4)




# 
# 
# 
# ##########################################################################################
# ############################### RAND and permWAG analysis ################################
# ##########################################################################################
# 
# 
# rates.fake %>%
#     ungroup() %>%
#     filter(rv == "No", dataset %in% repr, MLE >= 1e-8) %>%
#     select(-rv, -Upper, -Lower, -LogL_global, -LogL_local) %>%
#     group_by(type) %>%
#     spread(model, MLE) %>%
#     gather(model, MLE, WAG:RAND) -> wag.vs.all2
# wag.vs.all2$model <- factor(wag.vs.all2$model, levels = model.levels2[c(-1)])
# 
# 
# wag.vs.all2 %>%
#     ggplot(aes(x = LG, y = MLE)) + 
#         geom_point(size=1) + 
#         geom_abline(color = "blue") + 
#         facet_grid(type~model) + 
#         panel_border() + 
#         scale_x_log10() + scale_y_log10() + 
#         xlab("WAG MLE") + ylab("Model MLE") -> scatter.grid.wag.vs.all
# ggsave(paste0(figdir, "scatterplot_grid_permWAG_RAND.pdf"), scatter.grid.wag.vs.all, width=7, height=3.5)
# 
# 
# 
# 
# ##### Rank correlations between each pair of models, for no rv only
# 
# rates.fake %>%
#     ungroup() %>%
#     filter(rv == "No") %>%
#     select(-rv, -Upper, -Lower, -LogL_global, -LogL_local) %>%
#     group_by(dataset, type) %>% 
#     spread(model, MLE) -> spread.rates2
# 
# corr.between.models2 <- data.frame("dataset" = character(), "type" = character(), "rho" = numeric(), "model1" = character(), "model2" = character())
# for (m1 in model.levels2){
#     for (m2 in model.levels2) {
#     
#         if (m1 == m2) {
#             rhos <- data.frame(dataset = spread.rates2$dataset,
#                            type = spread.rates2$type,
#                            rho = 1,
#                            "model1" = m1, 
#                            "model2" = m2)
#         } else{
#             spread.rates2 %>% 
#                 select(dataset, type, m1, m2) %>%
#                 ungroup() -> subdat
#             names(subdat) <- c("dataset", "type", "m1", "m2")
#             subdat %>% 
#                 group_by(type, dataset) %>%
#                 summarize( rho = cor(m1, m2, method = "spearman")) %>%
#                 mutate("model1" = m1, "model2" = m2) %>% 
#                 as.data.frame() -> rhos
#         }
#       corr.between.models2 <- rbind( corr.between.models2, rhos) 
#     }
# }
# 
# 
# 
# #### mean ####
# corr.between.models2 %>%
#     ungroup() %>%
#     group_by(type, model1, model2) %>%
#     summarize(x = mean(rho)) -> meanrho
# ecorr <- extract.corr.type(meanrho, "Enzyme", model.levels2)
# ccorr <- extract.corr.type(meanrho, "Chloroplast", model.levels2)
# mcorr <- extract.corr.type(meanrho, "Mitochondria", model.levels2)
# ecorr %>% 
#     rbind(ccorr) %>% 
#     rbind(mcorr) %>% 
#     mutate(model2 = rep(model.levels2,3)) %>%
#     gather(model1, x, WAG:RAND) %>%
#     rename(meanrho = x) -> all.meanrho
# 
# ### standard deviation ###
# corr.between.models2 %>%
#     ungroup() %>%
#     group_by(type, model1, model2) %>%
#     summarize(x = sd(rho)) -> sdrho
# ecorr <- extract.corr.type(sdrho, "Enzyme", model.levels2)
# ccorr <- extract.corr.type(sdrho, "Chloroplast", model.levels2)
# mcorr <- extract.corr.type(sdrho, "Mitochondria", model.levels2)
# ecorr %>% 
#     rbind(ccorr) %>% 
#     rbind(mcorr) %>%
#     mutate(model2 = rep(model.levels2, 3)) %>%
#     gather(model1, x, WAG:RAND) %>%
#     rename(sdrho = x) -> all.sdrho
#     
# ### join up mean and sd and format the label ###
# left_join(all.meanrho, all.sdrho) %>%
#     mutate(label = paste0( round(meanrho, 2), "\n(",  round(sdrho, 3), ")") ) %>%
#     mutate(label = ifelse(label == "NA\n(NA)", "", label)) -> all.corr2
#     
#     
# ### factor order
# all.corr2$model1 <- factor(all.corr$model1, levels=model.levels2)
# all.corr2$model2 <- factor(all.corr$model2, levels=model.levels2)
# all.corr2$type <- factor(all.corr$type, levels=c("Enzyme", "Chloroplast", "Mitochondria"))
# 
# 
# ### PLOT!
# all.corr2 %>%
#     ggplot(aes(x = model1, y = model2, fill = meanrho)) + 
#     geom_tile(color="white") + 
#     geom_text(aes(label = label), size=2) + 
#     scale_fill_gradient(name = "Rank correlation", low = "red", high = "yellow", na.value = "grey90") + 
#     xlab("Model 1") + ylab("Model 2") +
#     facet_grid(~type)  -> heatmap.rho.between.models
# ggsave(paste0(figdir, "heatmap_correlations_permWAG_RAND.pdf"), heatmap.rho.between.models, width=10, height=3)    
# 
# 
# 
# 
# ############## Custom matrix analysis ############
# 
# custom.rates <- read_csv(paste0(datadir, "rate_inferences_custom.csv"))
# custom.rates %>% 
#     mutate(model = "custom") %>%
#     mutate(MLE = ifelse(MLE>=1e3, 1e3, MLE)) -> custom.rates
#     
#     
#     
# custom.datasets <- unique(custom.rates$dataset)
# rates %>% 
#     filter(model %in% c("WAG", "LG"), rv == "No", dataset %in% custom.datasets) %>%
#     select(-rv) %>%
#     rbind(custom.rates) %>%
#     filter(MLE >= 1e-8) %>%
#     select(-Lower, -Upper, -LogL_global, -LogL_local, -type) %>%
#     group_by(dataset) %>%
#     spread(model, MLE) %>%
#     gather(model, MLE, LG:WAG) -> custom.rates2
# custom.rates2$model <- factor(custom.rates2$model, levels=c("LG", "WAG"))
# 
# 
# custom.rates2 %>% 
#     group_by(dataset, model) %>% 
#     summarize(rho = cor(custom, MLE, method="spearman")) -> custom.rho
# #   dataset  model `cor(custom, MLE, method = "spearman")`
# #     <chr> <fctr>                                   <dbl>
# # 1  1QAZ_A     LG                               0.9662724
# # 2  1QAZ_A    WAG                               0.9716017
# # 3    ATP8     LG                               0.8619204
# # 4    ATP8    WAG                               0.8413403
# # 5   rpoC1     LG                               0.9858021
# # 6   rpoC1    WAG                               0.9774610
# #        
#     
# custom.rates2 %>%
#     ggplot(aes(x = custom, y = MLE)) + 
#         geom_point(size=1) + 
#         geom_abline(color = "blue") + 
#         facet_grid(dataset~model) + 
#         scale_x_log10() + scale_y_log10() + 
#         panel_border() + 
#         xlab("Custom MLE") + ylab("Model MLE") +
#         geom_text(data = custom.rho, aes(x = 0.035, y = 200, label=paste("rho ==",round(rho, 2))), size=3, parse=T)
# 
# 




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
