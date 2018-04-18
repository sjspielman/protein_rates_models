###### Crux analysis file for rates
###### SJS
#source("load.R") ### shared plotting bits

lower <- 1e-8
boundrange <- 1e3

############################################################################################
###################### Create, save correlations between models ############################
############################################################################################

#### No RV - Gamma. 
rates %>%
    select(-Lower, -Upper, -LogL_global, -LogL_local) %>% 
    group_by(dataset, type, model) %>%
    spread(rv, MLE) %>%
    summarize( rho = cor(G, No, method = "spearman")) -> rv.correlations
write_csv(rv.correlations, paste0(datadir, "correlations_between_rv.csv"))


## Spread rates
rates %>%
    filter(rv == "No") %>%
    select(-rv, -Lower, -Upper, -LogL_global, -LogL_local) %>%
    group_by(dataset, type) %>% 
    spread(model, MLE) %>%
    ungroup() -> spread.rates

## Correlate pairwise
corr.between.models <- tibble("dataset" = character(), "type" = character(), "rho" = numeric(), "model1" = character(), "model2" = character())
for (m1 in model.labels){
    for (m2 in model.labels) {
    
        if (m1 == m2) {
            rhos <- data.frame(dataset = spread.rates$dataset,
                           type = spread.rates$type,
                           rho = 1.0,
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
        corr.between.models <- bind_rows(corr.between.models, rhos) 
    }
}
write_csv(corr.between.models, paste0(datadir, "correlations_between_models.csv"))
############################################################################################
############################################################################################




##########################################################################################
############################# Representative scatterplots ################################
##########################################################################################

rates %>% 
    filter(dataset %in% repr, MLE >= lower) %>%
    mutate(bounded = (Upper-Lower)<boundrange) %>%
    group_by(model, type) %>%
    select(-Lower, -Upper, -LogL_global, -LogL_local) %>% 
    spread(rv, MLE) -> rv.scatter.data
rv.scatter.data$type <- factor(rv.scatter.data$type, levels=type.labels, labels=repr.names)
rv.scatter.data %>%
    ggplot(aes(x = No, y = G)) + 
        geom_point(size=1, alpha=0.6, aes(color=bounded)) + 
        geom_abline(color = "blue") + 
        facet_grid(type~model) + 
        panel_border() + 
        scale_color_manual(values = c("grey60", "black")) + 
        scale_x_log10(labels = c(0.1, 1., 10., 100.)) + 
        scale_y_log10(labels = c(0.1, 1., 10., 100.)) + 
        xlab("MLE, no rate variation") + ylab("MLE, +G") +
        theme(legend.position="none", strip.text.y = element_text(size = 8)) -> scatter.grid.rv
ggsave(paste0(figdir, "scatterplot_grid_ratevariation_BOUNDED.pdf"), scatter.grid.rv, width=9, height=4)

rates %>% 
    filter(dataset %in% repr, MLE >= lower, rv == "No") %>%
    mutate(bounded = (Upper-Lower)<boundrange) %>%
    select(-Lower, -Upper, -LogL_global, -LogL_local) %>% 
    group_by(type) %>%
    spread(model, MLE) %>%
    gather(model, MLE, WAG:JC) -> model.scatter.data
model.scatter.data$type <- factor(model.scatter.data$type, levels=type.labels, labels=repr.names)    
model.scatter.data$model <- factor(model.scatter.data$model, levels = model.labels[c(-1)])

model.scatter.data %>% 
    ggplot(aes(x = LG, y = MLE)) + 
        geom_point(size=1, alpha=0.6, aes(color=bounded)) + 
        geom_abline(color = "blue") + 
        facet_grid(type~model) + 
        panel_border() + 
        scale_color_manual(values = c("grey60", "black")) + 
        scale_x_log10(labels = c(0.1, 1., 10., 100.)) + 
        scale_y_log10(labels = c(0.1, 1., 10., 100.)) + 
        xlab("MLE, LG") + ylab("MLE, compared model")  + 
        theme(legend.position="none", strip.text.y = element_text(size = 8))-> scatter.grid.model
ggsave(paste0(figdir, "scatterplot_grid_LG_vs_all_BOUNDED.pdf"),scatter.grid.model, width=9, height=4)
    
############################################################################################
############################################################################################



##########################################################################################
###################### Jitter plot: Correlations between no rv/gamma #####################
##########################################################################################

rv.correlations %>% 
    ggplot(aes(x = model, y = rho, fill = type))+
    geom_point(pch=21, position = position_jitterdodge())  + 
    ylab("Rank Correlation") + xlab("Model") + scale_fill_manual(name = "Dataset", values = ordered.colors) + 
    theme(legend.position = "bottom",legend.box.spacing = unit(0., "cm"), legend.box.margin = margin(0,0,0,0)) -> rv.jitter.plot
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
corr.between.models$model1 <- factor(corr.between.models$model1, levels=model.labels)
corr.between.models$model2 <- factor(corr.between.models$model2, levels=model.labels)
#### mean ####
corr.between.models %>%
    ungroup() %>%
    group_by(type, model1, model2) %>%
    summarize(x = mean(rho)) -> meanrho
ecorr <- extract.corr.type(meanrho, "Enzyme", model.labels)
ccorr <- extract.corr.type(meanrho, "Chloroplast", model.labels)
mcorr <- extract.corr.type(meanrho, "Mitochondria", model.labels)
gcorr <- extract.corr.type(meanrho, "GPCR", model.labels)
bind_rows(ecorr,ccorr, mcorr, gcorr) %>% 
    mutate(model2 = rep(model.labels,4)) %>%
    gather(model1, x, model.labels) %>%
    rename(meanrho = x) -> all.meanrho

### standard deviation ###
corr.between.models %>%
    ungroup() %>%
    group_by(type, model1, model2) %>%
    summarize(x = sd(rho)) -> sdrho
ecorr <- extract.corr.type(sdrho, "Enzyme", model.labels)
ccorr <- extract.corr.type(sdrho, "Chloroplast", model.labels)
mcorr <- extract.corr.type(sdrho, "Mitochondria", model.labels)
gcorr <- extract.corr.type(sdrho, "GPCR", model.labels)
bind_rows(ecorr,ccorr, mcorr, gcorr) %>% 
    mutate(model2 = rep(model.labels, 4)) %>%
    gather(model1, x, model.labels) %>%
    rename(sdrho = x) -> all.sdrho
    
### join up mean and sd and format the label ###
left_join(all.meanrho, all.sdrho) %>%
    mutate(label = paste0( round(meanrho, 2), "\n(",  round(sdrho, 3), ")") ) %>%
    mutate(label = ifelse(label == "NA\n(NA)", "", label)) -> all.corr
    
### factor order
all.corr$model1 <- factor(all.corr$model1, levels=model.labels)
all.corr$model2 <- factor(all.corr$model2, levels=model.labels)
all.corr$type <- factor(all.corr$type, levels=type.labels)


### Plot, save the heatmap
all.corr %>%
    ggplot(aes(x = model1, y = model2, fill = meanrho)) + 
    geom_tile(color="white") + 
    geom_text(aes(label = label), size=2.15) + 
    scale_fill_gradient(name = "Rank Correlation", low = "red", high = "yellow", na.value = "#f2f2f2") + 
    facet_grid(~type) + 
    xlab("") + ylab("") + 
    theme(axis.text.x = element_text(angle=15, size=8), axis.text.y = element_text(size=8.5)) -> heatmap.rho.between.models
ggsave(paste0(figdir, "heatmap_correlations.pdf"), heatmap.rho.between.models, width=11.5, height=2.8)    


############################################################################################
############################################################################################









