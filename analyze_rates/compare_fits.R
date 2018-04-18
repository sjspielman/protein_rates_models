#source("load.R")


model.fits <- read_csv(paste0(datadir, "model_fits.csv"))
model.fits %>% 
    mutate(AIC = 2*(nparam - logl), 
           BIC = (log(sequences*sites)*nparam) - 2*logl) %>%
    gather(criterion, value, AIC, AICc, BIC) -> model.fits2

model.fits2 %>% 
    group_by(dataset, criterion) %>%
    mutate(critrank = as.integer(rank(value))) -> model.ranks    
    
    
    
model.ranks %>%
    filter(critrank == 1) %>%
    select(dataset, sequences, sites, criterion, value) %>%
    rename(best = value) -> best.model

model.fits2 %>% left_join(best.model) -> compare.models
compare.models %>%
    ungroup() %>% 
    group_by(dataset, type, criterion) %>%
    mutate(weight = exp(-0.5*(value - best))) %>%
    ungroup() %>%
    group_by(dataset, criterion) %>% 
    mutate(normweight =weight/sum(weight)) %>%
    mutate(critrank = as.integer(rank(value))) %>%
    arrange(dataset) -> all.weights
all.weights$type <- factor(all.weights$type, levels=type.levels, labels = type.labels)
all.weights %>% filter(value == best, normweight < 0.5) %>% ungroup() %>% group_by(criterion, dataset, type) %>% tally()
# A tibble: 3 x 4
# Groups:   criterion, dataset [?]
#   criterion dataset type       n
#   <chr>     <chr>   <fct>  <int>
# 1 AIC       PSBI    Chloro     1
# 2 AICc      PSBI    Chloro     1
# 3 BIC       PSBI    Chloro     1



all.weights %>% 
    filter(critrank == 1, criterion == "AIC") %>% 
    ggplot(aes(x = type, y = normweight)) + 
    geom_jitter(aes(fill = type), pch=21) + 
    xlab("Dataset") + ylab("Weight of Best Model") +
    geom_hline(yintercept=0.5, color="grey50") +
    scale_fill_manual(values=ordered.colors) + 
    theme(legend.position = "none", axis.title.y = element_text(size=10))  -> weights.fig

model.fits2 %>% 
    separate(model, c("model", "rv"), sep = "-") %>%
    filter(rv == "No", criterion == "AIC") %>%
    group_by(dataset) %>%
    mutate(critrank = as.integer(rank(value))) %>%
    filter(critrank == 1) %>% 
    mutate(model = ifelse(model %in% c("JTT", "WAG", "LG"), "LG/WAG/JTT", model)) %>%
    group_by(type, model) %>% 
    tally() %>% 
    ungroup() %>%
    group_by(type) %>%
    mutate(p = n/sum(n), rv = "No rate variation") ->x.norv
model.fits2 %>% 
    separate(model, c("model", "rv"), sep = "-") %>%
    filter(rv == "G", criterion == "AIC") %>%
    group_by(dataset) %>%
    mutate(critrank = as.integer(rank(value))) %>%
    filter(critrank == 1) %>% 
    mutate(model = ifelse(model %in% c("JTT", "WAG", "LG"), "LG/WAG/JTT", model)) %>%
    group_by(type, model) %>% 
    tally() %>% 
    ungroup() %>%
    group_by(type) %>%
    mutate(p = n/sum(n), rv = "+G rate variation") %>%
    bind_rows(x.norv) -> x
x$model <- factor(x$model, levels=c("LG/WAG/JTT", "mtMet", "gcpREV", "HIVBm", "JC69"), labels=c("Gen.", "mtMet", "gcpREV", "HIVB", "JC69"))
x$type <- factor(x$type, levels=c("enzyme", "chloro", "mito", "gpcr"), labels=c("Enzyme", "Chloroplast", "Mitochondria", "GPCR"))
x$rv <- factor(x$rv, levels=c("No rate variation", "+G rate variation"))
ggplot(x, aes(fill = model, x=type, y=p)) +
  geom_bar(stat="identity") + 
  xlab("Dataset") + ylab("Percent Preferred") + 
  scale_fill_manual(name = "Model", values = ordered.colors) + 
  facet_wrap(~rv, ncol=1) + 
  theme(legend.position = "bottom", 
        legend.title = element_text(size=10), 
        legend.key.size = unit(0.4, "cm"), 
        legend.text = element_text(size = 7))  -> pref.models
        
plot_grid(pref.models, weights.fig, labels="auto", scale=0.98) -> modelfits
save_plot(paste0(figdir, "model_fits_weights.pdf"), modelfits, base_width=7.3, base_height=4)



# 
# all.weights %>% 
#     select(dataset, model, type, criterion, critrank) %>%
#     filter(critrank == 1) %>%
#     separate(model, c("model", "rv")) %>%
#     select(-model) %>%
#     group_by(type, criterion, rv) %>% 
#     tally() %>% 
#     ungroup() %>%
#     group_by(type, criterion) %>%
#     mutate(p = n/sum(n)) ->x 
# x$rv <- factor(x$rv, levels=c("G", "No"), labels=c("Gamma", "None"))
# #x$type <- factor(x$type, levels=c("enzyme", "chloro", "mito", "virus"), , labels=c("Enzyme", "Chloroplast", "Mitochondria", "Virus"))
# ggplot(x, aes(fill = rv, x=type, y=p)) +
#   geom_bar(stat="identity") + 
#   xlab("Dataset") + ylab("Percent Preferred") + 
#   scale_fill_discrete(name = "Rate variation") +  
#   facet_wrap(~criterion, ncol=1) + 
#   theme(legend.position = "bottom", 
#         legend.title = element_text(size=10), 
#         legend.key.size = unit(0.4, "cm"), 
#         legend.text = element_text(size = 9)) -> rv.plot
# pref.models <- plot_grid(model.plot, rv.plot, nrow=2, labels="auto", label_size=12, scale=0.95)
# save_plot(paste0(figdir, "preferred_models.pdf"), pref.models, base_width=5.25, base_height=4)
# pref.models <- plot_grid(model.plot, rv.plot, nrow=1, labels="auto", scale=0.98)
# save_plot(paste0(figdir, "pref_models.pdf"), pref.models, base_width=6.25, base_height=5.75)



##################################################################################################################
##################################################################################################################
##################################################################################################################

tl <- model.fits %>% 
    select(dataset, model, type, treelength) %>%
    separate(model, c("model", "rv"), sep="-")
tl$model <- factor(tl$model, levels = model.levels, labels = model.labels)
tl$type <- factor(tl$type, levels = type.levels, labels = type.labels)

### Normalize, for rv each, by LG tree length
tl %>%
    filter(model == "LG") %>%
    rename(LG = treelength) %>%
    select(-model) -> tl2
    
tl %>% 
    left_join(tl2) %>%
    mutate(normtreelength = treelength / LG) -> normtl

### Right censor RPS12
tl %>% mutate(treelength = ifelse(treelength > 500, 500, treelength)) -> tl
normtl %>% mutate(normtreelength = ifelse(normtreelength > 10, 10, normtreelength)) -> normtl


normtl %>%
    filter(dataset != "RPS12") %>% 
    ggplot(aes(x = model, y = normtreelength, group = dataset)) + 
        geom_point() + geom_line() + 
        geom_hline(yintercept=1, color="red") + 
        facet_wrap(~rv+type, scales="free", nrow=2) + 
        theme(axis.text.x =element_text(angle=30, size=8)) +
        background_grid() + panel_border() +
        xlab("Model") + ylab("Normalized tree length") -> plot.tl
save_plot(paste0(figdir,"SI_inferred_treelength.pdf"), plot.tl, base_width=11, base_height=4)

tl %>%
    spread(rv, treelength) %>%
    ggplot(aes(x = No, y = G)) + geom_point() + geom_smooth(method ="lm", se=FALSE) + 
        facet_wrap(~model, scales="free", nrow=1) +
        geom_abline(color="red", size=1)  + 
        panel_border() +
        xlab("Tree length, no rate variation") + ylab("Tree length, +G") -> rv.tl
save_plot(paste0(figdir,"SI_treelength_ratevariation.pdf"), rv.tl, base_width=10, base_height=2)
        

####### This is probably overkill ##########
#         
# normtl %>% 
#     mutate(relerr = abs(LG-treelength)/LG) %>% 
#     filter(model != "LG") %>% 
#     ggplot(aes(x = model, y= relerr)) + 
#         geom_jitter(size=1) + 
#         facet_wrap(~rv) + 
#         scale_y_log10(breaks=c(0.001, 0.01, 0.1, 1.0, 10., 100.))+
#         xlab("Model") + ylab("Relative error, compared to LG") -> blrelerr
# save_plot(paste0(figdir,"SI_treelength_relerror.pdf"), base_width=8, base_height=3, blrelerr)
#         
#         
        
        
        
        
        

        
        
        
        
