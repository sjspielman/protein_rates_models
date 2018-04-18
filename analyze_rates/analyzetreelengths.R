library(tidyverse)
library(cowplot)
figdir <- "figures/"
datadir <- "summarized_data/"
theme_set(theme_classic() + theme(strip.background = element_rect(fill = "grey90")))
model.levels  <- c("LG", "WAG", "JTT", "gcpREV", "mtMet", "HIVBm", "JC69")
model.labels  <- c("LG", "WAG", "JTT", "gcpREV", "mtMet", "HIVB", "JC69")
type.levels <- c("enzyme", "mito", "chloro", "gpcr")
type.labels <- c("Enzyme", "Mitochondria", "Chloroplast", "GPCR")


tl <- read_csv(paste0(datadir, "inferred_tree_lengths.csv"))
tl$model <- factor(tl$model, levels = model.levels)
tl$type <- factor(tl$type, levels = type.levels, labels = type.labels)

### Normalize, for rv each, by LG tree length
tl %>%
    filter(model == "LG") %>%
    rename(LG = treelength) %>%
    select(-model) -> tl2
    
tl %>% 
    left_join(tl2) %>%
    mutate(normtreelength = treelength / LG) -> normtl

# 1 rpl36   HIVBm G     Chloroplast      20144  76.1            265

normtl %>% filter(dataset !="rpl36") -> normtl2


normtl2 %>%
    ggplot(aes(x = model, y = normtreelength, group = dataset)) + 
        geom_point() + geom_line() + 
        geom_hline(yintercept=1, color="red") + 
        facet_wrap(~rv+type, scales="free", nrow=2) + 
        theme(axis.text.x =element_text(angle=30)) +
        background_grid() + panel_border() -> plot.tl
save_plot("inferred_bl.pdf", plot.tl, base_width=11, base_height=4)

tl %>%
    spread(rv, treelength) %>%
    filter(G <= 1000) %>%
    ggplot(aes(x = No, y = G)) + geom_point() + geom_smooth(method ="lm", se=FALSE) + 
        facet_wrap(~model, scales="free", nrow=1) +
        geom_abline(color="red", size=1)  + 
        panel_border() -> rv.tl
save_plot("tl_rv.pdf", rv.tl, base_width=10, base_height=2)
        
        
normtl2 %>% mutate(relerr = abs(LG-treelength)/LG) %>% ggplot(aes(x = model, y= relerr))+geom_jitter(size=1) + facet_wrap(~rv) + scale_y_continuous(limits=c(0,1)) -> blrelerr
save_plot("relerror_tl.pdf", blrelerr)

        
        
        
        
        
        
        

        
        
        
        