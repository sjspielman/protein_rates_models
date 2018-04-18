######### Shared libraries, settings, data to load ##########
require(cowplot)
require(tidyverse)
require(broom)
options(tibble.width = Inf)

mito <- c("ATP6", "ATP8", "COX1", "COX2", "COX3", "CYTB", "ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6") 
figdir <- "figures/"
datadir <- "summarized_data/"
theme_set(theme_classic() + theme(strip.background = element_rect(fill = "grey90")))

model.levels  <- c("LG", "WAG", "JTT", "gcpREV", "mtMet", "HIVBm", "JC69")
model.labels  <- c("LG", "WAG", "JTT", "gcpREV", "mtMet", "HIVb", "JC")
type.levels <- c("enzyme", "mito", "chloro", "gpcr")
type.labels <- c("Enzyme", "Mitochondria", "Chloroplast", "GPCR")
type.labels.abbr <- c("Enzyme", "Mito", "Chloro", "GPCR")
repr <- c("1AJ8_A", "ENSG00000196639", "MATK", "CYTB") ###, "Porcine_reproductive_and_respiratory_syndrome_virus_E_AA") ## Representative datasets for scatterplots
repr.names <- c("CS", "HRH1", "matK", "CYTB")

# These are the default colors but reordered because I want mito to be pink and chloro to be green, because of science.
## enz, mito, chloro, gpcr
ordered.colors <- c("#00BFC4", "#F8766D", "#7CAE00", "#C77CFF")

raw.rates1 <- read_csv(paste0(datadir, "rate_inferences_enzyme.csv"))
raw.rates2 <- read_csv(paste0(datadir,"rate_inferences_organelle.csv"))
raw.rates3 <- read_csv(paste0(datadir,"rate_inferences_gpcr.csv"))
bind_rows(raw.rates1, raw.rates2, raw.rates3) %>% 
    group_by(dataset, model, rv) %>%
    mutate(MLE  = ifelse(MLE >= 1e3, 1e3, MLE)) %>%
    ungroup() -> rates
rates$model <- factor(rates$model, levels = model.levels, labels = model.labels)
rates$type <- factor(rates$type, levels = type.levels, labels = type.labels)