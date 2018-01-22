######### Shared libraries, settings to load ##########
require(cowplot)
require(tidyverse)
require(broom)
options(tibble.width = Inf)

mito <- c("ATP6", "ATP8", "COX1", "COX2", "COX3", "CYTB", "ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6") 
figdir <- "figures/"
datadir <- "summarized_data/"
theme_set(theme_classic() + theme(strip.background = element_rect(fill = "grey90")))
all.models    <- c("LG", "WAG", "JTT", "gcpREV",  "mtVer", "JC69", "permWAG", "RAND")
model.levels  <- c("LG", "WAG", "JTT", "gcpREV",  "mtVer", "JC69")
model.levels2 <- c("LG", "WAG", "permWAG", "RAND")
type.levels <- c("enzyme", "chloro", "mito", "virus")
type.labels <- c("Enzyme", "Chloroplast", "Mitochondria", "Virus")
repr <- c("1A50_B", "matK", "CYTB", "Porcine_reproductive_and_respiratory_syndrome_virus_E_AA") ## Representative datasets for scatterplots
repr.labels <- c("Tryp Synthase", "matK", "CYTB", "PRRSV")
