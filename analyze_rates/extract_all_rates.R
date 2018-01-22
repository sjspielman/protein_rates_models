###############################################################################################################################
### SJS
### This script generates the two files, `rate_inferences_enzymes.csv` and `rate_inferences_organelles.csv` which together contain all inferred rates.
### Note that two files are used, rather than one, because of GitHub file size constraints.
###############################################################################################################################



require(tidyverse)
mito <- c("ATP6", "ATP8", "COX1", "COX2", "COX3", "CYTB", "ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6") 

outpath <- "summarized_data/"

###### Read in and save virus rates #######
fpath = "../rate_inference/virus_data-inference/"
files <- dir(path = fpath, pattern = "*csv")
data_frame(filename = files) %>%
  mutate(file_contents = map(filename,
           ~ read_csv(file.path(fpath, .)))) %>%
  unnest() %>%
  separate(filename, c("dataset", "byebye1", "model", "byebye2", "byebye3"), sep="\\.") %>%
  select(dataset, model, site, MLE, Lower, Upper, LogL_global, LogL_local) %>%
  separate(model, c("model", "rv"), "-") %>%
  mutate(type = "virus") -> raw.virus.rates
write_csv(raw.virus.rates, paste0(outpath, "rate_inferences_virus.csv"))


###### Read in and save organelle rates #######
fpath = "../rate_inference/organelle_data-inference/"
files <- dir(path = fpath, pattern = "*csv")
data_frame(filename = files) %>%
  mutate(file_contents = map(filename,
           ~ read_csv(file.path(fpath, .)))) %>%
  unnest() %>%
  separate(filename, c("dataset", "byebye1", "model", "byebye2", "byebye3"), sep="\\.") %>%
  select(dataset, model, site, MLE, Lower, Upper, LogL_global, LogL_local) %>%
  separate(model, c("model", "rv"), "-") %>%  
  mutate(type = ifelse(dataset %in% mito, "mito", "chloro")) -> raw.organelle.rates
write_csv(raw.organelle.rates, paste0(outpath,"rate_inferences_organelles.csv"))


###### Read in and save enzyme rates #######
fpath = "../rate_inference/enzyme_data-inference/"
files <- dir(path = fpath, pattern = "*csv")
data_frame(filename = files) %>%
  mutate(file_contents = map(filename,
           ~ read_csv(file.path(fpath, .)))) %>%
  unnest() %>%
  separate(filename, c("dataset", "byebye1", "model", "byebye2", "byebye3"), sep="\\.") %>%
  select(dataset, model, site, MLE) %>%
  separate(model, c("model", "constraint"), "-") %>%  
  mutate(type = "enzyme") -> raw.enzyme.rates
write_csv(raw.enzyme.rates, paste0(outpath,"rate_inferences_enzymes.csv"))



###### Read in and save dataset-specific rates #######
fpath = "../rate_inference/custom_matrix-inference/"
files <- dir(path = fpath, pattern = "*csv")
data_frame(filename = files) %>%
  mutate(file_contents = map(filename,
           ~ read_csv(file.path(fpath, .)))) %>%
  unnest() %>%
  separate(filename, c("dataset", "byebye1", "model", "byebye2", "byebye3"), sep="\\.") %>%
  select(dataset, model, site, MLE, Lower, Upper, LogL_global, LogL_local) %>%
  mutate(type = "custom") -> raw.custom.rates
write_csv(raw.custom.rates, paste0(outpath,"rate_inferences_custom.csv"))


