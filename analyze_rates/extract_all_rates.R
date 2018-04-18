###############################################################################################################################
### SJS
### This script generates the three files, `rate_inferences_<enzyme,organelle,gpcr.csv` which together contain all inferred rates.
### Note that three files are used, rather than one, because of GitHub file size constraints.
###############################################################################################################################



require(tidyverse)
mito <- c("ATP6", "ATP8", "COX1", "COX2", "COX3", "CYTB", "ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6") 

outpath <- "summarized_data/"

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
write_csv(raw.organelle.rates, paste0(outpath,"rate_inferences_organelle.csv"))

###### Read in and save enzyme rates #######
fpath = "../rate_inference/enzyme_data-inference/"
files <- dir(path = fpath, pattern = "*csv")
data_frame(filename = files) %>%
  mutate(file_contents = map(filename,
           ~ read_csv(file.path(fpath, .)))) %>%
  unnest() %>%
  separate(filename, c("dataset", "byebye1", "model", "byebye2", "byebye3"), sep="\\.") %>%
  select(dataset, model, site, MLE, Lower, Upper, LogL_global, LogL_local) %>%
  separate(model, c("model", "rv"), "-") %>%  
  mutate(type = "enzyme") -> raw.enzyme.rates
write_csv(raw.enzyme.rates, paste0(outpath,"rate_inferences_enzyme.csv"))



###### Read in and save gpcr rates #######
fpath = "../rate_inference/gpcr_data-inference/"
files <- dir(path = fpath, pattern = "*csv")
data_frame(filename = files) %>%
  mutate(file_contents = map(filename,
           ~ read_csv(file.path(fpath, .)))) %>%
  unnest() %>%
  separate(filename, c("dataset", "byebye1", "model", "byebye2", "byebye3"), sep="\\.") %>%
  select(dataset, model, site, MLE, Lower, Upper, LogL_global, LogL_local) %>%
  separate(model, c("model", "rv"), "-") %>%  
  mutate(type = "gpcr") -> raw.gpcr.rates
write_csv(raw.gpcr.rates, paste0(outpath,"rate_inferences_gpcr.csv"))

