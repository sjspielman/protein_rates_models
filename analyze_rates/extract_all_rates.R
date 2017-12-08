###############################################################################################################################
### SJS
### This script generates the two files, `rate_inferences_enzymes.csv` and `rate_inferences_organelles.csv` which together contain all inferred rates.
### Note that two files are used because of GitHub file size constraints.
###############################################################################################################################



require(tidyverse)
require(purrr)
mito <- c("ATP6", "ATP8", "COX1", "COX2", "COX3", "CYTB", "ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6") 


###### Read in organelle rates #######
fpath = "../rate_inference/organelle_data-inference/"
files <- dir(path = fpath, pattern = "*csv")
data_frame(filename = files) %>%
  mutate(file_contents = map(filename,
           ~ read_csv(file.path(fpath, .)))) %>%
  unnest() %>%
  separate(filename, c("dataset", "byebye1", "model", "byebye2", "byebye3"), sep="\\.") %>%
  select(dataset, model, site, MLE, Lower, Upper) %>%
  separate(model, c("model", "rv"), "-") %>%  
  mutate(type = ifelse(dataset %in% mito, "mito", "chloro")) -> raw.organelle.rates
write_csv(raw.organelle.rates, "rate_inferences_organelles.csv")


###### Read in enzyme rates and merge with organelle to create single tibble of rates #######
fpath = "../rate_inference/enzyme_data-inference/"
files <- dir(path = fpath, pattern = "*csv")
data_frame(filename = files) %>%
  mutate(file_contents = map(filename,
           ~ read_csv(file.path(fpath, .)))) %>%
  unnest() %>%
  separate(filename, c("dataset", "byebye1", "model", "byebye2", "byebye3"), sep="\\.") %>%
  select(dataset, model, site, MLE, Lower, Upper) %>% 
  separate(model, c("model", "rv"), "-") %>%  
  mutate(type = "enzyme") -> raw.enzyme.rates
write_csv(raw.enzyme.rates, "rate_inferences_enzyme.csv")

