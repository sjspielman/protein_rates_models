###############################################################################################################################
### SJS
### This script generates the file `rate_inferences_all.csv`, which contains all inferred rates in one easy place.
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
  mutate(type = "enzyme") %>%
  rbind(raw.organelle.rates) -> raw.rates

##### SAVE RATES TO rate_inferences_all.csv
write_csv(raw.rates, "rate_inferences_all.csv")
