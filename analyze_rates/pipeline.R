### To reproduce analyses, first run these two python scripts in this order:
### 1. `python parse_leisr.py`
### 2. `python parse_inferred_tree_lengths.py`

### Then, execute this R script from **within an active session**

## NOTE: Reproducibility of `detect_differing_rates.R` is plagued by some crazy inefficient code, see this file for more details. 
#### Currently the long chunk is commented out to avoid you going insane and values are precomputed in a csv summarized_data/sites_outside_ci.csv

source("extract_all_rates.R")      #### This data is currently in summarized_data as tarball'd csv's
source("load.R")                   #### Load all the data and set up factors, some plotting settings, etc
source("compare_fits.R")           #### Compare fits among models and make according figures
source("detect_differing_sites.R") #### Assess comparability of site rates inferred with different models and make according figures
source("compare_rates.R")          #### Compare inferred rates among models and make according figures
