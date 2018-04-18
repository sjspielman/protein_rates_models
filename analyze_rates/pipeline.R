## Run in an active R session ##
source("extract_all_rates.R")

source("load.R")
source("compare_fits.R")
source("detect_differing_sites.R") ### Note that reproducibility of this step is plagued by some crazy inefficient code, see this file for more details.
source("compare_rates.R")
