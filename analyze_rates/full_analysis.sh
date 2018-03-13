#! /bin/bash

 
###### Convert each JSON to a CSV (in directory `../rate_inference/<type>_data-inference` and create `summarized_data/model_fits.csv` #####
python parse_leisr.py


##### All the R bits, documentary forthcoming.
Rscript pipeline.R