### Rate analysis pipeline

All LEISR rate inferences are in `../rate_inferences/<organelle or enzyme>_data-inference/`.

Run in this order:

**bold** 
*emph*

```
	inline code chunk
```
[link](link)





+ The python script `parse_leisr.py` extracts a CSV for each JSON (sends to same directory as JSON) and creates the file `model_fits.csv`, which has per dataset and inference:
	+ LogL, AICc, number of sequences, number of sites, model, datatype (enzyme, chloro, mito), gamma parameter (if Gamma model)

+ The R script `extract_all_rates.R` creates the file `rate_inferences_all.csv`, which contains ALL rate inferences in one place for downstream convenience

+ The R script `compare_fits.R` determines best-fit model per dataset from the file `model_fits.csv`. All mitochondrial are best fit by mtVer, chloroplasts are 50/50 cpREV and <other>, and enzymes are largely WAG with some cpREV and JTT. This is fairly consistent with these models as land plants were only a small subset of cpREV training data.

+ The R script `compare_rates.R` calculates correlations across models and produces figures

> this is a block quote

