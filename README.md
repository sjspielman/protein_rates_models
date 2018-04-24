Repository for [*Relative evolutionary rates in proteins are largely insensitive to the substitution model*](https://doi.org/10.1101/304758)

+ [Stephanie J. Spielman](http://sjspielman.org) 
	+ \* Corresponding: `stephanie.spielman <at> temple.edu`
+ [Sergei L. Kosakovsky Pond](http://hyphy.org)

-----------------------------------------------------

### Contents

+ `rate_inference/` contains all code and results for inferring relative evolutionary rates using the [LEISR](https://peerj.com/articles/4339/) implementation

+ `data/` contains all sequence alignments and phylogenies, including generating code, analyzed in `rate_inference/`

+ `count_subs/` contains all code and results for substitution counting

+ `analyze_rates/` contains all code used to analyze and process data, as well as MS figures and processed data CSVs
	+ Reproduce using the R script `analyze_rates/pipeline.R`. See comments in this file for some key notes, including Python scripts which need to be run first.