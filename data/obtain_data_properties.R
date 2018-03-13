library(tidyverse)
library(seqinr)
library(ape)

datasets <- read_csv("datasets.csv")

outframe <- tibble(dataset = as.character(),
                   alnlen  = as.numeric(),
                   nseq    = as.numeric(),
                   treelength = as.numeric())

for (d in datasets$dataset) { 

    aln    <- read.fasta(paste0("fasta/", d, ".fasta"))
    alnlen <- length(aln[[names(aln)[1]]])
    nseq   <- length(aln)
    
    tree <- read.tree(paste0("trees/", d, ".tre"))
    tl <- sum(tree$edge.length)
    
    temp <- tibble("dataset" = d,
                   "alnlen"  = alnlen,
                   "nseq"    = nseq,
                   "treelength" = tl)
    outframe <- bind_rows(outframe, temp)
}

outframe <- left_join(outframe, datasets)
write_csv(outframe,"dataset_properties.csv")

##### paper table ########
outframe %>% 
    group_by(type) %>%
    summarize(n(), mean(alnlen), mean(nseq), mean(treelength))
# 1 chloro    79            483        333                 9.29
# 2 enzyme   100            657        248                85.8 
# 3 gpcr     227            455         22.5               2.08
# 4 mito      13            475        342                82.0 
    


