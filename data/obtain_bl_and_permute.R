####################################################
##
## For all datasets, get mean, median, sum branch length.
## Permute branch lengths as well in say 10 ways each
## 
####################################################
require(ape)
require(tidyverse)
set.seed(101188) # I'm extra creative

datapath <- paste0("trees/")
files <- dir(path = datapath, pattern = "*tre")

permpath <- "permuted_branch_lengths/trees_permuted_bl/"


bldf <- data.frame(dataset = as.character(), meanbl = as.numeric(), medianbl = as.numeric(), treelength = as.numeric())
for (file in files){
    print(file)
    name <- str_split(file, "\\.")[[1]][1]
    tree <- read.tree(paste0(datapath, file))
    
    ##### for branch length csv #####
    temp <- data.frame(dataset = name, 
                       meanbl = mean(tree$edge.length),
                       medianbl = median(tree$edge.length),
                       treelength = sum(tree$edge.length))
    bldf <- rbind(bldf, temp)
    
    
    ##### permute branch lengths 10x each #####
    original.bl <- tree$edge.length
    for (i in 1:10){
        outtree <- paste0(permpath, name, "_permuted_bl_", i, ".tre")
        tree$edge.length <- sample(original.bl)
        write.tree(tree, outtree)
    }
}

write_csv(bldf, "tree_branch_lengths.csv")
        
    
    



