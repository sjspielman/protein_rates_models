############
## Pilot round 2
###########
require(tidyverse)
require(cowplot)
theme_set(theme_classic() + theme(strip.background = element_rect(fill = "grey90")))

pilot.datasets <- c("13PK_A", "1C0K_A", "1E6E_A", "1GP1_A", "1MQW_A", "1REQ_A", "1YON_A", "1A16_A", "1C2T_A", "1E7Q_A")
plot.path <- "pilot_randomization_plots/"

fpath = "rate_inference_random_trees/" #1YON_A_random_tree_9.fna.LG.LEISR.json
files <- dir(path = fpath, pattern = "*csv")
data_frame(filename = files) %>%
  mutate(file_contents = map(filename,
           ~ read_csv(file.path(fpath, .)))) %>%
  unnest() %>%
  separate(filename, c("stuff", "byebye1", "model", "byebye2", "byebye3"), sep="\\.", convert=TRUE) %>%
  select(stuff, model, site, MLE) %>%
  separate(stuff, c("data1", "data2", "byebye1", "byebye2", "x"), sep="_", convert=TRUE) %>%
  unite(dataset, c("data1", "data2")) %>%
  select(-byebye1, -byebye2) %>%
  mutate(type = "random_trees") -> random.tree.rates


######################### Forthcoming ##########################
# fpath = "rate_inference_permuted_branch_lengths/"
# files <- dir(path = fpath, pattern = "*csv")
# data_frame(filename = files) %>%
#   mutate(file_contents = map(filename,
#            ~ read_csv(file.path(fpath, .)))) %>%
#   unnest() %>%
#   separate(filename, c("stuff", "byebye1", "model", "byebye2", "byebye3"), sep="\\.", convert=TRUE) %>%
#   select(stuff, model, site, MLE) %>%
#   separate(stuff, c("data1", "data2", "byebye", "byebye2", "x"), sep="_", convert=TRUE) %>%
#   unite(dataset, c("data1", "data2")) %>%
#   select(-byebye, -byebye2) %>%
#   mutate(type = "permuted_bl") -> permuted.bl.rates

######################### Forthcoming ##########################
# fpath = "rate_inference_random_EFV/"
# files <- dir(path = fpath, pattern = "*csv")
# data_frame(filename = files) %>%
#   mutate(file_contents = map(filename,
#            ~ read_csv(file.path(fpath, .)))) %>%
#   unnest() %>%
#   separate(filename, c("stuff", "byebye1", "model", "byebye2", "byebye3"), sep="\\.", convert=TRUE) %>%
#   select(stuff, model, site, MLE) %>%
#   separate(stuff, c("data1", "data2", "byebye", "byebye2", "x"), sep="_", convert=TRUE) %>%
#   unite(dataset, c("data1", "data2")) %>%
#   select(-byebye, -byebye2) %>%
#   mutate(type = "randomEFV") -> permuted.bl.rates



fpath = "rate_inference_branch_length_constrain/meanbl_constraint_inference/" #1EXP_A.fna-meanbl.JC69.LEISR.json
files <- dir(path = fpath, pattern = "*csv")
data_frame(filename = files) %>%
  mutate(file_contents = map(filename,
           ~ read_csv(file.path(fpath, .)))) %>%
  unnest() %>%
  separate(filename, c("dataset", "byebye1", "model", "byebye2", "byebye3"), sep="\\.", convert=TRUE) %>%
  select(dataset, model, site, MLE) %>%
  filter(dataset %in% pilot.datasets) %>%
  mutate(type = "meanbl") -> meanbl.rates


### this might need some replicates ###
fpath = "rate_inference_branch_length_constrain/value_constraint_inference/" #1EXP_A.fna-<gamma,one,five>.JC69.LEISR.json
files <- dir(path = fpath, pattern = "*-gamma.LEISR.csv")
data_frame(filename = files) %>%
  mutate(file_contents = map(filename,
           ~ read_csv(file.path(fpath, .)))) %>%
  unnest() %>%
  separate(filename, c("dataset", "byebye1", "model", "byebye2", "byebye3"), sep="\\.", convert=TRUE) %>%
  separate(model, c("model", "byebye"), sep = "-") %>%
  select(dataset, model, site, MLE) %>%
  filter(dataset %in% pilot.datasets) %>%
  mutate(type = "gammabl") -> gammabl.rates



fpath = "rate_inference_randomized_amino_acids/" #1A16_A_randomized_aa_7.fna.JC69.LEISR.json
files <- dir(path = fpath, pattern = "*csv")
data_frame(filename = files) %>%
  mutate(file_contents = map(filename,
           ~ read_csv(file.path(fpath, .)))) %>%
  unnest() %>%
  separate(filename, c("stuff", "byebye1", "model", "byebye2", "byebye3"), sep="\\.", convert=TRUE) %>%
  select(stuff, model, site, MLE) %>%
  separate(stuff, c("data1", "data2", "byebye", "byebye2", "x"), sep="_", convert=TRUE) %>%
  unite(dataset, c("data1", "data2")) %>%
  select(-byebye, -byebye2) %>%
  mutate(type = "randomizedAA") -> randomized.AA.rates



fpath = "rate_inference_randomized_sequences_over_tree/" #1E6E_A_randomized_sequences_1.fna.LG.LEISR.json
files <- dir(path = fpath, pattern = "*csv")
data_frame(filename = files) %>%
  mutate(file_contents = map(filename,
           ~ read_csv(file.path(fpath, .)))) %>%
  unnest() %>%
  separate(filename, c("stuff", "byebye1", "model", "byebye2", "byebye3"), sep="\\.", convert=TRUE) %>%
  select(stuff, model, site, MLE) %>%
  separate(stuff, c("data1", "data2", "byebye", "byebye2", "x"), sep="_", convert=TRUE) %>%
  unite(dataset, c("data1", "data2")) %>%
  select(-byebye, -byebye2) %>%
  mutate(type = "randomized_sequences") -> randomized.sequences.rates



rates <- read_csv("analyze_rates/summarized_data/rate_inferences_enzymes.csv")
rates %>% 
    filter(dataset %in% pilot.datasets, model %in% c("LG", "WAG", "JC69"), rv == "No") %>% 
    select(dataset, model, site, MLE) %>%
    mutate(type = "real") -> real.rates

#### There must be a better way to do the following (repeat 10 times)
real.rates %>% mutate(x = 1) -> real.rates.x
for (i in 2:10){
    real.rates %>% mutate(x = i) %>% rbind(real.rates.x) -> real.rates.x
}
real.rates.x %>% arrange(x) -> real.rates.x
real.rates.x %>% group_by(dataset, model, x) %>% spread(type, MLE) -> real.rates.x.spread
    

## "meanbl.rates"                no replicates
## "gammabl.rates"               no replicates
## "permuted.bl.rates"           10 replicates
## "random.tree.rates"           10 replicates
## "randomized.AA.rates"         10 replicates   
## "randomized.sequences.rates"  10 replicates   
## "real.rates"                  truf


randomized.AA.rates %>% 
    spread(type, MLE) %>%
    left_join(real.rates.x.spread) %>% 
    group_by(dataset, model, x) %>%
    summarize(rho = cor(randomizedAA, real, method = "spearman")) %>%
    ggplot(aes(x = dataset, fill = model, y = rho)) + 
        geom_point(pch=21, position = position_jitterdodge()) + 
        ylab("Rank correlation") + xlab("Dataset, arbitrary order") +
        background_grid() -> randomized.aa.correlations
ggsave(paste0(plot.path, "randomized_AA_correlations.pdf"), randomized.aa.correlations)
        
    
randomized.sequences.rates %>% 
    spread(type, MLE) %>%
    left_join(real.rates.x.spread) %>% 
    group_by(dataset, model, x) %>%
    summarize(rho = cor(randomized_sequences, real, method = "spearman")) %>%
    ggplot(aes(x = dataset, fill = model, y = rho)) + 
        geom_point(pch=21, position = position_jitterdodge()) + 
        ylab("Rank correlation") + xlab("Dataset, arbitrary order") +
        background_grid() -> randomized.sequences.correlations
ggsave(paste0(plot.path, "randomized_sequences_over_tree_correlations.pdf"), randomized.sequences.correlations)
 


random.tree.rates %>% 
    spread(type, MLE) %>%
    left_join(real.rates.x.spread) %>% 
    group_by(dataset, model, x) %>%
    summarize(rho = cor(random_trees, real, method = "spearman")) %>%
    ggplot(aes(x = dataset, fill = model, y = rho)) + 
        geom_point(pch=21, position = position_jitterdodge()) + 
        ylab("Rank correlation") + xlab("Dataset, arbitrary order") +
        background_grid() -> random.trees.correlations
ggsave(paste0(plot.path, "random_trees_correlations.pdf"), random.trees.correlations)
 



gammabl.rates %>% 
    rbind(real.rates) %>% 
    group_by(model, dataset) %>%
    spread(type, MLE) %>%
    summarize(rho = cor(gammabl, real, method = "spearman")) %>%
    ggplot(aes(x = dataset, fill = model, y = rho)) + 
        geom_bar(stat = "identity", position = position_dodge()) + 
        coord_cartesian(ylim=c(0.9,1)) +
        ylab("Rank correlation, NOTE SCALE") + xlab("Dataset, arbitrary order") + 
        ggtitle("***single replicate***") -> gamma.bl.correlations
ggsave(paste0(plot.path, "gamma_bl_correlations.pdf"), gamma.bl.correlations)
 



meanbl.rates %>% 
    rbind(real.rates) %>% 
    group_by(model, dataset) %>%
    spread(type, MLE) %>%
    summarize(rho = cor(meanbl, real, method = "spearman")) %>%
    ggplot(aes(x = dataset, fill = model, y = rho)) + 
        geom_bar(stat = "identity", position = position_dodge()) + 
        coord_cartesian(ylim=c(0.9,1)) +
        ylab("Rank correlation, NOTE SCALE") + xlab("Dataset, arbitrary order") + 
        ggtitle("***single replicate (obviously, since mean)***") -> mean.bl.correlations
ggsave(paste0(plot.path, "mean_bl_correlations.pdf"), mean.bl.correlations)
 




























# 
# 
# 
# ######## PERTURBED: COMPLETELY IDENTICAL #########
# for (i in 1:10){
#     print(i)   
#     randomized.rates %>% filter(x==i) %>% select(-x)-> temp.rates
#     rbind(rates, temp.rates) %>%
#         mutate(MLE=ifelse(MLE>=100, 100, MLE)) %>%
#         spread(type, MLE) %>%
#         group_by(model) %>%
#         summarize(rho = cor(randomized, real,method="spearman")) %>% print.data.frame() ## low correlations, when occur, are driven seemingly entirely by HIGH rate outliers, where one gives 100 and another 600. 
# 
#     rbind(rates, temp.rates) %>%
#         filter(MLE > 0) %>%
#         mutate(MLE=ifelse(MLE>=100, 100, MLE)) %>%
#         spread(type, MLE) %>%
#         ggplot(aes(x = real, y = randomized)) + geom_text(aes(label=site), size=1) + scale_y_log10() + scale_x_log10() + geom_abline(color = "blue") +facet_wrap(~model) -> p
#         print(p) ### yeah this is a 1:1 line
#     readline()
# }
# 
# 
# 
# 
# 
# 
# 
# ######## PERTURBED: COMPLETELY IDENTICAL #########
# for (i in 1:10){
#     print(i)   
#     permuted.rates %>% filter(x==i) %>% select(-x)-> temp.rates
#     rbind(rates, temp.rates) %>%
#         mutate(MLE=ifelse(MLE>=100, 100, MLE)) %>%
#         spread(type, MLE) %>%
#         group_by(model) %>%
#         summarize(rho = cor(permuted, real)) %>% print.data.frame() ## low correlations, when occur, are driven seemingly entirely by HIGH rate outliers, where one gives 100 and another 600. 
# 
#     rbind(rates, temp.rates) %>%
#         filter(MLE > 0) %>%
#         mutate(MLE=ifelse(MLE>=50, 50, MLE)) %>%
#         spread(type, MLE) %>%
#         ggplot(aes(x = real, y = permuted)) + geom_point() + geom_abline(color = "blue") +facet_wrap(~model) -> p
#         print(p) ### yeah this is a 1:1 line
#     readline()
# }
# 
# 
# 
# 
# 
# 
# 
# 
# ######## SHUFFLED:  #########
# for (i in 1:10){
#     print(i)   
#     shuffled.rates %>% filter(x==i) %>% select(-x)-> temp.rates
#     rbind(rates, temp.rates) %>%
#         mutate(MLE=ifelse(MLE>=100, 100, MLE)) %>%
#         mutate(MLE=ifelse(MLE<= 1e-8, 1e-8, MLE)) %>%
#         spread(type, MLE) %>%
#         group_by(model) %>%
#         mutate(relerror = abs(real-shuffled)/real) %>% 
#         filter(real > 1e-3, model == "WAG")
#         summarize(sum(relerror <= 1), sum(relerror >1)) %>% print.data.frame()}
#         
#         
#         
#         
#         
#         
#         
#         
#         ggplot(aes(x = relerror)) + geom_histogram() + facet_grid(~model, scales="free") + scale_x_continuous(limits=c(0, 10)) -> p
#     print(p)
#     readline()}
#         
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
#     rbind(rates, temp.rates) %>%
#         spread(type, MLE) %>%
#         group_by(model) %>%
#         summarize(rho = cor(shuffled, real)) %>% print.data.frame() ## low correlations, when occur, are driven seemingly entirely by HIGH rate outliers, where one gives 100 and another 600. 
# }
# 
# 
# 
#     rbind(rates, temp.rates) %>%
#         filter(MLE > 0) %>%
#         mutate(MLE=ifelse(MLE>=100, 100, MLE)) %>%
#         spread(type, MLE) %>%
#         ggplot(aes(x = real, y = shuffled)) + geom_point() + scale_x_log10() + scale_y_log10() + geom_abline(color = "blue") +facet_wrap(~model) -> p
#         print(p) ### yeah this is a 1:1 line
# }
# 
# 
# 
# 







# 
# [1] 1
#         rho
# 1 0.9999993
# [1] 2
#         rho
# 1 0.9999993
# [1] 3
#         rho
# 1 0.9999995
# [1] 4
#     
#         rho
# 1 0.9999993
# [1] 5
#         rho
# 1 0.9999986
# [1] 6
#         rho
# 1 0.9999971
# [1] 7
#         rho
# 1 0.9999993
# [1] 8
#         rho
# 1 0.9999991
# [1] 9
#         rho
# 1 0.9999986
# [1] 10
#         rho
# 1 0.9999992























