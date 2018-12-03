library(readr)
library(dplyr)
library(tidyr)
library(magrittr)
library(here)

# read in expression data
expression <- readRDS(here('data', 'raw', 'metacells.BrainCellAtlas_Saunders_version_2018.04.01.RDS'))
# convert to tibble while creating column with rownames
expression <- as_tibble(expression, rownames = 'gene_symbol')
dim(expression)

# read in metadata
meta <- readRDS(here('data', 'raw', 'annotation.BrainCellAtlas_Saunders_version_2018.04.01.RDS'))
tissue_levels <- list("CB" = "Cerebellum", "FC" = "Frontal Cortex", "PC" = "Posterior Cortex", 
                      "ENT" = "Entopeduncular", "GP" = "Globus Pallidus", "HC" = "Hippocampus", 
                      "STR" = "Striatum", "SN" = "Substantia Nigra", "TH" = "Thalamus")
meta %<>% mutate(tissue_name = recode(meta$tissue, !!!tissue_levels)) 
meta %<>% select(-tissue) 
meta_subset <- meta %>% select(tissue_subcluster, class, tissue_name)
#test <- expression[, names(expression) %in% meta$tissue_subcluster]
#dim(test)

tidy <- gather(expression, key = tissue_subcluster, value = expression, -gene_symbol) 
tidy %<>% mutate(log1Exp=log(1+expression))
tidy %<>% group_by(gene_symbol) %>% 
  mutate(log1ExpZ = (log1Exp - mean(log1Exp)) / sd(log1Exp))

print(dim(tidy))
tidy %<>% filter(!is.na(log1ExpZ))
print(dim(tidy))

#no genes are duplicated -- not needed
#tidy %<>% group_by(gene_symbol, tissue_subcluster) %>% summarise(log1ExpZ2 = mean(log1ExpZ))
#print(dim(tidy))



saunders_ranks_matrix <- tidy %>% 
  group_by(tissue_subcluster) %>% 
  mutate(log1ExpZRank = rank(log1ExpZ)) %>% 
  select(-expression, -log1Exp, -log1ExpZ) %>% 
  spread(tissue_subcluster, value = log1ExpZRank)

dir.create(here('data', 'processed'))
write_csv(saunders_ranks_matrix, here('data', 'processed', 'processed_saunders_ranks.csv'))
saveRDS(saunders_ranks_matrix, file = here('data', 'processed', 'saunders_ranks_matrix.RDS'))
write_csv(meta, here('data', 'processed', 'metadata.csv'))
saveRDS(meta, file = here('data', 'processed', 'meta.RDS'))

# to create separate table for AUC/pvalues by class
tidy <- left_join(tidy, meta_subset, by='tissue_subcluster')
tidy_by_class <- tidy %>% 
  group_by(class, gene_symbol) %>% 
  summarise(exp = mean(expression)) %>% 
  mutate(log1Exp=log(1+exp))
  
tidy_by_class %<>% 
  group_by(gene_symbol) %>% 
  mutate(log1ExpZ = (log1Exp - mean(log1Exp)) / sd(log1Exp))

dim(tidy_by_class)
tidy_by_class %<>% filter(!is.na(log1ExpZ))

class_ranks_matrix <- tidy_by_class %>% 
  group_by(class) %>% 
  mutate(log1ExpZRank = rank(log1ExpZ)) %>% 
  select(-exp, -log1Exp, -log1ExpZ) %>% 
  spread(class, value = log1ExpZRank)

dim(class_ranks_matrix)

write_csv(class_ranks_matrix, here('data', 'processed', 'processed_saunders_classes_ranks.csv'))
saveRDS(class_ranks_matrix, file = here('data', 'processed', 'saunders_classes_ranks_matrix.RDS'))
