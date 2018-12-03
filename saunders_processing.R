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

saunders_ranks_matrix <- tidy %>% group_by(tissue_subcluster) %>% 
  mutate(log1ExpZRank = rank(log1ExpZ)) %>% select(-expression, -log1Exp, -log1ExpZ) %>% 
  spread(tissue_subcluster, value = log1ExpZRank)

dir.create(here('data', 'processed'))
write_csv(saunders_ranks_matrix, here('data', 'processed', 'processed_saunders_ranks.csv'))
saveRDS(saunders_ranks_matrix, file = here('data', 'processed', 'saunders_ranks_matrix.RDS'))
