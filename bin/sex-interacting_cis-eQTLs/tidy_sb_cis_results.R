# Upload libraries
library(data.table)
library(dplyr)
library(qvalue)

setwd("/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/sex_specific_eQTL_analysis/results/2022-09-27_sex-biased-eQTL-analysis")

files <- list.files(path=".", pattern="simple_and_interaction_models_*")
dataset <- lapply(files, function(x) { tryCatch(fread(x) , error=function(e) NULL)})

# Get file names and add to the list
filenames <- sub("simple_and_interaction_models_*", "", files) %>% sub('*.tsv', '', .)
names(dataset) <- filenames
dataset2 <- bind_rows(dataset)

df_interaction <- dataset2 %>% filter(model=='interaction') 
sb_cis_eQTLs <- df_interaction[df_interaction$"p.value_genotype:sex" < 0.05,]
sb_cis_eQTLs %>% 
    group_by(cell_type) %>% 
    count() %>% 
    data.frame()

library(tidyverse)

df_interaction$pvalues <- df_interaction$'p.value_genotype:sex'

df_bf <- df_interaction %>% 
    group_by(cell_type, geneid) %>% 
    mutate(adjusted_p = p.adjust(pvalues, method="bonferroni"))

df_bf_top <- df_bf %>% 
    group_by(cell_type,geneid) %>%
    slice(which.min(adjusted_p))
    
df_c_fdr <- df_bf_top %>% 
    filter(!cell_type %in% c("HSPC", "Platelet")) %>%
    group_by(cell_type) %>% 
    mutate(fdr = qvalue(adjusted_p)$lfdr)

df_25 <- df_c_fdr[df_c_fdr$fdr < 0.25,] %>% data.frame()

fwrite(df_25,"sb_cis_eQTLs_at_25perc_fdr_20220928.tsv", sep="\t")