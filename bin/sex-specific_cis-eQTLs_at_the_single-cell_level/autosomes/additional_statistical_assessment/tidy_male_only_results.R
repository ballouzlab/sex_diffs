#!/usr/bin/env Rscript

setwd("/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/sex_specific_eQTL_analysis/results/2022-09-27_male-only")

library("dplyr")
library("data.table")
library("tidyverse")

celltypes <- c("B_intermediate", "B_memory", "B_naive", 
    "CD14_Mono", "CD16_Mono", 
    "CD4_CTL", "CD4_Naive", "CD4_TCM","CD4_TEM",
    "CD8_Naive", "CD8_TCM", "CD8_TEM", "DC",
    #"cDC", "cDC1", "cDC2", "DC", "pDC", "pDC1",
    "dnT","gdT", "MAIT", "Treg",
    "NK", "NK_CD56bright", "NK_Proliferating", 
    "Plasmablast"
    #"Platelet"
    )
dir.create("male_eqtls")

one_file_per_celltype_male <- function (i) {
    m = i
    files <- list.files(path=sprintf("./%s/round1", m), pattern="*_significant_correlation_results_male.tsv", full.names=T)
    print(length(files))
    print(files)
    dataset <- lapply(files, function(x) { tryCatch(fread(x) , error=function(e) NULL)})
    dataset <- Filter(function(x) dim(x)[1] > 0, dataset)
    dataset2 <- bind_rows(dataset, .id="listnumber")
    fwrite(dataset2, sprintf("./male_eqtls/%s_male_eqtls.tsv", m), sep="\t")
    # Identify the top eSNP for each eGene 
    eSNP1 <- dataset2 %>%
    group_by(geneid) %>%
    arrange(qvalue) %>%
    filter(row_number()==1)
    fwrite(eSNP1,sprintf("./male_eqtls/%s_male_top_eqtls.tsv", m),sep="\t",quote=F)
}

sapply(celltypes, one_file_per_celltype_male)

files <- list.files(path="./male_eqtls", pattern="*_male_top_eqtls.tsv",full.names=T)
print(length(files))
print(files)
celltypes <- gsub("*_male_top_eqtls.tsv","", files) %>% gsub("./male_eqtls/*", "", .)
dataset <- lapply(files, function(x) { tryCatch(fread(x) , error=function(e) NULL)})
names(dataset) <- celltypes
dataset <- Filter(function(x) dim(x)[1] > 0, dataset)
dataset2 <- bind_rows(dataset, .id="celltype") %>% select(-c("listnumber"))
fwrite(dataset2,"./male_eqtls/male_top_eqtls.tsv",sep="\t",quote=F)

