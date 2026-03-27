#!/usr/bin/env Rscript

# Upload libraries
library(data.table)
library(dplyr)

setwd("/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/sex_specific_eQTL_analysis/results")

# Cell types
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

extract_indv <- function(x) {
    # x=celltypes[1]
    dir.create(sprintf("frequencies/%s",x))
    males <- fread(sprintf("2022-03-14_conditional-analysis-stratified/%s/round1/%s_chr22_log_residuals_male.tsv", x, x))
    males_sampleids <- males %>% select(sampleid) %>% mutate(sex="M")
    females <- fread(sprintf("2022-03-14_conditional-analysis-stratified/%s/round1/%s_chr22_log_residuals_female.tsv", x, x))
    females_sampleids <- females %>% select(sampleid) %>% mutate(sex="F")
    samples <- rbind(males_sampleids, females_sampleids)
    samples$fid <- "0"
    samples <- samples %>% select(fid, sampleid, sex)
    fwrite(samples, sprintf("frequencies/%s/%s_samples.tsv",x,x), col.names=F, sep="\t")
}

lapply(celltypes, extract_indv)

setwd("/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/sex_specific_eQTL_analysis/results/pi_analysis/p_outputs")
ExpFileList <- list.files(path = ".", pattern = "_pvalues.tsv")

extract_snps <- function(y) {
    #y=ExpFileList[22]
    fname <- sub("*.pvalues.tsv", "", y)
    celltype <- sub("female_vs_male*.","", fname)
    pi0_results <- fread(y)
    pi0_results_df <- pi0_results %>% 
        select("snpid") %>% 
        mutate(snp=gsub("_.*", "", snpid)) %>% 
        mutate(chr=gsub(":.*","", snp))
    chr_list <- unique(pi0_results_df$chr)
    for (i in chr_list) {
        sub_df <- pi0_results_df %>% 
            filter(chr==i) %>% 
            select("snp")
        fwrite(sub_df,sprintf("/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/sex_specific_eQTL_analysis/results/frequencies/%s/%s_chr%s_snps.tsv", celltype, celltype, i), col.names=F, sep="\t")
    }
}

lapply(ExpFileList[22:42],extract_snps)