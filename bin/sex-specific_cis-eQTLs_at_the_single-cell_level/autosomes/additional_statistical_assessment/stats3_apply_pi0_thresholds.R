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
    "dnT","gdT", "MAIT", "Treg",
    "NK", "NK_CD56bright", "NK_Proliferating", 
    "Plasmablast"
    )

### Female vs Male ###
pi_results <- fread("2022-10-18_pi-analysis/female_vs_male_pi0_values.tsv")
female_eqtls <- fread("2022-10-12_female-only/female_eqtls/female_top_eqtls.tsv")
male_eqtls <- fread("2022-09-27_male-only/male_eqtls/male_top_eqtls.tsv")
joint_eqtls <- fread("/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/sex_specific_eQTL_analysis/results/2022-09-27_joint_analysis_one_round/joint_analysis_one_round_results.tsv")
joint_eqtls <- joint_eqtls %>% 
    mutate(snpid=paste0(Chromosome,":",Position,"_", SNP_assessed_allele))
    colnames(joint_eqtls) [1:2] <- c("celltype", "geneid")

rm_false_positives <- function (x) { 
    # x=celltypes[1]
    female_results_df <- female_eqtls %>% filter(celltype==x)
    dim(female_results_df)
    male_results_df <- male_eqtls %>% filter(celltype==x)
    joint_results_df <- joint_eqtls %>% filter(celltype==x)
    shared_eqtls <- female_results_df[female_results_df$snpid %in% male_results_df$snpid,]
    dim(shared_eqtls)
    shared_eqtls <- shared_eqtls [!shared_eqtls$snpid %in% joint_results_df$snpid,]
    dim(shared_eqtls)
    female_only_eqtls <-female_results_df[!female_results_df$snpid %in% joint_results_df$snpid,]
    dim(female_only_eqtls)
    female_only_eqtls <- female_only_eqtls[!female_only_eqtls$snpid %in% male_results_df$snpid,]
    dim(female_only_eqtls)
    not_male_eqtl_df <- fread(sprintf("/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/sex_specific_eQTL_analysis/results/2022-10-18_pi-analysis/p_outputs/female_vs_male_%s_pvalues.tsv", x))
    dim(not_male_eqtl_df)

    if (x %in% pi_results$celltypes) {
        pi_cut_off = pi_results$pi0[pi_results$celltype==x]
        df <- not_male_eqtl_df %>% 
            arrange(qvalue) %>% 
            slice_head(prop=(1-pi_cut_off))
    } else {
        pi_cut_off = pi_results$pi0[pi_results$celltype=="All"]
        df <- not_male_eqtl_df %>% 
            arrange(qvalue) %>% 
            slice_head(prop=(1-pi_cut_off))
    }
    dim(df)

    not_shared_eqtls <- female_only_eqtls[!female_only_eqtls$snpid %in% df$snpid,]
    print(x)
    print(dim(not_shared_eqtls))
     shared_eqtls$uniqueness <- "common"
    not_shared_eqtls$uniqueness <- "female_only"
    rbind(shared_eqtls, not_shared_eqtls)
}

culled_results <- lapply(celltypes, rm_false_positives)
culled_results_df <- bind_rows(culled_results)
fwrite(culled_results_df, "2022-10-18_pi-analysis/female_eQTLs_after_pi0_threshold.tsv", sep="\t")

### Male vs Female ###
rm(list=ls())

setwd("/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/sex_specific_eQTL_analysis/results")

celltypes <- c("B_intermediate", "B_memory", "B_naive", 
    "CD14_Mono", "CD16_Mono", 
    "CD4_CTL", "CD4_Naive", "CD4_TCM","CD4_TEM",
    "CD8_Naive", "CD8_TCM", "CD8_TEM", "DC",
    "dnT","gdT", "MAIT", "Treg",
    "NK", "NK_CD56bright", "NK_Proliferating", 
    "Plasmablast"
    )

pi_results <- fread("2022-10-18_pi-analysis/male_vs_female_pi0_values.tsv")
male_eqtls <- fread("2022-09-27_male-only/male_eqtls/male_top_eqtls.tsv")
female_eqtls <- fread("2022-10-12_female-only/female_eqtls/female_top_eqtls.tsv")
joint_eqtls <- fread("/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/sex_specific_eQTL_analysis/results/2022-09-27_joint_analysis_one_round/joint_analysis_one_round_results.tsv")
joint_eqtls <- joint_eqtls %>% 
    mutate(snpid=paste0(Chromosome,":",Position,"_", SNP_assessed_allele))
    colnames(joint_eqtls) [1:2] <- c("celltype", "geneid")

rm_false_positives <- function (x) { 
    # x=celltypes[1]
    female_results_df <- female_eqtls %>% filter(celltype==x)
    dim(female_results_df)
    male_results_df <- male_eqtls %>% filter(celltype==x)
    joint_results_df <- joint_eqtls %>% filter(celltype==x)
    shared_eqtls <- male_results_df[male_results_df$snpid %in% female_results_df$snpid,]
    dim(shared_eqtls)
    shared_eqtls <- shared_eqtls [!shared_eqtls$snpid %in% joint_results_df$snpid,]
    dim(shared_eqtls)
    male_only_eqtls <- male_results_df[!male_results_df$snpid %in% joint_results_df$snpid,]
    dim(male_only_eqtls)
    male_only_eqtls <-  male_only_eqtls[! male_only_eqtls$snpid %in% female_results_df$snpid,]
    dim(male_only_eqtls)
    not_female_eqtl_df <- fread(sprintf("/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/sex_specific_eQTL_analysis/results/2022-10-18_pi-analysis/p_outputs/male_vs_female_%s_pvalues.tsv", x))
    dim(not_female_eqtl_df)

    if (x %in% pi_results$celltypes) {
        pi_cut_off = pi_results$pi0[pi_results$celltype==x]
        df <- not_female_eqtl_df %>% 
            arrange(qvalue) %>% 
            slice_head(prop=abs(1-pi_cut_off))
    } else {
        pi_cut_off = pi_results$pi0[pi_results$celltype=="All"]
        df <- not_female_eqtl_df %>% 
            arrange(qvalue) %>% 
            slice_head(prop=abs(1-pi_cut_off))
    }
    dim(df)

    not_shared_eqtls <- male_only_eqtls[!male_only_eqtls$snpid %in% df$snpid,]
    print(x)
    print(dim(not_shared_eqtls))
    shared_eqtls$uniqueness <- "common"
    not_shared_eqtls$uniqueness <- "male_only"
    rbind(shared_eqtls, not_shared_eqtls)
}

culled_results <- lapply(celltypes, rm_false_positives)
culled_results_df <- bind_rows(culled_results)
fwrite(culled_results_df, "2022-10-18_pi-analysis/male_eQTLs_after_pi0_threshold.tsv", sep="\t")
