##############################################################################
# Script information                                                      
# Title: Conditional cis-eQTL mapping - round 1
# Author: Seyhan Yazar
# Date: 2020-12-23
# Description: This R script was written to run sex-biased eQTLs by adding an interaction 
# term. Only first round of conditional analysis outputs are used.
##############################################################################

# Set working directory
setwd("/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/sex_specific_eQTL_analysis")

# Import libraries ------------------------------------
library(tidyverse)
library(dsLib)

# Set output ------------------------------------------
output <- set_output("2022-09-27", "sex-biased-eQTL-analysis")

library(data.table)
library(dplyr)
library(tidymodels)

args = commandArgs(trailingOnly=TRUE)

chunk_start <- args[1]
chunk_end <- args[2]

cis_eqtls <- fread("results/2022-09-27_joint_analysis_one_round/cis_assoc_to_test_for_sex_bias.tsv")

cis_eqtls <- cis_eqtls[chunk_start:chunk_end,]

# cis_eqtls <- cis_eqtls[14501:14510,]


identify_sb_cis_eqtl <- function(rownumber){
    # rownumber="1"
    # print(rownumber)
    cis_df <- cis_eqtls %>% filter(row (cis_eqtls) == rownumber)
    celltype <- cis_df$cell_type
    gene <- cis_df$geneid
    snp <- cis_df$snpid
    chr <- cis_df$chr
    
    expression_filename <- sprintf("/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/cell_specific_eQTL_analysis_June21/outputs/filtered_matrices/%s_expression_genes_removed.tsv", celltype)

    expression_df <- fread(expression_filename)
    expression_df <- expression_df %>% select ('sampleid', all_of(gene))

    covariate_df <- fread(sprintf("/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/cell_specific_eQTL_analysis_June21/outputs/peer_factors/%s_peer_factors_after_removing_genes.tsv", celltype))
    
    genotype_df <- fread(sprintf("/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/cell_specific_eQTL_analysis_April20/Genotype_Files/genotype_chr%s.tsv", chr))
    genotype_df <- genotype_df %>% select('sampleid', all_of(snp))
    
    lm_df <- purrr::reduce(list(expression_df,genotype_df, covariate_df), left_join, by = "sampleid")
    colnames(lm_df)[2] <- 'expression'
    colnames(lm_df)[3] <- 'genotype'
    
    simple_model <- lm(log(expression + 1) ~ genotype + sex + age + pc1 + pc2 + pc3 + pc4 + pf1 + pf2 , data=lm_df)
    simple_model_results <- tidy(simple_model) %>% filter(term %in% c("genotype", "sex"))
    simple_model_results$model <- 'simple' 
    simple_model_results <- simple_model_results %>% pivot_wider(names_from=term, values_from=c(estimate,std.error,statistic,p.value))
    simple_model_results <- cbind(cis_df, simple_model_results)

    interaction_model <- lm(expression ~ genotype + sex + genotype*sex + age + pc1 + pc2 + pc3 + pc4 + pf1 + pf2, data=lm_df)
    interaction_model_results <- tidy(interaction_model) %>% filter(term %in% c("genotype", "sex", "genotype:sex"))
    interaction_model_results$model <- 'interaction' 
    interaction_model_results <- interaction_model_results %>% pivot_wider(names_from=term, values_from=c(estimate,std.error,statistic,p.value))
    interaction_model_results <- cbind(cis_df, interaction_model_results)

    output <- list(simple_model_results,interaction_model_results)
    return(output)
}

both_models_cis_eqtl <- lapply(1:nrow(cis_eqtls), identify_sb_cis_eqtl)

both_models_cis_eqtl_df <- bind_rows(both_models_cis_eqtl)

fwrite(both_models_cis_eqtl_df, paste0("/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/sex_specific_eQTL_analysis/results/2022-09-27_sex-biased-eQTL-analysis/simple_and_interaction_models_",chunk_start,"-",chunk_end,".tsv"), sep="\t")