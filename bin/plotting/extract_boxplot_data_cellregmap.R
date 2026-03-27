#!/usr/bin/env Rscript

# Set working directory
setwd("/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/sex_specific_eQTL_analysis")

# Import libraries ------------------------------------
library(tidyverse)
library(dsLib)

# Set output ------------------------------------------
output <- set_output("2022-07-13", "sex-specific-eQTL-expr-geno-plots")

library("dplyr")
library("data.table")
library("tidyverse")

# Call the arguments
args = commandArgs(trailingOnly=TRUE)

genename <- args[1]
snp <- args[2]

# snp <- "12:56435929_G"
# genename <- "RPS26"

# Cell Labels
cellLabels <- c("B_intermediate", "B_memory", "B_naive", "CD14_Mono", "CD16_Mono", "CD4_CTL", "CD4_Naive", "CD4_TCM", "CD4_TEM",
    "CD8_Naive", "CD8_TCM", "CD8_TEM" , "dnT", "gdT", "MAIT", "NK", "NK_CD56bright", "NK_Proliferating", "DC", "Plasmablast", "Treg")

chr <- gsub(":.*", "", snp)

# Expression filename
get_expressions <- function (x) {
    # x="CD4_Naive"
    cDir <- sprintf("/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/sex_specific_eQTL_analysis/results/2022-03-14_conditional-analysis-stratified/%s",x)
    #expression_filename <- sprintf("%s/round1/%s_expression.tsv", cDir, cellLabel)
    # For females:
    expression_filename_f <- sprintf("%s/round1/%s_chr%s_log_residuals_female.tsv", cDir,x, chr)
    expression_df_f <- fread(expression_filename_f)
    if(genename %in% colnames(expression_df_f)){
    expression_df2_f <- expression_df_f %>% select (geneid=all_of(genename), sampleid)
    expression_df2_f$sex <- "female" }
    else {       
    }
    # For males:
    expression_filename_m <- sprintf("%s/round1/%s_chr%s_log_residuals_male.tsv", cDir,x, chr)
    expression_df_m <- fread(expression_filename_m)
    if(genename %in% colnames(expression_df_m)){
    expression_df2_m <- expression_df_m %>% select (geneid=all_of(genename), sampleid)
    expression_df2_m$sex <- "male"}
    else {
    }

    expression_df3 <- bind_rows(expression_df2_f, expression_df2_m)
    expression_df3
}

combine_expressions <- lapply (cellLabels, get_expressions)
names(combine_expressions) <- cellLabels
l <- lapply(combine_expressions, function(x) if(is.null(x)) data.frame(geneid = NA, sampleid= NA) else x)
combine_expressions_df <- bind_rows(l, .id="cellType")
dim(combine_expressions_df)

# Set the directory path
genotype_dir <- '/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/cell_specific_eQTL_analysis_April20/Genotype_Files'

# Prepare genotype file
genotype <- fread(sprintf("%s/genotype_chr%s.tsv", genotype_dir, chr))
genotype_df <- genotype %>% select("sampleid", all_of(snp))
colnames(genotype_df)[2] <- "genotype"

expr_geno_df <- left_join(combine_expressions_df, genotype_df, by="sampleid")
colnames(expr_geno_df)[2] <- "expression"

df_snps <- fread(sprintf("/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/imputed_data/decompressed/snps_with_maf_greaterthan0.05/chr%s.SNPs.txt", chr))
eSNP_letter <- gsub("_.*", "", snp)
df_snps <- df_snps %>% filter(SNP==eSNP_letter)
A1 <- df_snps$'REF(0)'
A2 <- df_snps$'ALT(1)'
expr_geno_df$genotype <- factor(expr_geno_df$genotype, labels=c(paste0(A1,A1), paste0(A1,A2), paste0(A2,A2)))

# Get the snpid
df_hrc <- readRDS ("/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/cell_specific_eQTL_analysis_April20/hrc_ids_all.rds")
eSNP_df <- as.data.frame(df_hrc) %>% filter(snpid==snp)
rsid <- eSNP_df$ID

fwrite(expr_geno_df, sprintf("%s/%s_%s_plot_data.tsv", output,genename, rsid), sep="\t", quote=F) 

print("JOB IS DONE!")

quit()
