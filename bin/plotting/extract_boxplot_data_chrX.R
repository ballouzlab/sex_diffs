#!/usr/bin/env Rscript

# Set working directory
setwd("/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/sex_specific_eQTL_analysis")

# Import libraries ------------------------------------
library(tidyverse)
library(dsLib)

# Set output ------------------------------------------
# output <- set_output("2022-11-14", "sex-interacting-eQTL-expr-geno-plots")
# output <- set_output("2022-11-28", "overlap-sex-egenes-plots-females")
output <- set_output("2022-12-14", "expr-geno-plots")

library("dplyr")
library("data.table")
library("tidyverse")

# Call the arguments
args = commandArgs(trailingOnly=TRUE)

celltype <- args[1]
genename <- args[2]
snp <- args[3]

# snp <- "X:119505884:C:G"
# genename <- "SEPT6"
# celltype <- "CD4_Naive"
get_expressions <- function (celltype) {
    #x="NK"
    cDir <-"/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/sex_specific_eQTL_analysis/results/2022-12-09_sex-chr-analysis"
    #expression_filename <- sprintf("%s/round1/%s_expression.tsv", cDir, cellLabel)
    # For females:
    expression_filename_f <- sprintf("%s/NonPAR-XaXi/%s/%s_nonPAR_XaXi_log_residuals.tsv", cDir, celltype,celltype)
    expression_df_f <- fread(expression_filename_f)
    if(genename %in% colnames(expression_df_f)){
    expression_df2_f <- expression_df_f %>% select (geneid=all_of(genename), sampleid)
    expression_df2_f$sex <- "female"} else {
    expression_df2_f<-setNames(data.frame(matrix(ncol = 3, nrow = 0)), c(genename, "sampleid", "sex")) 
    }

     # For males:
    expression_filename_m <- sprintf("%s/NonPAR-X-ig/%s/%s_nonPAR_X_ig_log_residuals.tsv",cDir, celltype,celltype)
    expression_df_m <- fread(expression_filename_m)
    if(genename %in% colnames(expression_df_m)){
    expression_df2_m <- expression_df_m %>% select (geneid=all_of(genename), sampleid)
    expression_df2_m$sex <- "male"} else {
    expression_df2_m<-setNames(data.frame(matrix(ncol = 3, nrow = 0)), c(genename, "sampleid", "sex")) 
    }

    # For joint:
    expression_filename_j <- sprintf("%s/NonPAR-XaXiX/%s/%s_nonPAR_XaXiX_log_residuals.tsv",cDir, celltype,celltype)
    expression_df_j <- fread(expression_filename_j)
    if(genename %in% colnames(expression_df_j)){
    expression_df2_j <- expression_df_j %>% select (geneid=all_of(genename), sampleid)
    expression_df2_j$sex <- "joint"
    }else {
    expression_df2_j<-setNames(data.frame(matrix(ncol = 3, nrow = 0)), c(genename, "sampleid", "sex")) 
    }
    expression_df3 <- bind_rows(expression_df2_f, expression_df2_m)
    expression_df4 <- bind_rows(expression_df3, expression_df2_j)
    expression_df4
} 

expression_res <- get_expressions(celltype)

# Set the directory path
genotype_dir <- '/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/cell_specific_eQTL_analysis_April20/Genotype_Files'

# Prepare genotype file
genotype <- fread(sprintf("./results/2022-12-09_sex-chr-analysis/NonPAR-XaXiX/%s/%s_nonPAR_XaXiX_genotypes.tsv", celltype,celltype))
genotype_df <- genotype %>% select("sampleid", all_of(snp))
colnames(genotype_df)[2] <- "genotype"

expr_geno_df <- left_join(expression_res , genotype_df, by="sampleid")
colnames(expr_geno_df)[2] <- "expression"

A1 <- unlist(stringr::str_split(snp, ":"))[3]
A2 <- unlist(stringr::str_split(snp, ":"))[4]
expr_geno_df$genotype <- factor(expr_geno_df$genotype, labels=c(paste0(A1,A1), paste0(A1,A2), paste0(A2,A2)))

fwrite(expr_geno_df, sprintf("%s/plot_data_%s_%s_%s.tsv", output,celltype,genename,snp), sep="\t", quote=F) 

print("JOB IS DONE!")

quit()
