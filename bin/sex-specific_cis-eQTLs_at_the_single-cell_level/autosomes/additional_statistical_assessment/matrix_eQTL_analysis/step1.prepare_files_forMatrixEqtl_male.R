##############################################################################
# Script information                                                      
# Title: Matrix eQTL analysis
# Author: Seyhan Yazar
# Date: 2021-09-01
# Description: 
##############################################################################
# Call the cell type
args = commandArgs(trailingOnly=TRUE)
cellType <- args[1]
# cellType <- "B_intermediate"
print(cellType)
chrNumber <- args[2]
# chrNumber <- "22"
print(chrNumber)

# Set working directory
setwd("/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/sex_specific_eQTL_analysis")

# Import libraries ------------------------------------
library(tidyverse)
library(dsLib)

# Set output ------------------------------------------
output <- set_output("2022-10-04", "matrix-eqtl-analysis-male")

out.dir <- paste0(output, '/', cellType)
dir.create(out.dir)

library("dplyr")
library("data.table")
library("BEDMatrix")

# rs ids
# df_hrc <- readRDS ("/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/cell_specific_eQTL_analysis_April20/hrc_ids_all.rds")

# Set the directory path
data.dir <- '/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/cell_specific_eQTL_analysis_June21/'

# Prepare expression file 
expression <- fread(sprintf("%s/outputs/filtered_matrices/%s_expression_genes_removed.tsv", data.dir, cellType))
# Subset expression data for the individuals tested in original analysis
residuals <- fread(sprintf("/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/sex_specific_eQTL_analysis/results/2022-09-27_male-only/%s/round1/%s_chr22_log_residuals_male.tsv", cellType,cellType))
expression <- expression %>% filter(sampleid %in% residuals$sampleid)
dim(expression)
sampleids <- expression$sampleid
t_expression <- data.frame(t(expression[,-c("sampleid")]), check.names=F)
colnames(t_expression) <- sampleids

# Prep gene location file
my_genes <- fread('/directflow/SCCGGroupShare/projects/SeyhanYazar/covid19/matrix_eqtl_27_cells/results/2021-09-01_gene-locations/gene_locations_all.tsv')
my_genes <- my_genes %>% filter(gene_name %in% row.names(t_expression))
my_genes <- my_genes %>% select(-'strand') %>% filter(chr==paste0("chr",chrNumber))
fwrite(my_genes, paste0(out.dir, '/gene_chr',chrNumber, ".tsv"), quote=F,row.names=F, sep="\t")

# subset expression dataset
t_expression  <- t_expression[c(row.names(t_expression) %in% my_genes$gene_name),]
t_expression$id <- rownames(t_expression)
t_expression <- t_expression %>% select(id, everything())
fwrite(t_expression, paste0(out.dir, '/expression_chr',chrNumber, ".tsv"), quote=F,row.names=F, sep="\t")

# Prep covariate file
covariate_df <- fread(sprintf("%s/outputs/peer_factors/%s_peer_factors_after_removing_genes.tsv", data.dir, cellType)) %>% 
    filter(sampleid %in% residuals$sampleid) %>%
    select(-c(sex,pf3:pf10))
colnames(covariate_df)[1] <- "id"
t_covariate <- data.frame(t(covariate_df))
t_covariate$vars <- rownames(t_covariate)
t_covariate <- t_covariate %>% select(vars, everything())

fwrite(t_covariate, paste0(out.dir, '/covariate_chr',chrNumber, ".tsv"), quote=F,row.names=F, col.names=F, sep="\t")

# Prepare Genotype file
chr <- BEDMatrix(paste0('/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/imputed_data/decompressed/filter_vcf_r08_maf005/plink_chr',chrNumber,'.bed',collapse=''))
chrMatrix <- t(as.matrix (chr))
n <- gsub("^0_*","", colnames(chrMatrix))
colnames (chrMatrix) <- n
chrMatrix <- chrMatrix [,c(sampleids)]
chrMatrix <- data.frame(chrMatrix, check.names=F)
chrMatrix$snpid <- rownames(chrMatrix)
genotype <- chrMatrix %>% select (snpid, everything())
colnames(genotype)[1] <- "id"
fwrite(genotype, paste0(out.dir, '/SNPs_chr',chrNumber, ".tsv"), quote=F,row.names=F, sep="\t")

# Prepare snp location file
# snp_df <- fread(sprintf("/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/cell_specific_eQTL_analysis_April20/%s/linear_regression/input_files/snpsloc_chr%s.tsv", cellType, chrNumber))
snps<- data.frame(chrMatrix$snpid)
snps <- cbind(snps, read.table(text = as.character(snps$chrMatrix.snpid), sep = ":"))
snps <- cbind(snps, read.table(text = as.character(snps$V2), sep = "_"))
snps <- snps [,c(1,2,4)]
snps$V1 <- paste0('chr',chrNumber)
colnames(snps) <- c("snpid","chr","pos")
colnames(snps)[1] <- "id"
fwrite(snps, paste0(out.dir, '/snpsloc_chr',chrNumber, ".tsv"), quote=F,row.names=F, sep="\t")
