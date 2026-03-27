##############################################################################
# Script information                                                      
# Title: sex specific eqtl analysis
# Author: Seyhan Yazar
# Date: 2020-10-18
# Description: 
##############################################################################

# Set working directory
setwd("/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/sex_specific_eQTL_analysis")

# Import libraries ------------------------------------
library(tidyverse)
library(dsLib)

# Set output ------------------------------------------
output <- set_output("2022-12-08", "female-male-common-egenes")

library("data.table")
library("dplyr")

common_egenes <- fread("/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/sex_specific_eQTL_analysis/data/female_male_egene_common.csv", header=T)

snps_to_test <- common_egenes %>% select(f_snpid, m_snpid) %>% unique()
snps <- pivot_longer(snps_to_test, cols=c("f_snpid", "m_snpid")) %>% select(value) 
snps$value <- gsub("*_.", "", snps$value)

fwrite(snps, file=paste0(output,"/snp_list.txt"), col.names=F)

## NK cell results for manhattan plot

