##############################################################################
# Script information                                                      
# Title: Tidy PAR1 region results
# Author: Seyhan Yazar
# Date: 2021-09-01
# Description: 
##############################################################################

library(data.table)
library(dplyr)

# Set working directory
setwd("/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/sex_specific_eQTL_analysis/results/2022-12-09_sex-chr-analysis")

# Extract PAR2_XXY results

file_list <- list.files(path="PAR2-XXY", pattern="*_par2_significant_correlation_results.tsv", full.name=TRUE, recursive=TRUE)
celltypes <- sub("PAR2-XXY/*", "", file_list) %>% sub("/.*","", .)

df_list <- lapply(file_list, function(x) { tryCatch(fread(x) , error=function(e) NULL)})
names(df_list) <- celltypes
df_list <- Filter(function(x) dim(x)[1] > 0, df_list)
df <- bind_rows(df_list, .id="celltype")
eSNP1_XXY <- df %>%
    group_by(celltype, geneid) %>%
    arrange(qvalue) %>%
    filter(row_number()==1)

colnames(eSNP1_XXY)<-paste(colnames(eSNP1_XXY),"XXY",sep="_")

df_xx <- list()
for (i in 1:nrow(eSNP1_XXY)) {
    ctype <- eSNP1_XXY$celltype_XXY[i]
    gene <- eSNP1_XXY$geneid_XXY[i]
    snp <- eSNP1_XXY$snpid_XXY[i]
    df_xx[[i]] <- fread(sprintf("PAR2-XX/%s/%s_female_par2_correlation_results.tsv", ctype, ctype)) %>%
        filter(geneid==gene) %>% 
        filter(snpid==snp) %>% 
        mutate(celltype=ctype) %>%
        data.frame()         
    } 
df_xx <- bind_rows(df_xx)
colnames(df_xx)<-paste(colnames(df_xx),"XX",sep="_")

df_xy <- list()
for (i in 1:nrow(eSNP1_XXY)) {
    ctype <- eSNP1_XXY$celltype_XXY[i]
    gene <- eSNP1_XXY$geneid_XXY[i]
    snp <- eSNP1_XXY$snpid_XXY[i]
    if (file.exists(sprintf("PAR2-XY/%s/%s_male_par2_correlation_results.tsv", ctype, ctype))) {
        df_xy[[i]] <- fread(sprintf("PAR2-XY/%s/%s_male_par2_correlation_results.tsv", ctype, ctype)) %>%
        filter(geneid==gene) %>% 
        filter(snpid==snp) %>% 
        mutate(celltype=ctype) %>%
        data.frame() }
    else {
        print('Does not exist')
    }   
} 
df_xy <- bind_rows(df_xy)
colnames(df_xy)<-paste(colnames(df_xy),"XY",sep="_")

eSNP1 <- left_join(eSNP1_XXY, df_xx, by=c("celltype_XXY"="celltype_XX", "geneid_XXY"="geneid_XX", "snpid_XXY"="snpid_XX"))
eSNP1 <- left_join(eSNP1, df_xy, by=c("celltype_XXY"="celltype_XY", "geneid_XXY"="geneid_XY", "snpid_XXY"="snpid_XY"))
eSNP1 <- eSNP1 %>% select(-c("method_XXY","alternative_XXY","method_XX", "alternative_XX", "method_XY","alternative_XY"))

fwrite(eSNP1,sprintf("PAR2-XXY/PAR2-XXY_top_eqtls.tsv"),sep="\t",quote=F)

rm(list=ls())

# Extract PAR1_XX results

file_list <- list.files(path="PAR2-XX", pattern="*_par2_significant_correlation_results.tsv", full.name=TRUE, recursive=TRUE)
celltypes <- sub("PAR1-XX/*", "", file_list) %>% sub("/.*","", .)

df_list <- lapply(file_list, function(x) { tryCatch(fread(x) , error=function(e) NULL)})
names(df_list) <- celltypes
df_list <- Filter(function(x) dim(x)[1] > 0, df_list)
df <- bind_rows(df_list, .id="celltype")
eSNP1_XX <- df %>%
    group_by(celltype, geneid) %>%
    arrange(qvalue) %>%
    filter(row_number()==1)

colnames(eSNP1_XX)<-paste(colnames(eSNP1_XX),"XX",sep="_")

df_xy <- list()
for (i in 1:nrow(eSNP1_XX)) {
    ctype <- eSNP1_XX$celltype_XX[i]
    gene <- eSNP1_XX$geneid_XX[i]
    snp <- eSNP1_XX$snpid_XX[i]
    df_xy[[i]] <- fread(sprintf("PAR2-XY/%s/%s_male_par2_correlation_results.tsv", ctype, ctype)) %>%
        filter(geneid==gene) %>% 
        filter(snpid==snp) %>% 
        mutate(celltype=ctype) %>%
        data.frame()         
    } 
df_xy <- bind_rows(df_xy)
colnames(df_xy)<-paste(colnames(df_xy),"XY",sep="_")

df_xxy <- list()
for (i in 1:nrow(eSNP1_XX)) {
    ctype <- eSNP1_XX$celltype_XX[i]
    gene <- eSNP1_XX$geneid_XX[i]
    snp <- eSNP1_XX$snpid_XX[i]
    df_xxy[[i]] <- fread(sprintf("PAR2-XXY/%s/%s_joint_par2_correlation_results.tsv", ctype, ctype)) %>%
        filter(geneid==gene) %>% 
        filter(snpid==snp) %>% 
        mutate(celltype=ctype) %>%
        data.frame()         
    } 
df_xxy <- bind_rows(df_xxy)
colnames(df_xxy)<-paste(colnames(df_xxy),"XXY",sep="_")

eSNP1 <- left_join(eSNP1_XX, df_xy, by=c("celltype_XX"="celltype_XY", "geneid_XX"="geneid_XY", "snpid_XX"="snpid_XY"))
eSNP1 <- left_join(eSNP1, df_xxy, by=c("celltype_XX"="celltype_XXY", "geneid_XX"="geneid_XXY", "snpid_XX"="snpid_XXY"))
eSNP1 <- eSNP1 %>% select(-c("method_XXY","alternative_XXY","method_XX", "alternative_XX", "method_XY","alternative_XY"))

fwrite(eSNP1,sprintf("PAR2-XX/PAR2-XX_top_eqtls.tsv"),sep="\t",quote=F)

rm(list=ls())

# Extract PAR1_XY results

file_list <- list.files(path="PAR2-XY", pattern="*_par2_significant_correlation_results.tsv", full.name=TRUE, recursive=TRUE)
celltypes <- sub("PAR1-XY/*", "", file_list) %>% sub("/.*","", .)

df_list <- lapply(file_list, function(x) { tryCatch(fread(x) , error=function(e) NULL)})
names(df_list) <- celltypes
df_list <- Filter(function(x) dim(x)[1] > 0, df_list)
df <- bind_rows(df_list, .id="celltype")
eSNP1_XY <- df %>%
    group_by(celltype, geneid) %>%
    arrange(qvalue) %>%
    filter(row_number()==1)

colnames(eSNP1_XY)<-paste(colnames(eSNP1_XY),"XY",sep="_")

df_xx <- list()
for (i in 1:nrow(eSNP1_XY)) {
    ctype <- eSNP1_XY$celltype_XY[i]
    gene <- eSNP1_XY$geneid_XY[i]
    snp <- eSNP1_XY$snpid_XY[i]
    df_xx[[i]] <- fread(sprintf("PAR2-XX/%s/%s_female_par2_correlation_results.tsv", ctype, ctype)) %>%
        filter(geneid==gene) %>% 
        filter(snpid==snp) %>% 
        mutate(celltype=ctype) %>%
        data.frame()         
    } 
df_xx <- bind_rows(df_xx)
colnames(df_xx)<-paste(colnames(df_xx),"XX",sep="_")

df_xxy <- list()
for (i in 1:nrow(eSNP1_XY)) {
    ctype <- eSNP1_XY$celltype_XY[i]
    gene <- eSNP1_XY$geneid_XY[i]
    snp <- eSNP1_XY$snpid_XY[i]
    df_xxy[[i]] <- fread(sprintf("PAR2-XXY/%s/%s_joint_par2_correlation_results.tsv", ctype, ctype)) %>%
        filter(geneid==gene) %>% 
        filter(snpid==snp) %>% 
        mutate(celltype=ctype) %>%
        data.frame()         
    } 
df_xxy <- bind_rows(df_xxy)
colnames(df_xxy)<-paste(colnames(df_xxy),"XXY",sep="_")

eSNP1 <- left_join(eSNP1_XY, df_xx, by=c("celltype_XY"="celltype_XX", "geneid_XY"="geneid_XX", "snpid_XY"="snpid_XX"))
eSNP1 <- left_join(eSNP1, df_xxy, by=c("celltype_XY"="celltype_XXY", "geneid_XY"="geneid_XXY", "snpid_XY"="snpid_XXY"))
eSNP1 <- eSNP1 %>% select(-c("method_XXY","alternative_XXY","method_XX", "alternative_XX", "method_XY","alternative_XY"))

fwrite(eSNP1,sprintf("PAR2-XY/PAR2-XY_top_eqtls.tsv"),sep="\t",quote=F)