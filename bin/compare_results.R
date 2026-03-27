#!/usr/bin/env Rscript

setwd("/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/sex_specific_eQTL_analysis/results/2022-03-14_conditional-analysis-stratified")

library("dplyr")
library("data.table")
library("tidyverse")
library("UpSetR")
library("grid")

# Identify lead SNPs in females
files <- list.files(path="./female_eqtls", pattern="*.tsv", full.names=T)
dataset <- lapply(files, function(x) { tryCatch(fread(x) , error=function(e) NULL)})
dataset <- Filter(function(x) dim(x)[1] > 0, dataset)
celltype_order <- gsub("./female_eqtls/*", "", files) 
celltype_order <- gsub("*_female_eqtls.tsv", "", celltype_order)
names(dataset) <- celltype_order
dataset2 <- bind_rows(dataset, .id="celltype")
fleads <- dataset2 %>% 
    group_by(celltype, geneid) %>%
    arrange(qvalue) %>%
    filter(row_number()==1) %>%
    ungroup() %>% data.frame()
fleads$chr <- as.integer(sub(":.*","", fleads$snpid))
fleads <- fleads %>% drop_na()
## fwrite(fleads, "female_top_eQTLs_26May22.tsv", sep="\t")

# Identify lead SNPs in males
files <- list.files(path="./male_eqtls", pattern="*.tsv", full.names=T)
dataset <- lapply(files, function(x) { tryCatch(fread(x) , error=function(e) NULL)})
dataset <- Filter(function(x) dim(x)[1] > 0, dataset)
celltype_order <- gsub("./male_eqtls/*", "", files) 
celltype_order <- gsub("*_male_eqtls.tsv", "", celltype_order)
names(dataset) <- celltype_order
dataset2 <- bind_rows(dataset, .id="celltype")
mleads <- dataset2 %>% 
    group_by(celltype, geneid) %>%
    arrange(qvalue) %>%
    filter(row_number()==1) %>%
    ungroup() %>% data.frame()
mleads$chr <- as.integer(sub(":.*","", mleads$snpid))
mleads <- mleads %>% drop_na()
## fwrite(mleads, "male_top_eQTLs_26May22.tsv", sep="\t")

# Read joint analysis results
# j_chrx <- fread("/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/sex_specific_eQTL_analysis/results/chrX_cis_eQTLs_20220314_v2.tsv")
j_auto <- fread("/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/cell_specific_eQTL_analysis_June21/results/2022-03-03_tidied_results/OneK1K_eQTLs_Results_20220310.tsv")
# jleads <- rbind(j_chrx, j_auto)
jleads <- j_auto %>% filter(eSNP_rank=="eSNP1")

# Combine all three groups in one

mleads <- mleads %>% select ("celltype","geneid", "chr", "snpid", "estimate","p.value","localFDR") %>% ungroup() 
mleads$group <- "male"
fleads <- fleads %>% select ("celltype","geneid", "chr", "snpid", "estimate","p.value","localFDR") %>% ungroup() 
fleads$group <- "female"
jleads$snpid <- paste0(jleads$Chromosome,":",jleads$Position,"_",jleads$SNP_assessed_allele)
jleads <- jleads %>% select ("cell_type", "GeneID", "Chromosome", "snpid", "rho_correlation_coefficient", "pvalue", "FDR")
colnames(jleads) <- c("celltype","geneid", "chr", "snpid", "estimate","p.value","localFDR")
jleads$group <- "joint"

df_combined <- rbind(fleads,mleads)
df_combined <- rbind(df_combined,jleads)

# Remove chr X to analyse separately
df_combined <- df_combined %>% filter(chr!="23") %>% drop_na()

# Number of eQTLs by celltype in each group
df_combined %>% ungroup %>% 
    group_by(celltype, group) %>% 
    count() %>% spread(group, n) %>% data.frame()


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

# eQTL upset plots
icell=celltypes[21]
df_icell <- df_combined %>% filter(celltype==icell)
df_icell_wide <- df_icell %>% 
    group_by(geneid, snpid) %>%
    count(group) %>% 
    mutate(eqtl = 1, Total = sum(n)) %>%
        pivot_wider(
        id_cols = c(geneid, snpid, Total),
        names_from = group,
        values_from = eqtl,
        values_fill = list(eqtl = 0)) %>% ## fill empty cells with 0
        data.frame()

png(sprintf("/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/sex_specific_eQTL_analysis/results/Figures_and_Tables/eQTL_upsets/%s_upset_plot.png", icell),
    width = 12*300, height = 8*300,res = 300)   
    upset(df_icell_wide, 
    #sets = c("male", "female", "joint"),
    keep.order= TRUE,
    sets.bar.color = c("#2C788F", "#671D5B", "#F7A51D"),
    # order.by = "freq",
    point.size = 3.5,
    text.scale = c(2, 2, 2, 1.5, 2, 1.75),
    mainbar.y.label = "Number of cis associations",
    sets.x.label = "Total no. of eQTLs"
    ) 
    grid.text(icell, x = 0.65, y = 0.95, gp = gpar(fontsize = 20))   
dev.off()

# eGenes upset plots
icell=celltypes[21]
df_icell <- df_combined %>% filter(celltype==icell)
df_icell_wide <- df_icell %>% 
    group_by(geneid) %>%
    count(group) %>% 
    mutate(egenes = 1, Total = sum(n)) %>%
        pivot_wider(
        id_cols = c(geneid, Total),
        names_from = group,
        values_from = egenes,
        values_fill = list(egenes = 0)) %>% ## fill empty cells with 0
        data.frame()

png(sprintf("/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/sex_specific_eQTL_analysis/results/Figures_and_Tables/eGenes_upsets/%s_egene_upset_plot.png", icell),
    width = 12*300, height = 8*300,res = 300)   
    upset(df_icell_wide, 
    #sets = c("male", "female", "joint"),
    keep.order= TRUE,
    sets.bar.color = c("#2C788F", "#671D5B", "#F7A51D"),
    # order.by = "freq",
    point.size = 3.5,
    text.scale = c(2, 2, 2, 1.5, 2, 1.75),
    mainbar.y.label = "Number of eGenes",
    sets.x.label = "Total no. of eGenes"
    ) 
    grid.text(icell, x = 0.65, y = 0.95, gp = gpar(fontsize = 20))   
dev.off()







meGenes <- mleads %>% ungroup() %>% select("geneid") %>% unique()
meGenes$no1 <- '1'
meGenes$no1 <- as.numeric(meGenes$no1)
feGenes <- fleads %>% ungroup() %>% select("geneid") %>% unique()
feGenes$no10 <- '10'
feGenes$no10 <- as.numeric(feGenes$no10)
jeGenes <- jleads %>% select("GeneID") %>% unique()
colnames(jeGenes) <- "geneid"
jeGenes$no100 <- "100"
jeGenes$no100 <- as.numeric(jeGenes$no100)

all_genes <- full_join(meGenes, jeGenes)
all_genes <- full_join(all_genes, feGenes)
all_genes[is.na(all_genes)] <- 0

all_genes$sum <- as.factor (as.numeric(all_genes$no1) + as.numeric(all_genes$no10) + as.numeric(all_genes$no100))
