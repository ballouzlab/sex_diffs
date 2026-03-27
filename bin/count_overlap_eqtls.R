library(data.table)
library(dplyr)

setwd("/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/sex_specific_eQTL_analysis/results/")

joint_df <- fread("2022-09-27_joint_analysis_one_round/joint_analysis_one_round_results.tsv")
joint_eGenes <- joint_df %>% select("cell_type","GeneID")
colnames(joint_eGenes) <- c("celltype","geneid")

female_df <- fread("2022-10-19_sex-specific-eQTLs-autosomes/sex_specific_female_eQTLs_not_in_joint_results_and_passed_both_tests.tsv")
female_specific_eGenes <- female_df %>% select("celltype","geneid")

male_df <- fread("2022-10-19_sex-specific-eQTLs-autosomes/sex_specific_male_eQTLs_not_in_joint_results_and_passed_both_tests.tsv")
male_specific_eGenes <- male_df %>% select("celltype","geneid")

celltype_list <- unique(joint_eGenes$celltype) 
celltype_list <- celltype_list[celltype_list!="HSPC"]

# Count of eQTLs

# 1) Joint - no female or male specific

celltype_specific_overlaps_1 <- function(c) {
   fm_df <- joint_df %>% filter(cell_type==c)
   f_df <- female_df %>% filter(celltype==c)
   m_df <- male_df %>% filter(celltype==c)
   new_df <- fm_df %>% filter(!fm_df$GeneID %in% f_df$geneid) 
   new_df2 <- new_df %>% filter(!new_df$GeneID %in% m_df$geneid)
   new_df2
}

list1 <- lapply(celltype_list,celltype_specific_overlaps_1)

# 2) Joint - female specific

celltype_specific_overlaps_2 <- function(c) {
   fm_df <- joint_df %>% filter(cell_type==c)
   f_df <- female_df %>% filter(celltype==c)
   m_df <- male_df %>% filter(celltype==c)
   new_df <- fm_df %>% filter(fm_df$GeneID %in% f_df$geneid)
   new_df2 <- new_df %>% filter(!new_df$GeneID %in% m_df$geneid)
   new_df2
}

list2 <- lapply(celltype_list,celltype_specific_overlaps_2)

# 3) Joint - male specific

celltype_specific_overlaps_3 <- function(c) {
   fm_df <- joint_df %>% filter(cell_type==c)
   m_df <- male_df %>% filter(celltype==c)
   f_df <- female_df %>% filter(celltype==c)
   new_df <- fm_df %>% filter(fm_df$GeneID %in% m_df$geneid)
   new_df2 <- new_df %>% filter(!new_df$GeneID %in% f_df$geneid)
   new_df2
}

list3 <- lapply(celltype_list,celltype_specific_overlaps_3)

# 4) No joint - female specific

celltype_specific_overlaps_4 <- function(c) {
   fm_df <- joint_df %>% filter(cell_type==c)
   f_df <- female_df %>% filter(celltype==c)
   m_df <- male_df %>% filter(celltype==c)
   new_df <- f_df %>% filter(!f_df$geneid %in% fm_df$GeneID)
   new_df2 <- new_df %>% filter(!new_df$geneid %in% m_df$geneid)
   new_df2
}

list4 <- lapply(celltype_list,celltype_specific_overlaps_4)

# 5) No joint - male specific

celltype_specific_overlaps_5 <- function(c) {
   fm_df <- joint_df %>% filter(cell_type==c)
   m_df <- male_df %>% filter(celltype==c)
   new_df <- m_df %>% filter(!m_df$geneid %in% fm_df$GeneID)
   new_df
}

list5 <- lapply(celltype_list,celltype_specific_overlaps_5)

### Count eQTLs in each condition 
count_eqtls_per_celltype <- function(x) {nrow(x)}

count_eqtls_per_list <- function (y) {
    # y=list5
    df <- lapply(y,count_eqtls_per_celltype)
    names(df) <- celltype_list
    bind_rows(df, .id="celltype")
}

group_list <- list(list1,list2,list3,list4,list5)
table <- lapply(group_list,count_eqtls_per_list) %>% bind_rows(., .id="group")
table$group <- factor(table$group, labels=c("Present in joint analysis only",
   "Present in joint and female-specific analysis",
   "Present in joint and male-specific analysis",
   "Present in female-specific analysis only",
   "Present in male-specific analysis only"))
fwrite(data.frame(table), "Figures_and_Tables/tables/eqtl_count_table.tsv", sep="\t")

summary_stats_tables <- function(x) { 
   list_df <- group_list[x]
   df <- bind_rows(list_df)
   fwrite(df, sprintf("Figures_and_Tables/tables/summary_stats_group%s.tsv", x), sep="\t")}
summary_stat_list <- lapply(1:5, summary_stats_tables)

sexint <- fread("/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/sex_specific_eQTL_analysis/results/2022-09-27_sex-biased-eQTL-analysis/sb_cis_eQTLs_at_25perc_fdr_20220928.tsv")


list1[[2]] %>% filter(GeneID=="TSPAN3")