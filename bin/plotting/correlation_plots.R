#!/usr/bin/env Rscript

# Upload libraries
library(data.table)
library(dplyr)
library(ggpubr)
library(hrbrthemes)

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

female_eqtls <- fread("2022-10-19_sex-specific-eQTLs-autosomes/sex_specific_female_eQTLs_not_in_joint_results_and_passed_both_tests.tsv")
female_eqtls$group <- "female"

correlations <- function (x) { 
    #x=celltypes[13]
    setwd("/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/sex_specific_eQTL_analysis/results")
    female_results_df <- female_eqtls %>% filter(celltype==x)
    dim(female_results_df)

    male_results_df <- fread(sprintf("/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/sex_specific_eQTL_analysis/results/2022-10-18_pi-analysis/p_outputs/female_vs_male_%s_pvalues.tsv", x))
    
    not_in_male_df <- female_results_df[female_results_df$snpid %in% male_results_df$snpid,]
    dim(not_in_male_df)

    female_df <- not_in_male_df %>% 
        select("geneid","snpid","estimate","localFDR") %>% 
        rename(., f.estimate=estimate, f.FDR=localFDR)
    male_df <-  male_results_df %>% 
        select("geneid","snpid","estimate","localFDR") %>% 
        rename(., m.estimate=estimate, m.FDR=localFDR)
    
    df <-left_join(female_df, male_df, by=c("geneid", "snpid"))
    df  
}

corr_data <- lapply(celltypes,correlations) 
names(corr_data) <- celltypes
corr_data <- bind_rows(corr_data, .id="celltype")

load("../data/palette2.Rdata")
colpals2$L2_cells <-  sub(" ", "_", colpals2$V9)
colpal <- colpals2 %>% filter(L2_cells %in% celltypes)
colnames(colpal)[10:11]<- c("colour", "celltype")
colpal <- colpal[,10:11]

corr_data <- left_join(corr_data,colpal, by="celltype")
colors <- distinct(corr_data, celltype, colour)
pal <- colors$colour
names(pal) <- colors$celltype
pal


pdf(paste0("Figures_and_Tables/correlations/female_specific_estimate_corr_plot.pdf"), width=10, height=7)
png(paste0("Figures_and_Tables/correlations/female_specific_estimate_corr_plot.png"), width=10, height=7,res=300, units='in')
ggplot(corr_data, aes(x=f.estimate, y=m.estimate)) +
        geom_point(aes(size = f.FDR, colour = celltype), alpha=0.5) +
        geom_text(aes(label=ifelse(f.estimate>0.3,geneid,'')), size=3) +
        geom_text(aes(label=ifelse(f.estimate<(-0.3),geneid,'')), size=3) +
        scale_color_manual(values=pal, breaks=colpal$celltype) +
        scale_size(trans = 'reverse') +
        geom_vline(xintercept = 0 ) +
        geom_hline(yintercept = 0) +
        theme_minimal () +
        xlab("female rho estimate") + xlim(-0.45,0.45) +
        ylab("male rho estimate") + ylim(-0.3, 0.3)
dev.off()

rm(list=ls())

### Male specific

# Cell types
celltypes <- c("B_intermediate", "B_memory", "B_naive", 
    "CD14_Mono", "CD16_Mono", 
    "CD4_CTL", "CD4_Naive", "CD4_TCM","CD4_TEM",
    "CD8_Naive", "CD8_TCM", "CD8_TEM", "DC",
    "dnT","gdT", "MAIT", "Treg",
    "NK", "NK_CD56bright", "NK_Proliferating", 
    "Plasmablast"
    )

male_eqtls <- fread("2022-10-19_sex-specific-eQTLs-autosomes/sex_specific_male_eQTLs_not_in_joint_results_and_passed_both_tests.tsv")
male_eqtls$group <- "male"

correlations <- function (x) { 
    #x=celltypes[13]
    setwd("/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/sex_specific_eQTL_analysis/results")
    male_results_df <- male_eqtls %>% filter(celltype==x)
    dim(male_results_df)

    female_results_df <- fread(sprintf("/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/sex_specific_eQTL_analysis/results/2022-10-18_pi-analysis/p_outputs/male_vs_female_%s_pvalues.tsv", x))
    
    not_in_female_df <- male_results_df[male_results_df$snpid %in% female_results_df$snpid,]
    dim(not_in_female_df)

    male_df <- not_in_female_df %>% 
        select("geneid","snpid","estimate","localFDR") %>% 
        rename(., m.estimate=estimate, m.FDR=localFDR)
    female_df <-  female_results_df %>% 
        select("geneid","snpid","estimate","localFDR") %>% 
        rename(., f.estimate=estimate, f.FDR=localFDR)
    
    df <-left_join(male_df, female_df, by=c("geneid", "snpid"))
    df  
}

corr_data <- lapply(celltypes,correlations) 
names(corr_data) <- celltypes
corr_data <- bind_rows(corr_data, .id="celltype")

load("../data/palette2.Rdata")
colpals2$L2_cells <-  sub(" ", "_", colpals2$V9)
colpal <- colpals2 %>% filter(L2_cells %in% celltypes)
colnames(colpal)[10:11]<- c("colour", "celltype")
colpal <- colpal[,10:11]

corr_data <- left_join(corr_data,colpal, by="celltype")
colors <- distinct(corr_data, celltype, colour)
pal <- colors$colour
names(pal) <- colors$celltype
pal

png(paste0("Figures_and_Tables/correlations/male_specific_estimate_corr_plot.png"), width=10, height=7,res=300, units='in')
pdf(paste0("Figures_and_Tables/correlations/male_specific_estimate_corr_plot.pdf"), width=10, height=7)
ggplot(corr_data, aes(x=m.estimate, y=f.estimate)) +
        geom_point(aes(size = m.FDR, colour = celltype), alpha=0.5) +
        geom_text(aes(label=ifelse(m.estimate>0.3,geneid,'')), size=3) +
        geom_text(aes(label=ifelse(m.estimate<(-0.3),geneid,'')), size=3) +
        scale_color_manual(values=pal, breaks=colpal$celltype) +
        scale_size(trans = 'reverse') +
        geom_vline(xintercept = 0) +
        geom_hline(yintercept = 0) +
        theme_minimal () +
        xlab("male rho estimate") +
        ylab("female rho estimate") + 
        xlim(-0.45,0.45) +
        ylim(-0.3, 0.3)
dev.off()


## MAF correlations
    setwd(sprintf("/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/sex_specific_eQTL_analysis/results/frequencies/%s", x))
    freq_list <-  list.files(path = ".", pattern = "*.frq.strat")
    freq_df <-lapply(freq_list, function(x) { tryCatch(fread(x) , error=function(e) NULL)})
    freq_df <- bind_rows(freq_df)
    maf_plot <- freq_df %>% 
        mutate(., SNPID=paste(SNP,A1,sep="_")) %>%
        select("SNPID","CLST","MAF") %>% 
        tidyr::pivot_wider(names_from = CLST, values_from = MAF) %>%
        ggplot(., aes(x=F, y=M)) +
        geom_point(shape=23) +
        theme_minimal()
    
    plots <- ggarrange (est_plot, pval_plot, maf_plot, ncol=3, nrow=1,
        labels = c("A", "B", "C"))

    png(sprintf("/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/sex_specific_eQTL_analysis/results/Figures_and_Tables/correlations/%s_plots.png",x), height=600, width=1200, res=150) # Open a new pdf file
    print(annotate_figure(plots, top = text_grob(print(x), 
               color = "purple", face = "bold", size = 14)))
    dev.off() # Close the file

