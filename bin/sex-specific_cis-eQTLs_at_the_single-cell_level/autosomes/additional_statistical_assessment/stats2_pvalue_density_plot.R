#!/usr/bin/env Rscript

# Upload libraries
library(data.table)
library(dplyr)
library(ggplot2)
library(qvalue)
library(hrbrthemes)
library(stringr)

# Cell types
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

#### Histogram of p-values ####
setwd("/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/sex_specific_eQTL_analysis/results/2022-10-18_pi-analysis/p_outputs")
ExpFileList <- list.files(path = ".", pattern = "_pvalues.tsv")

#### For female vs male (ie not eQTL in males) #####

calc_pi0 <- function (dfile) {
    pi_df <- fread(dfile)
    if (nrow(pi_df) <= 20) {
       stop ("too small") 
       } else { #continue the script
       print("Script did NOT end!")   
    }
    dname <- sub("*.pvalues.tsv", "", dfile)
    celltype <- sub("female_vs_male*.","", dname)
    qobj <- qvalue(pi_df$p.value,)
   
    piplot <- ggplot(pi_df, aes(x=p.value)) +
        geom_histogram(aes(y=..density..), fill="white", color="black") +
        theme_ipsum() +
        ggtitle(sprintf("female vs male %s", celltype)) +
        geom_hline(aes(yintercept=mean(qobj$pi0.lambda)), color="dark red",
             linetype="dashed") 
    setEPS()
    postscript(sprintf("../histograms/female_vs_male_%s_histogram.eps", str_replace(celltype, " ", "")))
    print(piplot)
    dev.off()
    pi0=mean(qobj$pi0.lambda)
    data.frame(celltype, pi0)
}

pi0_results <- lapply(ExpFileList[1:21], function(x) tryCatch(calc_pi0(x), error=function(e) NULL))
pi0_results_df <- bind_rows(pi0_results)

pi_list <- lapply(ExpFileList[1:21], function(x) { tryCatch(fread(x) , error=function(e) NULL)})
pi_list_df <- bind_rows(pi_list)
qobj_all <- qvalue(pi_list_df$p.value)
piplot_all <- ggplot(pi_list_df, aes(x=p.value)) +
        geom_histogram(aes(y=..density..), fill="white", color="black") +
        theme_ipsum() +
        ggtitle(sprintf("female vs male")) +
        geom_hline(aes(yintercept=mean(qobj_all$pi0.lambda)), color="dark red",
             linetype="dashed") 
    setEPS()
    postscript(sprintf("../histograms/female_vs_male_all_histogram.eps"))
    print(piplot_all)
    dev.off()
pi0_all=mean(qobj_all$pi0.lambda)
pi0_results_df <- pi0_results_df %>% add_row(celltype="All", pi0=pi0_all)

fwrite(pi0_results_df, "/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/sex_specific_eQTL_analysis/results/2022-10-18_pi-analysis/female_vs_male_pi0_values.tsv", sep="\t")

#### For male vs female (ie not eQTL in females) #####

calc_pi0 <- function (dfile) {
    pi_df <- fread(dfile)
    # if (nrow(pi_df) <= 20) {
    #   stop ("too small") 
    #   } else { #continue the script
    #   print("Script did NOT end!")   
    # }
    dname <- sub("*.pvalues.tsv", "", dfile)
    celltype <- sub("male_vs_female*.","", dname)
    qobj <- qvalue(pi_df$p.value)
   
    piplot <- ggplot(pi_df, aes(x=p.value)) +
        geom_histogram(aes(y=..density..), fill="white", color="black") +
        theme_ipsum() +
        ggtitle(sprintf("male vs female %s", celltype)) +
        geom_hline(aes(yintercept=mean(qobj$pi0.lambda)), color="dark red",
             linetype="dashed") 
    setEPS()
    postscript(sprintf("../histograms/male_vs_female_%s_histogram.eps", str_replace(celltype, " ", "")))
    print(piplot)
    dev.off()
    pi0=mean(qobj$pi0.lambda)
    data.frame(celltype, pi0)
}

pi0_results <- lapply(ExpFileList[22:42], function(x) tryCatch(calc_pi0(x), error=function(e) NULL))
pi0_results_df <- bind_rows(pi0_results)

pi_list <- lapply(ExpFileList[22:42], function(x) { tryCatch(fread(x) , error=function(e) NULL)})
pi_list_df <- bind_rows(pi_list)
qobj_all <- qvalue(pi_list_df$p.value)
piplot_all <- ggplot(pi_list_df, aes(x=p.value)) +
        geom_histogram(aes(y=..density..), fill="white", color="black") +
        theme_ipsum() +
        ggtitle(sprintf("female vs male")) +
        geom_hline(aes(yintercept=mean(qobj_all$pi0.lambda)), color="dark red",
             linetype="dashed") 
    setEPS()
    postscript(sprintf("../histograms/male_vs_female_all_histogram.eps"))
    print(piplot_all)
    dev.off()
pi0_all=mean(qobj_all$pi0.lambda)
pi0_results_df <- pi0_results_df %>% add_row(celltype="All", pi0=pi0_all)

fwrite(pi0_results_df, "/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/sex_specific_eQTL_analysis/results/2022-10-18_pi-analysis/male_vs_female_pi0_values.tsv", sep="\t")

