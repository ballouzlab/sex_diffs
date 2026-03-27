# Upload libraries
library(data.table)
library(dplyr)

setwd("/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/cell_specific_eQTL_analysis_June21/results/2021-08-16_conditional-analysis/all_celltypes")

# Get all the result files (significant eQTls)
files <- list.files(path=".", pattern="*_eSNP1.tsv")

dataset <- lapply(files[c(3:14,18:25,27,29)], function(x) { tryCatch(fread(x) , error=function(e) NULL)})
# Get file names and add to the list
filenames <- sub("*.tsv", "", files[c(3:14,18:25,27,29)])
names(dataset) <- filenames

# Filter empty dataframes
dataset <- Filter(function(x) dim(x)[1] > 0, dataset)
dataset <- bind_rows(dataset, .id='filename')
dataset$cell_type <- gsub("(.*?)\\s*_\\s*[^_]*$", "\\1", dataset$filename)

dataset2 <- dataset %>% select(cell_type, geneid, snpid) %>% unique()
dataset2$chr <- gsub(':.*',"", dataset2$snpid)

fwrite(dataset2, "/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/sex_specific_eQTL_analysis/results/2022-09-27_joint_analysis_one_round/cis_assoc_to_test_for_sex_bias.tsv", sep="\t")

# Read HRC ids and add to the dataframes
hrc_ids <- fread("/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/reference/HRC.r1-1.GRCh37.wgs.mac5.sites.tab")
head(hrc_ids)
hrc_ids$snpid <- sprintf("%s:%s_%s", hrc_ids$'#CHROM', hrc_ids$POS, hrc_ids$ALT)
hrc_ids <- hrc_ids %>% select(snpid, '#CHROM', POS, ID)

# Add the HRC ids
df_eqtl <- left_join(dataset, hrc_ids, by="snpid")
nrow(df_eqtl)

# Read MAF files to extract allele info 
maf_files <- list.files(path="/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/imputed_data/decompressed/snps_with_maf_greaterthan0.05", 
    pattern=".txt", full.name=TRUE)
mafs_lst <- lapply(maf_files, function(x) fread(x))
mafs_df <- bind_rows(mafs_lst)
mafs_df$snpid <- sprintf("%s_%s", mafs_df$SNP, mafs_df$'ALT(1)')
mafs_df <- mafs_df %>% select(snpid, 'REF(0)','ALT(1)')
colnames(mafs_df) <- c("snpid", "other allele", "assessed allele")

# Add the allele info 
df_eqtl <- left_join(df_eqtl, mafs_df, by="snpid")

# Gene names 
genes <- fread("/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/reference/genes.tsv")
genes <- genes %>% select(Geneid, GeneSymbol)
colnames(genes) <- c("ensemblid", "geneid")
genes <- genes[!duplicated(genes$geneid),]

df_eqtl <- left_join(df_eqtl, genes, by="geneid")
str(df_eqtl) 

# df_eqtl_final  <- df_eqtl %>% select(geneid, ensemblid, ID, `assessed allele`, new_names, estimate, statistic, p.value, qvalue, localFDR)
df_eqtl_final  <- df_eqtl %>% select(cell_type, geneid, ensemblid, ID, `#CHROM`, POS, `assessed allele`, estimate, statistic, p.value, qvalue, localFDR)
str(df_eqtl_final)
colnames (df_eqtl_final) <- c("cell_type", "GeneID", "Gene_EnsemblID", "rsID", "Chromosome", "Position", "SNP_assessed_allele", "rho_correlation_coefficient", "S-statistics","pvalue", "qvalue","FDR")
nrow(df_eqtl_final)

specify_decimal <- function(x) trimws(format(round(x, 3), nsmall=3))
df_eqtl_final[,c("rho_correlation_coefficient")] <- apply(df_eqtl_final[,c("rho_correlation_coefficient")], 2, specify_decimal)
df_eqtl_final[,c("S-statistics")] <- apply(df_eqtl_final[,c("S-statistics")], 2, specify_decimal)

scientific_notation <- function (x) formatC(x, format = "e", digits = 3)
df_eqtl_final[,c("pvalue")] <- apply(df_eqtl_final[,c("pvalue")], 2,scientific_notation) 
df_eqtl_final[,c("qvalue")] <- apply(df_eqtl_final[,c("qvalue")], 2,scientific_notation) 
df_eqtl_final[,c("FDR")] <- apply(df_eqtl_final[,c("FDR")], 2,scientific_notation) 

fwrite(df_eqtl_final, "/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/sex_specific_eQTL_analysis/results/2022-09-27_joint_analysis_one_round/joint_analysis_one_round_results.tsv", sep="\t", quote=FALSE)

