# Upload libraries
library(data.table)
library(dplyr)

setwd("/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/cell_specific_eQTL_analysis_June21/results/2021-08-16_conditional-analysis/all_celltypes")

# Get all the result files (significant eQTls)
files <- list.files(path=".", pattern="*.tsv")
dataset <- lapply(files, function(x) { tryCatch(fread(x) , error=function(e) NULL)})
# Get file names and add to the list
filenames <- sub("*.tsv", "", files)
names(dataset) <- filenames

# Filter empty dataframes
dataset <- Filter(function(x) dim(x)[1] > 0, dataset)

dataset2 <- bind_rows(dataset, .id='filename')
dataset2$cell_type <- gsub("(.*?)\\s*_\\s*[^_]*$", "\\1", dataset2$filename)
dataset2$eSNP_rank <- gsub(".*_","",dataset2$filename)

dataset2 <- dataset2 %>% select(cell_type, geneid, snpid) %>% unique()
dataset2$chr <- gsub(':.*',"", dataset2$snpid)

cis_to_test <- fread("/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/sex_specific_eQTL_analysis/results/cis_eQTLs_to_test.tsv")
cis_to_test$chr <- as.character(cis_to_test$chr)

dataset3 <- anti_join(dataset2, cis_to_test)
fwrite(dataset3, "/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/sex_specific_eQTL_analysis/results/cis_eQTLs_to_test_additional.tsv", sep="\t")

cis_add <- fread("/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/cell_specific_eQTL_analysis_June21/results/2022-03-03_tidied_results/additional_eqtls")
cis_add$snpid <- paste0(cis_add$V5,":", cis_add$V6, "_", cis_add$V7)
cis_add_test <- cis_add %>% select("V1", "V2", "snpid", "V5")
colnames(cis_add_test) <- c("cell_type", "geneid", "snpid", "chr")
fwrite(cis_add_test, "/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/sex_specific_eQTL_analysis/results/cis_eQTLs_to_test_additional2.tsv", sep="\t")

### tidy ChrX results ###
setwd("/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/sex_specific_eQTL_analysis/results/2022-03-08_conditional-analysis/eSNPs")
files <- list.files(path=".", pattern="*.tsv")
dataset <- lapply(files, function(x) { tryCatch(fread(x) , error=function(e) NULL)})
# Get file names and add to the list
filenames <- sub("*.tsv", "", files)
names(dataset) <- filenames
dataset <- Filter(function(x) dim(x)[1] > 0, dataset)
dataset2 <- bind_rows(dataset,.id='filename')
dataset2$eSNP_rank <- gsub(".*_","",dataset2$filename)
dataset2$cell_type <- gsub("(.*?)\\s*_\\s*[^_]*$", "\\1", dataset2$filename) %>% gsub('*_chrX','',.)

# Read HRC ids and add to the dataframes
hrc_ids <- fread("/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/reference/HRC.r1-1.GRCh37.wgs.mac5.sites.tab")
head(hrc_ids)
hrc_ids$snpid <- sprintf("%s:%s_%s", hrc_ids$'#CHROM', hrc_ids$POS, hrc_ids$ALT)
hrc_ids <- hrc_ids %>% select(snpid, '#CHROM', POS, ID)
colnames(hrc_ids) <- c("ID", "CHR","POS", "snpid")

# Add the HRC ids
df_eqtl <- left_join(dataset2, hrc_ids, by="snpid")
nrow(df_eqtl)

# Gene names 
setwd("/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/CPG/gene_data/")

genes <- fread("gencode.v38.annotation.gtf")
setnames(genes, names(genes), c("chr","source","type","start","end","score","strand","phase","attributes") )

# [optional] focus, for example, only on entries of type "gene", 
# which will drastically reduce the file size
genes <- genes[type == "gene"]

extract_attributes <- function(gtf_attributes, att_of_interest){
  att <- strsplit(gtf_attributes, "; ")
  att <- gsub("\"","",unlist(att))
  if(!is.null(unlist(strsplit(att[grep(att_of_interest, att)], " ")))){
    return( unlist(strsplit(att[grep(att_of_interest, att)], " "))[2])
  }else{
    return(NA)}
}

genes$gene_id <- unlist(lapply(genes$attributes, extract_attributes, "gene_id"))
genes$gene_name <- unlist(lapply(genes$attributes, extract_attributes, "gene_name"))

genes <- genes %>% select(gene_id, gene_name, chr, start, end, strand)
my_genes <- genes %>% select(gene_id, gene_name)
colnames(my_genes) <- c("ensemblid", "geneid")

my_genes <- my_genes[!duplicated(my_genes$geneid),]

df_eqtl <- left_join(df_eqtl, my_genes, by="geneid")
str(df_eqtl) 

df_eqtl$`assessed allele` <- gsub(".*_", "", df_eqtl$ID)

# df_eqtl_final  <- df_eqtl %>% select(geneid, ensemblid, ID, `assessed allele`, new_names, estimate, statistic, p.value, qvalue, localFDR)
df_eqtl_final  <- df_eqtl %>% select(cell_type, geneid, ensemblid, ID, CHR, POS, `assessed allele`, eSNP_rank, estimate, statistic, p.value, qvalue, localFDR)
str(df_eqtl_final)

colnames (df_eqtl_final) <- c("GeneID", "Gene_EnsemblID", "SNP_rsID", "SNP_assessed_allele", "cell_type", "rho_correlation_coefficient", "s-statistics", "pvalue", "qvalue","FDR")
colnames (df_eqtl_final) <- c("cell_type", "GeneID", "Gene_EnsemblID", "rsID", "Chromosome", "Position", "SNP_assessed_allele",  "eSNP_rank", "rho_correlation_coefficient", "S-statistics","pvalue", "qvalue","FDR")
nrow(df_eqtl_final)
specify_decimal <- function(x) trimws(format(round(x, 3), nsmall=3))
# df_eqtl_final[,c("rho_correlation_coefficient","S-statistics")] <- apply(df_eqtl_final[,c("rho_correlation_coefficient","S-statistics")], 2, specify_decimal)
df_eqtl_final[,c("rho_correlation_coefficient")] <- apply(df_eqtl_final[,c("rho_correlation_coefficient")], 2, specify_decimal)
df_eqtl_final[,c("S-statistics")] <- apply(df_eqtl_final[,c("S-statistics")], 2, specify_decimal)

scientific_notation <- function (x) formatC(x, format = "e", digits = 3)
df_eqtl_final[,c("pvalue")] <- apply(df_eqtl_final[,c("pvalue")], 2,scientific_notation) 
df_eqtl_final[,c("qvalue")] <- apply(df_eqtl_final[,c("qvalue")], 2,scientific_notation) 
df_eqtl_final[,c("FDR")] <- apply(df_eqtl_final[,c("FDR")], 2,scientific_notation) 

df_eqtl_final <- df_eqtl_final %>% filter(!cell_type %in% c('pDC','cDC2','cDC'))

fwrite(df_eqtl_final, "/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/sex_specific_eQTL_analysis/results/chrX_cis_eQTLs_20220314_v2.tsv", sep="\t")