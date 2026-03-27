library(data.table)
library(dplyr)

setwd("/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/reference/CPG_ref/")

genes <- fread("gencode.v38.annotation.gtf.gz")
setnames(genes, names(genes), c("chr","source","type","start","end","score","strand","phase","attributes") )

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
genes$gene_type <- unlist(lapply(genes$attributes, extract_attributes, "gene_type"))
genes$gene_id2 <- sub("\\..*","", genes$gene_id)
genes <- genes %>% select(gene_id, gene_name, chr, start, end, strand)

setwd("/share/ScratchGeneral/vikgna/data/reference_data/refdata-cellranger-GRCh38/genes/")

problematic_genes <- fread("genes.gtf.gz")
setnames(problematic_genes, names(problematic_genes), c("chr","source","type","start","end","score","strand","phase","attributes") )
problematic_genes <- problematic_genes[type == "gene"]

problematic_genes$gene_id <- unlist(lapply(problematic_genes$attributes, extract_attributes, "gene_id"))
problematic_genes$gene_name <- unlist(lapply(problematic_genes$attributes, extract_attributes, "gene_name"))
problematic_genes <- problematic_genes %>% select(gene_id, gene_name, chr, start, end, strand)

problematic_genes$gene_id %in% genes$gene_id2

table(genes$gene_id2 %in% problematic_genes$gene_id)
genes2 <- genes %>% filter(genes$gene_id2 %in% problematic_genes$gene_id)
table(problematic_genes$gene_name %in% genes2$gene_name)