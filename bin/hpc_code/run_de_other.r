args <- commandArgs(trailingOnly = TRUE)
i  <- as.numeric(args[1])
print(i) 
library(Seurat)
library(DESeq2) 

in_indiv_phen_file = "data/phen_final.Rdata"
in_meta_file = "data/metadata.Rdata"
in_pcs_file = "onek1k_peer/pcs.Rdata"
load(in_meta_file) 
load(in_indiv_phen_file) 
load(in_pcs_file) 


freq = plyr::count( metadata_all$predicted.celltype.l3 ) 
celltypes = freq[,1]
celltypes = celltypes[-c(17,29)]


   celltype = celltypes[i] 
   print(celltype)
   temp = readRDS(file=paste0(celltype, "_seurat.Rds"))
   pseudo_temp <- AggregateExpression(temp, assays = "RNA", return.seurat = T, group.by = c("individual", "sex"))
   pseudo_temp$individual = gsub("g", "", pseudo_temp$individual)
   pseudo_temp$individual = gsub("-", "_", pseudo_temp$individual)
   m = match(pseudo_temp$individual, phens$sampleid)

   pseudo_temp[["age"]] = phens$age
   pseudo_temp[["pool_id"]] = phens$pool_id

   m = match(pseudo_temp$individual, pcs$sampleid)

   pseudo_temp@meta.data = cbind(pseudo_temp@meta.data, pcs[m,3:8])
#   bulk.temp.de <- FindMarkers(object = pseudo_temp, ident.1="1", ident.2="2",  group.by="sex", test.use = "DESeq2")
   bulk.temp.wilcox <- FindMarkers(object = pseudo_temp, ident.1="1", ident.2="2",  group.by="sex", test.use = "wilcox")


   data  = pseudo_temp@assays$RNA$counts
   coldata  = pseudo_temp@meta.data[,2:11]

   # dds <- DESeqDataSetFromMatrix(countData=data, colData=coldata, design= ~ age + pc1 + pc2 + pc3 + pc4 + pool_id + sex)
   dds <- DESeqDataSetFromMatrix(countData=data, colData=coldata, design= ~ age +  sex)
   dds <- DESeq(dds)
   # res_age <- results(dds, contrast=c("sex","1","2"))
   res_age <- results(dds, contrast=c("sex","1","2"))

   # dds <- DESeqDataSetFromMatrix(countData=data, colData=coldata, design= ~ sex)
   # dds <- DESeq(dds)
   # res_sex <- results(dds)

   
   # save(bulk.temp.de , coldata, res, res_sex, file=paste0(celltype, "_pseudo_deseq2.Rdata") )
   save(bulk.temp.wilcox , coldata, res_age,   file=paste0(celltype, "_pseudo_wilcox.Rdata") )

