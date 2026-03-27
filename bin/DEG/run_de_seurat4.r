args <- commandArgs(trailingOnly = TRUE)
i  <- as.numeric(args[1])
print(i) 
library(Seurat)
library(MAST)
in_seurat_file = "data/cell_type.RDS"
in_indiv_phen_file = "data/phen_final.Rdata"
in_meta_file = "data/metadata.Rdata"
in_pcs_file = "onek1k_peer/pcs.Rdata"
load( in_meta_file) 
load( in_indiv_phen_file) 
load( in_pcs_file) 
freq = plyr::count( metadata_all$predicted.celltype.l3 ) 
celltypes = freq[,1]
celltypes = celltypes[-c(17,29)]


celltype = celltypes[i] 
    print(celltype)
    temp = readRDS(file=paste0(celltype, "_seurat.Rds"))
    temp@active.assay = "SCT"
    test = FindMarkers(temp,  latent.var=c("age","sampleid", "pool_id"), verbose=F, ident.1="1", ident.2="2", logfc.threshold=-Inf, test.use="MAST", group.by="sex")
    save(test, file = paste0("mast_DE_v6_", celltype, ".Rdata"))
