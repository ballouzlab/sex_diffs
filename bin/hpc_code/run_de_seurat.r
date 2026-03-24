args <- commandArgs(trailingOnly = TRUE)
i  <- as.numeric(args[1])
print(i) 
library(Seurat)
in_seurat_file = "data/cell_type.RDS"
in_indiv_phen_file = "data/phen_final.Rdata"
in_meta_file = "data/metadata.Rdata"
in_pcs_file = "onek1k_peer/pcs.Rdata"
load( in_meta_file) 
load( in_indiv_phen_file) 

freq = plyr::count( metadata_all$predicted.celltype.l3 ) 
celltypes = freq[,1]
celltypes = celltypes[-c(17,29)]


celltype = celltypes[i] 
    print(celltype)
    temp = readRDS(file=paste0(celltype, "_seurat.Rds"))
    temp <- NormalizeData(temp, normalization.method = "LogNormalize", scale.factor = 10000)
    temp <- FindVariableFeatures(temp, selection.method = "vst", nfeatures = 2000)
    all.genes <- rownames(temp)
    temp <- ScaleData(temp, features = all.genes, vars.to.regress="age")

    test = FindMarkers(temp,  slot="scale.data", ident.1="1", ident.2="2", logfc.threshold=-Inf, test.use="wilcox", group.by="sex")
    save(test, file = paste0("wilcox_DE_v4_", celltype, ".Rdata"))
    
    test = FindMarkers(temp,  slot="scale.data", ident.1="1", ident.2="2", logfc.threshold=-Inf, test.use="MAST", group.by="sex")
    save(test, file = paste0("mast_DE_v4_", celltype, ".Rdata"))




