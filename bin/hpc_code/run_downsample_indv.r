args <- commandArgs(trailingOnly = TRUE)
i  <- as.numeric(args[1])
n  <- as.numeric(args[2])
r  <- as.numeric(args[3])


print(c(i,n))

library(Seurat)
library(MAST) 

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
temp@active.assay = "SCT"
N = dim(temp)[2]

if( n > N ) {
 exit
} 
 

subset = list() 
test_list = list() 
for( i in 1:r){ 
   subset[[i]] = sample(N, n)
}

for( i in 1:r){ 
   sub = temp[, subset[[i]] ]
   print(i)
   test_list[[i]] = FindMarkers(sub, ident.1="1", ident.2="2", logfc.threshold=-Inf ,  group.by="sex", test.use = "wilcox", verbose=F)
}
save(test_list, subset, n, N, file=paste0(celltype,"_", n, "_downsample_Seuratv4.Rdata"))

 
