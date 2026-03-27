library(data.table)
library(dplyr)
library(BEDMatrix)

setwd("/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/sex_specific_eQTL_analysis")

# Extract Genotypes from chr X
file <- "./data/X/chrX.filtered.bed"
bed_df <- BEDMatrix(file, simple_names=T)
chrMatrix <- as.matrix (bed_df)
sampleid <- gsub("^0_*","", rownames(chrMatrix))
rownames (chrMatrix) <- sampleid
chrMatrix <- data.frame(chrMatrix, check.names=F)
chrMatrix$sampleid <- rownames(chrMatrix)
genotype <- chrMatrix %>% select (sampleid, everything())
dim(genotype)
print(genotype[1:5,1:5])

fwrite(genotype, file='data/X/chrX_genotypes.tsv', quote=F,row.names=F, sep="\t")

par_region <- c("PAR1","PAR2")
# Extract Genotypes for PAR regions
for (i in 1:2) {
    p <- par_region[i]
    # Extract Genotypes from chr X
    file <- paste0("./data/XY/chrX_", p ,".bed")
    bed_df <- BEDMatrix(file, simple_names=T)
    chrMatrix <- as.matrix (bed_df)
    sampleid <- gsub("^0_*","", rownames(chrMatrix))
    rownames (chrMatrix) <- sampleid
    chrMatrix <- data.frame(chrMatrix, check.names=F)
    chrMatrix$sampleid <- rownames(chrMatrix)
    genotype <- chrMatrix %>% select (sampleid, everything())
    dim(genotype)
    print(genotype[1:5,1:5])
    write.table(genotype, file=paste0('data/XY/chrX_',p,'_genotypes.tsv'), quote=F,row.names=F, sep="\t")
}



