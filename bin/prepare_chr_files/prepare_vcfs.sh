##SGE SETTINGS
#$ -S /bin/bash
#$ -wd '/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/sex_specific_eQTL_analysis'
#$ -o '/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/sex_specific_eQTL_analysis/stdout'
#$ -q short.q
#$ -r yes
#$ -l mem_requested=128G,tmp_requested=128G,tmpfree=128G
#$ -N prep_vcfs
#$ -pe smp 16

export PATH=$PATH:/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/tools
export PATH=$PATH:/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/tools/bcftools

dataDir="/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/sex_specific_eQTL_analysis/data"

# plink2 --vcf ${dataDir}/X/chrX.dose.vcf --const-fid --maf 0.05 --exclude-if-info 'R2<0.8' --make-bed --out ${dataDir}/X/chrX
# plink2 --vcf ${dataDir}/XY/chrX_PAR1.dose.vcf --maf 0.05 --exclude-if-info 'R2<0.8' --const-fid --make-bed --out ${dataDir}/XY/chrX_PAR1 
# plink2 --vcf ${dataDir}/XY/chrX_PAR2.dose.vcf --const-fid --maf 0.05 --exclude-if-info 'R2<0.8' --make-bed --out ${dataDir}/XY/chrX_PAR2 

bcftools view -i 'R2>.8 & MAF>.05' -Oz ${dataDir}/fixed_2022_10_10/chrX.dose.vcf.gz > ${dataDir}/X/chrX.filtered.vcf.gz

plink2 --vcf ${dataDir}/X/chrX.filtered.vcf.gz --const-fid --make-bed --out ${dataDir}/X/chrX.filtered