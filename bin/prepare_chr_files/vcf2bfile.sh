## SGE SETTINGS
#$ -cwd
#$ -S /bin/bash
#$ -q short.q
#$ -r yes
#$ -l mem_requested=8G
#$ -l tmp_requested=8G
#$ -N vcf2bfile

cd $SGE_0_WORKDIR
BIN='/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/tools'
export PATH=$PATH:$BIN

dir='/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/sex_specific_eQTL_analysis/data/'

vcf_file='/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/sex_specific_eQTL_analysis/data/fixed_2022_10_10/chrX.dose.vcf.gz'

mkdir $dir/data

plink2 --vcf $vcf_file --make-bed --keep-allele-order --out $dir/chrX_20221010

# plink2 --bfile onek1k_chr23 --freq
