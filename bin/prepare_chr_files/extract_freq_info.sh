## SGE SETTINGS
#$ -S /bin/bash
#$ -wd '/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/sex_specific_eQTL_analysis'
#$ -o '/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/sex_specific_eQTL_analysis/stdout'
#$ -N extraction
#$ -q short.q
#$ -l mem_requested=24G,tmp_requested=24G,tmpfree=24G
#$ -r yes
#$ -j y

# debug
set -x

# Clear the environment
. /etc/profile.d/modules.sh

TOOLDIR="/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/tools"
IMPUTEDIR="/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/imputed_data/decompressed/filter_vcf_r08_maf005"

# set the variables
CTYPEFILE="/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/cell_specific_eQTL_analysis_June21/scripts/celltypes.txt" 

CELLLABEL=`head -n $SGE_TASK_ID $CTYPEFILE | tail -n 1 | awk '{print $1}'`
# CELLLABEL="DC"

cd "./results/frequencies/${CELLLABEL}"

for file in *_snps.tsv; do
    # file="DC_chr13_snps.tsv"
    chr=`echo $file | awk -F '[_.]' '{print $2}'| sed 's/^chr//'`
    ${TOOLDIR}/plink --bfile ${IMPUTEDIR}/plink_chr${chr} --extract $file --freq --within ${CELLLABEL}_samples.tsv --out chr_${chr} ;
done
