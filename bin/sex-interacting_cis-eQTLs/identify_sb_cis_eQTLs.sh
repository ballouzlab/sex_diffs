##############################################################################
# Script information   

# Title: Identify sex-biased cis-eQTLs (autosomal only)
# Author: Seyhan Yazar
# Date: 2020-12-23
# Description: This bash script was written to run an array job 
# per chromosome for 14 cell types using "identify_sex_biased_eQTLs" script. 

##############################################################################

## SGE SETTINGS
#$ -S /bin/bash
#$ -wd '/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/sex_specific_eQTL_analysis'
#$ -o '/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/sex_specific_eQTL_analysis/stdout'
#$ -N sb-cis-eqtl
#$ -l mem_requested=64G,tmp_requested=64G,tmpfree=64G
#$ -r yes
#$ -j y

cat /tmp/prolog_exec_"$JOB_ID"_"$SGE_TASK_ID".log
 
echo "JOB: $JOB_ID TASK: $SGE_TASK_ID"
echo "$HOSTNAME $tmp_requested $TMPDIR"

# debug
set -x

# Clear the environment
# . /etc/profile.d/modules.sh

# Chunks text files 
CHUNK_FILE="/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/sex_specific_eQTL_analysis/scripts/chunks.txt"

START=`head -n $SGE_TASK_ID $CHUNK_FILE | tail -n 1 | awk '{print $1}'`
END=`head -n $SGE_TASK_ID $CHUNK_FILE | tail -n 1 | awk '{print $2}'`

# R Path
R_PATH="/directflow/SCCGGroupShare/projects/SeyhanYazar/.conda/envs/r_env/bin"

R_SCRIPT_PATH="/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/sex_specific_eQTL_analysis/scripts/identify_sex_biased_cis_eQTLs.R" 

${R_PATH}/Rscript ${R_SCRIPT_PATH} ${START} ${END} > rout/sb_cis_eqtl_${START}.Rout

