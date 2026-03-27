##############################################################################
# Script information   

# Title: Matrix eQTL analysis
# Author: Seyhan Yazar
# Date: 2021-09-02
# Description:  

##############################################################################

## SGE SETTINGS
#$ -S /bin/bash
#$ -N meqtl_s2
#$ -q short.q
#$ -l mem_requested=48G
#$ -l tmp_requested=48G,tmpfree=48G
#$ -pe smp 2
#$ -r yes
#$ -j y
#$ -wd '/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/sex_specific_eQTL_analysis'
#$ -o '/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/sex_specific_eQTL_analysis/stdout'

# debug
set -x

# Clear the environment
. /etc/profile.d/modules.sh

# Upload the modules and paths
R_PATH="/directflow/SCCGGroupShare/projects/SeyhanYazar/.conda/envs/r_env/bin/"

# Master sample file
PARAMS="/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/sex_specific_eQTL_analysis/scripts/all_cells_array"

CHRNUMBER=`head -n $SGE_TASK_ID $PARAMS | tail -n 1 | awk '{print $1}'`
CELLTYPE=`head -n $SGE_TASK_ID $PARAMS | tail -n 1 | awk '{print $2}'`

# log file
LOG=logs/log_${CELLTYPE}

# Main script file
RSCRIPT="/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/sex_specific_eQTL_analysis/scripts/step2.run_MatrixEQTL_male.R"

# Do the main job
${R_PATH}/Rscript ${RSCRIPT} ${CELLTYPE} ${CHRNUMBER} > rout/step2.matrixeqtl_chr${CHRNUMBER}.Rout

# check exit status
STATUS=$?
if [[ $STATUS -eq 0 ]]; then
     # success, write MD5 verification file
     echo -e `date` "\t${JOB_ID}\t${SGE_TASK_ID}\t${CHRNUMBER}\t${CELLTYPE}\tOK" >> $LOG
else
     echo -e `date` "\t${JOB_ID}\t${SGE_TASK_ID}\t${CHRNUMBER}\tarrayjob.one\t${CELLTYPE}\tFAIL\t${STATUS}" >> $LOG
fi
exit $STATUS