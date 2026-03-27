##############################################################################
# Script information   

# Title: Conditional cis-eQTL mapping - round 1
# Author: Seyhan Yazar
# Date: 2020-12-23
# Description: This bash script was written to run an array job  

##############################################################################

## SGE SETTINGS
#$ -S /bin/bash
#$ -wd '/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/sex_specific_eQTL_analysis'
#$ -o '/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/sex_specific_eQTL_analysis/stdout'
#$ -N nonpar_eqtl
#$ -l mem_requested=48G,tmp_requested=48G,tmpfree=48G
#$ -pe smp 8
#$ -r yes
#$ -j y

# cat /tmp/prolog_exec_"$JOB_ID"_"$SGE_TASK_ID".log

echo "JOB: $JOB_ID TASK: $SGE_TASK_ID"
echo "$HOSTNAME $tmp_requested $TMPDIR"

# debug
set -x

# Clear the environment
# . /etc/profile.d/modules.sh

# R Path
R_PATH="/directflow/SCCGGroupShare/projects/SeyhanYazar/.conda/envs/r_env/bin"

# load cell type specific array.txt file
CTYPEFILE="/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/sex_specific_eQTL_analysis/scripts/celltypes.txt" 

CELLLABEL=`head -n $SGE_TASK_ID $CTYPEFILE | tail -n 1 | awk '{print $1}'`

# Main script file
RSCRIPT="/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/sex_specific_eQTL_analysis/scripts/NonPAR_XaXiX_analysis_part1.R"

# Do the main job
${R_PATH}/Rscript ${RSCRIPT} ${CELLLABEL} > rout/round1C_${CELLLABEL}.Rout

# log file
LOG=logs

# check exit status
STATUS=$?
if [[ $STATUS -eq 0 ]]; then
     # success, write MD5 verification file
     echo -e `date` "\t${JOB_ID}\t${SGE_TASK_ID}\t${CELLLABEL}\tOK" >> $LOG
else
     echo -e `date` "\t${JOB_ID}\t${SGE_TASK_ID}\t${CELLLABEL}\tFAIL\t${STATUS}" >> $LOG
fi
exit $STATUS