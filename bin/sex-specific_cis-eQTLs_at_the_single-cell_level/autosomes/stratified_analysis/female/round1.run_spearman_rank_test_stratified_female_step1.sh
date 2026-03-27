##############################################################################
# Script information   

# Title: Conditional cis-eQTL mapping - round 1
# Author: Seyhan Yazar
# Date: 2022-03
# Description: 

##############################################################################

## SGE SETTINGS
#$ -S /bin/bash
#$ -wd '/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/sex_specific_eQTL_analysis'
#$ -o '/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/sex_specific_eQTL_analysis/stdout'
#$ -N seqtl_f_s1
#$ -l mem_requested=24G,tmp_requested=24G,tmpfree=24G
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
RSCRIPT="/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/sex_specific_eQTL_analysis/scripts/round1.run_spearman_rank_test_stratified_female_step1.R"

# Do the main job
${R_PATH}/Rscript ${RSCRIPT} ${CELLLABEL} > rout/round1_${CELLLABEL}.Rout

# log file
LOG=logs/log_round1

# check exit status
STATUS=$?
if [[ $STATUS -eq 0 ]]; then
     # success, write MD5 verification file
     echo -e `date` "\t${JOB_ID}\t${SGE_TASK_ID}\t${CELLLABEL}\tOK" >> $LOG
else
     echo -e `date` "\t${JOB_ID}\t${SGE_TASK_ID}\t${CELLLABEL}\tFAIL\t${STATUS}" >> $LOG
fi
exit $STATUS