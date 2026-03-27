## SGE SETTINGS
#$ -S /bin/bash
#$ -wd '/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/sex_specific_eQTL_analysis'
#$ -o '/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/sex_specific_eQTL_analysis/stdout'
#$ -N z-test
#$ -q short.q
#$ -pe smp 2
#$ -l mem_requested=24G,tmp_requested=24G,tmpfree=24G
#$ -r yes
#$ -j y

cat /tmp/prolog_exec_"$JOB_ID"_"$SGE_TASK_ID".log
 
echo "JOB: $JOB_ID TASK: $SGE_TASK_ID"
echo "$HOSTNAME $tmp_requested $TMPDIR"

# debug
set -x

# Clear the environment
. /etc/profile.d/modules.sh

# Upload modules
module load /share/ClusterShare/Modules/modulefiles/contrib/evaben/gcc/gcc-7.3.0/7.3.0
R_PATH="/directflow/SCCGGroupShare/projects/SeyhanYazar/.conda/envs/r_env/bin/"

# set the variables
CTYPEFILE="/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/sex_specific_eQTL_analysis/scripts/celltypes.txt" 

CELLLABEL=`head -n $SGE_TASK_ID $CTYPEFILE | tail -n 1 | awk '{print $1}'`

# Main script file
RSCRIPT="/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/sex_specific_eQTL_analysis/scripts/stats4_calculate_zstat_female.R"

# Do the main job
${R_PATH}/Rscript ${RSCRIPT} ${CELLLABEL} > rout/zstat_${CELLLABEL1}.Rout

# log file
LOG=logs/zstat

# check exit status
STATUS=$?
if [[ $STATUS -eq 0 ]]; then
     # success, write MD5 verification file
     echo -e `date` "\t${JOB_ID}\t${SGE_TASK_ID}\t${PERCENT}\tOK" >> $LOG
else
     echo -e `date` "\t${JOB_ID}\t${SGE_TASK_ID}\t${PERCENT}\tFAIL\t${STATUS}" >> $LOG
fi
exit $STATUS
