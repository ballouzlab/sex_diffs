## SGE SETTINGS
#$ -S /bin/bash
#$ -wd '/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/sex_specific_eQTL_analysis'
#$ -o '/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/sex_specific_eQTL_analysis/stdout'
#$ -N uniformity
#$ -q short.q
#$ -pe smp 2
#$ -l mem_requested=48G,tmp_requested=48G,tmpfree=48G
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

if [ ! -d results/pi_analysis/p_outputs ]
then
     mkdir results/pi_analysis/p_outputs
else
     echo "Directory exists"
fi

# Main script file
RSCRIPT="/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/sex_specific_eQTL_analysis/scripts/stats1_pvalue_uniformity_female.R"

# Do the main job
${R_PATH}/Rscript ${RSCRIPT} ${CELLLABEL} > rout/uniform_female_vs_joint_${CELLLABEL}.Rout

# log file
LOG=logs/uniform

# check exit status
STATUS=$?
if [[ $STATUS -eq 0 ]]; then
     # success, write MD5 verification file
     echo -e `date` "\t${JOB_ID}\t${SGE_TASK_ID}\t${PERCENT}\tOK" >> $LOG
else
     echo -e `date` "\t${JOB_ID}\t${SGE_TASK_ID}\t${PERCENT}\tFAIL\t${STATUS}" >> $LOG
fi
exit $STATUS
