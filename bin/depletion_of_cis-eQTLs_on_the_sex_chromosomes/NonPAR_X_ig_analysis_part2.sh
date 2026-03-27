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
#$ -N np_x_ig_eqtl
#$ -l mem_requested=64G,tmp_requested=64G,tmpfree=64G
#$ -r yes
#$ -j y
#$ -q short.q

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
RSCRIPT="/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/sex_specific_eQTL_analysis/scripts/NonPAR_X_ig_analysis_part2.R"

# Do the main job
${R_PATH}/Rscript ${RSCRIPT} ${CELLLABEL} > rout/round1C_${CELLLABEL}.Rout
