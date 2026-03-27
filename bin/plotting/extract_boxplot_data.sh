## SGE SETTINGS
#$ -S /bin/bash
#$ -wd '/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/sex_specific_eQTL_analysis'
#$ -o '/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/sex_specific_eQTL_analysis/stdout'
#$ -N extraction
#$ -q short.q
#$ -l mem_requested=48G,tmp_requested=48G,tmpfree=48G
#$ -r yes
#$ -j y

# debug
set -x

# Clear the environment
. /etc/profile.d/modules.sh

# Upload modules
module load /share/ClusterShare/Modules/modulefiles/contrib/evaben/gcc/gcc-7.3.0/7.3.0
R_PATH="/directflow/SCCGGroupShare/projects/SeyhanYazar/.conda/envs/r_env/bin/"

PARAMS="/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/sex_specific_eQTL_analysis/scripts/interacting_eQTLs_to_plot.txt"
# PARAMS="/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/sex_specific_eQTL_analysis/data/sex_specific_examples.txt"
# PARAMS="/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/sex_specific_eQTL_analysis/data/sex_specific_examples_f.txt"

rsid=`head -n $SGE_TASK_ID $PARAMS | tail -n 1 | awk '{print $2}'`
genename=`head -n $SGE_TASK_ID $PARAMS | tail -n 1 | awk '{print $1}'`

# Main script file
RSCRIPT="/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/sex_specific_eQTL_analysis/scripts/extract_boxplot_data.R"

# Do the main job
$R_PATH/Rscript --verbose ${RSCRIPT} ${genename} ${rsid} > rout/extracted_${genename}_${rsid}.Rout
