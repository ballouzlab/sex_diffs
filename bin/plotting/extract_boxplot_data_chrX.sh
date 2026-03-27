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

PARAMS="/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/sex_specific_eQTL_analysis/results/2022-12-08_expr-geno-plots/plot_examples.txt"

cell=`head -n $SGE_TASK_ID $PARAMS | tail -n 1 | awk '{print $1}'`
genename=`head -n $SGE_TASK_ID $PARAMS | tail -n 1 | awk '{print $2}'`
snp=`head -n $SGE_TASK_ID $PARAMS | tail -n 1 | awk '{print $3}'`

# Main script file
RSCRIPT="/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/sex_specific_eQTL_analysis/scripts/extract_boxplot_data_chrX.R"

# Do the main job
$R_PATH/Rscript --verbose ${RSCRIPT} ${cell} ${genename} ${snp} > rout/extracted_${cell}_${genename}_${snp}.Rout
