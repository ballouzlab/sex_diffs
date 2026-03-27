## SGE SETTINGS
#$ -S /bin/bash
#$ -wd '/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/science_revision_Dec20'
#$ -o '/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/science_revision_Dec20/stdout'
#$ -N plotting
#$ -q short.q
#$ -l mem_requested=24G,tmp_requested=24G,tmpfree=24G
#$ -r yes
#$ -j y

# debug
set -x

# Clear the environment
. /etc/profile.d/modules.sh

# Upload modules
module load /share/ClusterShare/Modules/modulefiles/contrib/evaben/gcc/gcc-7.3.0/7.3.0
R_PATH="/directflow/SCCGGroupShare/projects/SeyhanYazar/.conda/envs/onek1kEnv/bin/"

PARAMS="/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/Scripts/abf3041_revision/genotype_expression_plots/extract_boxplot_rsID_gene_pairs.txt"

genename=`head -n $SGE_TASK_ID $PARAMS | tail -n 1 | awk '{print $1}'`
rsid=`head -n $SGE_TASK_ID $PARAMS | tail -n 1 | awk '{print $2}'`

# Main script file
RSCRIPT="/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/Scripts/abf3041_revision/genotype_expression_plots/generate_boxplot.R"

# Do the main job
$R_PATH/Rscript --verbose ${RSCRIPT} ${genename} ${rsid} > rout/extracted_${genename}_${rsid}.Rout
