## SGE SETTINGS
#$ -S /bin/bash
#$ -wd '/directflow/SCCGGroupShare/projects/SeyhanYazar/covid19/smr_analysis_210906'
#$ -o '/directflow/SCCGGroupShare/projects/SeyhanYazar/covid19/smr_analysis_210906/stdout'
#$ -N smr_s1
#$ -q short.q
#$ -pe smp 4
#$ -l mem_requested=12G,tmp_requested=12G,tmpfree=12G
#$ -r yes
#$ -j y

# Master sample file
PARAMS="/directflow/SCCGGroupShare/projects/SeyhanYazar/covid19/matrix_eqtl_27_cells/scripts/all_cells_array"

CHRNUMBER=`head -n $SGE_TASK_ID $PARAMS | tail -n 1 | awk '{print $1}'`
CELLTYPE=`head -n $SGE_TASK_ID $PARAMS | tail -n 1 | awk '{print $2}'`

# CHRNUMBER="10"
# CELLTYPE="CD4all"

smr_tool="/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/tools/smr_Linux"
input_dir="/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/sex_specific_eQTL_analysis/results/2022-03-14_conditional-analysis-stratified/${CELLTYPE}"
output_dir="/directflow/SCCGGroupShare/projects/SeyhanYazar/covid19/smr_analysis_210906/results/besd-files/${CELLTYPE}"

mkdir $output_dir

$smr_tool --eqtl-summary ${input_dir}/${CELLTYPE}_chr${CHRNUMBER}_cis_eqtls.tsv --matrix-eqtl-format --make-besd --out ${output_dir}/${CELLTYPE}_chr${CHRNUMBER}