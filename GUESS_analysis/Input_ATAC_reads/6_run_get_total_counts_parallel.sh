#PBS -S /bin/bash
#PBS -N counts
#PBS -l walltime=8:00:00
#PBS -l select=1:ncpus=32:mem=300gb
#PBS -o ~/HRJ_monocytes/AS_ATAC/WASP/mono_34reps_noMultiMapped/OU
#PBS -e ~/HRJ_monocytes/AS_ATAC/WASP/mono_34reps_noMultiMapped/ER
#PBS -J 1-34
##########################################################

TESTNAME=mono_34reps_noMultiMapped
DIR="~/HRJ_monocytes/AS_ATAC/WASP"
BAMS="${DIR}/${TESTNAME}/unique"
BED_DIR="${DIR}/${TESTNAME}/beds"
PEAKS=${DIR}/${TESTNAME}/snps/closest_peaks_5kb.bed

##########################################################
cd ${BAMS}
INFILE=$(head -n $PBS_ARRAY_INDEX samples.list | tail -n 1)

source activate DT_DPLYR # need R and argparser as well as data.table and dplyr

cd ${DIR}/${TESTNAME}/peak_counts

Rscript ~/HRJ_monocytes/AS_ATAC/scripts/6_get_total_counts_for_parallel.R \
	${BED_DIR}/${INFILE}.bed $PEAKS

conda deactivate
