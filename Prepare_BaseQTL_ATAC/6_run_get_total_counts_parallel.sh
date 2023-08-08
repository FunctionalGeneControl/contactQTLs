#PBS -S /bin/bash
#PBS -N counts
#PBS -l walltime=8:00:00
#PBS -l select=1:ncpus=32:mem=300gb
#PBS -o ~/HRJ_monocytes/AS_ATAC/BaseQTL/OU
#PBS -e ~/HRJ_monocytes/AS_ATAC/BaseQTL/ER
#PBS -J 1-3
##########################################################

DIR="~/HRJ_monocytes/AS_ATAC/BaseQTL"
BAMS="${DIR}/bams_dedup"
BED_DIR="${DIR}/bed"
PEAKS=${DIR}/snps/closest_peaks_5kb.bed

##########################################################
cd ${BAMS}
INFILE=$(head -n $PBS_ARRAY_INDEX redo.txt | tail -n 1)

source activate DT_DPLYR # need R and argparser as well as data.table and dplyr

cd ${DIR}/peak_counts

Rscript ~/HRJ_monocytes/AS_ATAC/BaseQTL/scripts/6_get_total_counts_for_parallel.R \
	${BED_DIR}/${INFILE}.bed $PEAKS

conda deactivate
