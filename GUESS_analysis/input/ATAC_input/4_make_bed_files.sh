#PBS -S /bin/bash
#PBS -N make_beds
#PBS -l walltime=12:00:00
#PBS -l select=1:ncpus=1:mem=12gb
#PBS -o ~/HRJ_monocytes/AS_ATAC/WASP/mono_34reps_noMultiMapped/OU
#PBS -e ~/HRJ_monocytes/AS_ATAC/WASP/mono_34reps_noMultiMapped/ER
#PBS -J 1-34

TESTNAME=mono_34reps_noMultiMapped

source activate WASP
DIR=~/HRJ_monocytes/AS_ATAC/WASP
bams=${DIR}/${TESTNAME}/unique

cd $bams
INFILE=$(head -n $PBS_ARRAY_INDEX samples.list | tail -n 1)

### Make bed files so that I can do total counts for peaks

## Get rid of X and Y and MT and make bed
cd $bams
bedtools bamtobed -i ${INFILE}.keep.merge.unique.sort.bam > ${DIR}/${TESTNAME}/beds/${INFILE}.bed.temp
cd ${DIR}/${TESTNAME}/beds
awk 'OFS="\t"{if (($1 != "X") && ($1 != "Y") && ($1 != "MT") && ($4 != "X") && ($4 != "Y") && ($4 != "MT")) print $0}' ${INFILE}.bed.temp > ${INFILE}.bed

rm ${INFILE}.bed.temp

conda deactivate
