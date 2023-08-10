#PBS -S /bin/bash
#PBS -N combine_bams
#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=1:mem=62gb
#PBS -o ~/HRJ_monocytes/CHiC/hicup_combinations/OU
#PBS -e ~/HRJ_monocytes/CHiC/hicup_combinations/ER

##########################################################
#######  CHECK THAT FORMAT OF NAMES FOR HICUP CAPTURED BAM FILES IS CORRECT
##########################################################

SAMPLES=~/HRJ_monocytes/CHiC/hicup_combinations
CHICAGO_BAM=~/HRJ_monocytes/CHiC/chicago/combinations_input/bam-files

module load bowtie2/2.1.0
module load samtools/1.3.1
module load anaconda3/personal

source activate Renv1

# This HiCUP path has the comb_dedup script
PATH=$PATH:~/spivakov/bin/HiCUP/

# Combine and deduplicate for CHiCAGO, putting into input dir
cd $SAMPLES
mysamples=${SAMPLES}/combine_samples.list

for INDIVIDUAL in $(cat $mysamples)
do
    echo $INDIVIDUAL

mkdir -p ${SAMPLES}/${INDIVIDUAL}_combined
cd ${SAMPLES}/${INDIVIDUAL}_combined
ln -s ../${INDIVIDUAL}_tech1/${INDIVIDUAL}_tech1_R1_2.fastq.gz.permuted.pair.bam.allocated.bam.prefiltered.hicup.captured.bam ./
ln -s ../${INDIVIDUAL}_tech2/${INDIVIDUAL}_tech2_R1_2.fastq.gz.permuted.pair.bam.allocated.bam.prefiltered.hicup.captured.bam ./

comb_dedup.pl \
  ${INDIVIDUAL}_tech1_R1_2.fastq.gz.permuted.pair.bam.allocated.bam.prefiltered.hicup.captured.bam \
  ${INDIVIDUAL}_tech2_R1_2.fastq.gz.permuted.pair.bam.allocated.bam.prefiltered.hicup.captured.bam

cd ${CHICAGO_BAM}
ln -s ${SAMPLES}/${INDIVIDUAL}_combined/*.dedup.bam ./${INDIVIDUAL}.hicup.captured.bam

done

# Link all other bams to input dir
cd $SAMPLES
mysamples_noncombine=${SAMPLES}/combine_samples2.list

for INDIVIDUAL in $(cat $mysamples_noncombine)
do
    echo $INDIVIDUAL

cd ${CHICAGO_BAM}
ln -s ${SAMPLES}/${INDIVIDUAL}/${INDIVIDUAL}_R1_2.fastq.gz.permuted.pair.bam.allocated.bam.prefiltered.hicup.captured.bam ./${INDIVIDUAL}.hicup.captured.bam

done


conda deactivate












