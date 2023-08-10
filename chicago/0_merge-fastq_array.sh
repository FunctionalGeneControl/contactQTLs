#PBS -S /bin/bash
#PBS -N merge_fastq
#PBS -l walltime=48:00:00
#PBS -l select=1:ncpus=8:mem=64gb
#PBS -o ~/HRJ_monocytes/CHiC/hicup_combinations/OU
#PBS -e ~/HRJ_monocytes/CHiC/hicup_combinations/ER
#PBS -J 1-16

BASE=~/HRJ_monocytes/CHiC/hicup_combinations
cd $BASE
mysample=$(head -n ${PBS_ARRAY_INDEX} combine_samples2.list | tail -n 1)

cd ${BASE}/${mysample}
cat ${mysample}_A_R1.fastq.gz  ${mysample}_B_R1.fastq.gz > ${mysample}_R1.fastq.gz
cat ${mysample}_A_R2.fastq.gz  ${mysample}_B_R2.fastq.gz > ${mysample}_R2.fastq.gz

## Actually note that for S026D1-21 there were actually three sets of fastq files. I'm now adding it into the non-array script.
