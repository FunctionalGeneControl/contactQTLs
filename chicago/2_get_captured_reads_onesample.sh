#PBS -S /bin/bash
#PBS -N get_captured_reads
#PBS -l walltime=12:00:00
#PBS -l select=1:ncpus=2:mem=12gb
#PBS -o ~/HRJ_monocytes/CHiC/hicup_combinations/OU
#PBS -e ~/HRJ_monocytes/CHiC/hicup_combinations/ER

##########################################################

#### CHANGE THIS
SAMPLES=~/HRJ_monocytes/CHiC/hicup_combinations
cd ${SAMPLES}
mysample=S026D1-21
####

module load samtools/1.3.1
module load anaconda3/personal
source activate Renv1
BAITS=~/HRJ_monocytes/CHiC/capture_efficiency/eCHiC_grch38_sorted.txt

cd ~/bin/HiCUP_combinations

perl ./Misc/get_captured_reads --baits ${BAITS} ${SAMPLES}/${mysample}/${mysample}_R1_2.fastq.gz.permuted.pair.bam.allocated.bam.prefiltered.hicup.bam
perl ./Misc/get_captured_reads --baits ${BAITS} ${SAMPLES}/${mysample}/${mysample}_R1_2.fastq.gz.permuted.pair.bam.allocated.bam.prefiltered.filt.bam

conda deactivate












