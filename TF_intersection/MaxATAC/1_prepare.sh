#PBS -S /bin/bash
#PBS -N prepare_maxatac
#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=8:mem=64gb
#PBS -o ~/HRJ_monocytes/ATAC-seq/maxatac/OU
#PBS -e ~/HRJ_monocytes/ATAC-seq/maxatac/ER
#PBS -J 1-34

DIR=~/HRJ_monocytes/ATAC-seq/maxatac
cd ${DIR}/bams_withchr
mysample=$(head -n $PBS_ARRAY_INDEX samples.list | tail -n 1)

source activate maxATAC
cd ${DIR}

### Prepare the maxatac bigwig from each bam file (using the clean bam file per biorep, with technical reps merged where appropriate). About the "-skip_dedup" flag: "Include this flag to perform PCR deduplication of the input BAM file if you know that it has not been deduplicated." but DO include it, because otherwise it tries to deduplicate! 
maxatac prepare -i ${DIR}/bams_withchr/${mysample}*.bam -skip_dedup -o ${DIR}/output/${mysample} -prefix ${mysample}


conda deactivate
