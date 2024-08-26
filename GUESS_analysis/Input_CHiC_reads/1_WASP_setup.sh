#PBS -S /bin/bash
#PBS -N WASP
#PBS -l walltime=1:00:00
#PBS -l select=1:ncpus=1:mem=12gb
#PBS -o ~/HRJ_monocytes/AS_CHiC/WASP/OU
#PBS -e ~/HRJ_monocytes/AS_CHiC/WASP/ER


##########################################################

############ This file sets up the directories needed for WASP. I have already linked the non-dedup, captured bam files into there.

TESTNAME=mono_34reps_allReads_using_eGene_filters

DIR="~/HRJ_monocytes/AS_CHiC/WASP"
SNPS="~/HRJ_monocytes/genotyping/imputed/hg38_filtered"

## Your installation of WASP.
PATH=$PATH:/rds/general/user/hrayjone/home/bin/WASP

### Make a list of samples (name and files)
cd ${DIR}/${TESTNAME}/bams
#ls -1 *.bam > bams.list
#cat bams.list | sed 's/.filt.captured.bam//g' > samples.list

### Make other directories. These will be the same for each WASP run.
mkdir ${DIR}/${TESTNAME}/find_intersecting_snps
mkdir ${DIR}/${TESTNAME}/map2
mkdir ${DIR}/${TESTNAME}/filter_remapped_reads
mkdir ${DIR}/${TESTNAME}/merge
mkdir ${DIR}/${TESTNAME}/unique
mkdir ${DIR}/${TESTNAME}/samples
mkdir ${DIR}/${TESTNAME}/bedpe
mkdir ${DIR}/${TESTNAME}/interactions
mkdir ${DIR}/${TESTNAME}/bams_filtered
mkdir ${DIR}/${TESTNAME}/ER
mkdir ${DIR}/${TESTNAME}/OU

