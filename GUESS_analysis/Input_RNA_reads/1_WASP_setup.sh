#PBS -S /bin/bash
#PBS -N WASP
#PBS -l walltime=1:00:00
#PBS -l select=1:ncpus=1:mem=12gb
#PBS -o ~/HRJ_monocytes/eqtls/WASP/OU
#PBS -e ~/HRJ_monocytes/eqtls/WASP/ER

##########################################################

############ This file sets up the directories needed for WASP. It links in the bam files and renames them as appropriate. Here, I have renamed them according to the names that are in the genotype files.

TESTNAME=mono_34reps

####### Get the bam files. These were already produced in the BaseQTL run.

BAMS="~/HRJ_monocytes/eqtls/BaseQTL/STAR"
DIR="~/HRJ_monocytes/eqtls/WASP"
SNPS="~/HRJ_monocytes/genotyping/imputed/hg38_filtered"

## Your installation of WASP.
PATH=$PATH:~/bin/WASP

### Make a new directory
cd ${DIR}
mkdir ${TESTNAME}
cd ${TESTNAME}
mkdir bams
cd bams

####### Bam files
mysamples=~/HRJ_monocytes/AS_CHiC/WASP/mono_34reps_allReads/bams/samples.list

for INDIVIDUAL in $(cat $mysamples)
do
    echo $INDIVIDUAL
    ln -s ${BAMS}/${INDIVIDUAL}/Aligned.sortedByCoord.out.bam ./${INDIVIDUAL}.bam
done

cp ${mysamples} ./samples.list
ls -1 *.bam > bams.list

mkdir ${DIR}/${TESTNAME}/find_intersecting_snps
mkdir ${DIR}/${TESTNAME}/map2
mkdir ${DIR}/${TESTNAME}/filter_remapped_reads
mkdir ${DIR}/${TESTNAME}/merge
mkdir ${DIR}/${TESTNAME}/unique
mkdir ${DIR}/${TESTNAME}/unique_targeted
mkdir ${DIR}/${TESTNAME}/samples
mkdir ${DIR}/${TESTNAME}/interactions
mkdir ${DIR}/${TESTNAME}/biostar214299
mkdir ${DIR}/${TESTNAME}/beds
mkdir ${DIR}/${TESTNAME}/allele_counts
mkdir ${DIR}/${TESTNAME}/for_elena
mkdir ${DIR}/${TESTNAME}/ER
mkdir ${DIR}/${TESTNAME}/OU


