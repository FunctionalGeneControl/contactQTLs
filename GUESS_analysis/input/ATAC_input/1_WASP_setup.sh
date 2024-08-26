#PBS -S /bin/bash
#PBS -N WASP
#PBS -l walltime=1:00:00
#PBS -l select=1:ncpus=1:mem=12gb
#PBS -o ~/HRJ_monocytes/AS_ATAC/WASP/OU
#PBS -e ~/HRJ_monocytes/AS_ATAC/WASP/ER


##########################################################

############ This file sets up the directories needed for WASP. It links in the bam files and renames them as appropriate. Here, I have renamed them according to the names that are in the genotype files. Probably use these names for ATAC-seq as well. 

TESTNAME=mono_34reps_noMultiMapped

####### These are the non dedup samlpes with the greatest number of reads. Should combine but not dedup to get the most reads. For ATAC-seq also need non-deduplicated BAMs. Will use sorted by read. NB these bam files do not include multimapped reads.

SAMPLE1=S025NM-01/S025NM-01_25K.filtered.dup.bam
RENAME1=S025NM-01
SAMPLE2=S025QG-02/S025QG-02_25K.filtered.dup.bam
RENAME2=S025QG-02
SAMPLE3=S025RE-03/S025RE-03_25K.filtered.dup.bam
RENAME3=S025RE-03
SAMPLE4=S025SC-04/S025SC-04_25K.filtered.dup.bam
RENAME4=S025SC-04
SAMPLE5=S025TA-05/S025TA-05_50K_1234.filtered.dup.bam
RENAME5=S025TA-05
SAMPLE6=S025U8-06/S025U8-06_25K.filtered.dup.bam
RENAME6=S025U8-06
SAMPLE7=S025V6-07/S025V6-07_25K.filtered.dup.bam
RENAME7=S025V6-07
SAMPLE8=S025W4-08.merged.filtered.dup.bam
RENAME8=S025W4-08
SAMPLE9=S025X2-09.merged.filtered.dup.bam
RENAME9=S025X2-09
SAMPLE10=S025Y0-10/S025Y0-10_100K_123.filtered.dup.bam
RENAME10=S025Y0-10
SAMPLE11=S0260R-11/S0260R-11_50K.filtered.dup.bam
RENAME11=S0260R-11
SAMPLE12=S0262N-12/S0262N-12_50K.filtered.dup.bam
RENAME12=S0262N-12
SAMPLE13=S0263L-13/S0263L-13_100K_1234.filtered.dup.bam
RENAME13=S0263L-13
SAMPLE14=S0264J-14/S0264J-14_50K_1234.filtered.dup.bam
RENAME14=S0264J-14
SAMPLE15=S0265H-15/S0265H-15_100K_1234.filtered.dup.bam
RENAME15=S0265H-15
SAMPLE16=S0267D-16/S0267D-16_50K_1234.filtered.dup.bam
RENAME16=S0267D-16
SAMPLE17=S0268B-17.merged.filtered.dup.bam
RENAME17=S0268B-17
SAMPLE18=S026A7-18/S026A7-18_100K_123.filtered.dup.bam
RENAME18=S026A7-18
SAMPLE19=S026B5-19/S026B5-19_50K_1234.filtered.dup.bam
RENAME19=S026B5-19
SAMPLE20=S026C3-20/S026C3-20_50K_1234.filtered.dup.bam
RENAME20=S026C3-20
SAMPLE21=S026D1-21/S026D1-21_50KP_1234.filtered.dup.bam
RENAME21=S026D1-21
SAMPLE22=S026FY-22/S026FY-22_50K_1234.filtered.dup.bam
RENAME22=S026FY-22
SAMPLE23=S026GW-23/S026GW-23_50K_1234.filtered.dup.bam
RENAME23=S026GW-23
SAMPLE24=S026HU-24/S026HU-24_50K_1234.filtered.dup.bam
RENAME24=S026HU-24
SAMPLE25=S026JQ-25/S026JQ-25_50K_1234.filtered.dup.bam
RENAME25=S026JQ-25
SAMPLE26=S026KO-26/S026KO-26_50K_1234.filtered.dup.bam
RENAME26=S026KO-26
SAMPLE27=S026MK-27/S026MK-27_50K_1234.filtered.dup.bam
RENAME27=S026MK-27
SAMPLE28=S026NI-28/S026NI-28_50K_1234.filtered.dup.bam
RENAME28=S026NI-28
SAMPLE29=S026S8-29/S026S8-29_50K_1234.filtered.dup.bam
RENAME29=S026S8-29
SAMPLE30=S026T6-30.merged.filtered.dup.bam
RENAME30=S026T6-30
SAMPLE31=S026XZ-31/S026XZ-31_100K_1234.filtered.dup.bam
RENAME31=S026XZ-31
SAMPLE32=S026YX-32/S026YX-32_50K_1234.filtered.dup.bam
RENAME32=S026YX-32
SAMPLE33=S0270N-33/S0270N-33_100K_1234.filtered.dup.bam
RENAME33=S0270N-33
SAMPLE34=S0271L-34/S0271L-34_50K_123.filtered.dup.bam
RENAME34=S0271L-34

BAMS="~/HRJ_monocytes/ATAC-seq/clean"
MERGED="~/HRJ_monocytes/ATAC-seq/clean_merged"
DIR="~/HRJ_monocytes/AS_ATAC/WASP"
SNPS="~/HRJ_monocytes/genotyping/imputed/hg38_filtered"

## Your installation of WASP.
PATH=$PATH:~/bin/WASP

### Make a new directory
cd ${DIR}
mkdir ${TESTNAME}
cd ${TESTNAME}
mkdir bams
cd bams

### Put bam files into directory with correct names and make a list.
ln -s ${BAMS}/${SAMPLE1} ./${RENAME1}.bam
ln -s ${BAMS}/${SAMPLE2} ./${RENAME2}.bam
ln -s ${BAMS}/${SAMPLE3} ./${RENAME3}.bam
ln -s ${BAMS}/${SAMPLE4} ./${RENAME4}.bam
ln -s ${BAMS}/${SAMPLE5} ./${RENAME5}.bam
ln -s ${BAMS}/${SAMPLE6} ./${RENAME6}.bam
ln -s ${BAMS}/${SAMPLE7} ./${RENAME7}.bam
ln -s ${MERGED}/${SAMPLE8} ./${RENAME8}.bam
ln -s ${MERGED}/${SAMPLE9} ./${RENAME9}.bam
ln -s ${BAMS}/${SAMPLE10} ./${RENAME10}.bam
ln -s ${BAMS}/${SAMPLE11} ./${RENAME11}.bam
ln -s ${BAMS}/${SAMPLE12} ./${RENAME12}.bam
ln -s ${BAMS}/${SAMPLE13} ./${RENAME13}.bam
ln -s ${BAMS}/${SAMPLE14} ./${RENAME14}.bam
ln -s ${BAMS}/${SAMPLE15} ./${RENAME15}.bam
ln -s ${BAMS}/${SAMPLE16} ./${RENAME16}.bam
ln -s ${MERGED}/${SAMPLE17} ./${RENAME17}.bam
ln -s ${BAMS}/${SAMPLE18} ./${RENAME18}.bam
ln -s ${BAMS}/${SAMPLE19} ./${RENAME19}.bam
ln -s ${BAMS}/${SAMPLE20} ./${RENAME20}.bam
ln -s ${BAMS}/${SAMPLE21} ./${RENAME21}.bam
ln -s ${BAMS}/${SAMPLE22} ./${RENAME22}.bam
ln -s ${BAMS}/${SAMPLE23} ./${RENAME23}.bam
ln -s ${BAMS}/${SAMPLE24} ./${RENAME24}.bam
ln -s ${BAMS}/${SAMPLE25} ./${RENAME25}.bam
ln -s ${BAMS}/${SAMPLE26} ./${RENAME26}.bam
ln -s ${BAMS}/${SAMPLE27} ./${RENAME27}.bam
ln -s ${BAMS}/${SAMPLE28} ./${RENAME28}.bam
ln -s ${BAMS}/${SAMPLE29} ./${RENAME29}.bam
ln -s ${MERGED}/${SAMPLE30} ./${RENAME30}.bam
ln -s ${BAMS}/${SAMPLE31} ./${RENAME31}.bam
ln -s ${BAMS}/${SAMPLE32} ./${RENAME32}.bam
ln -s ${BAMS}/${SAMPLE33} ./${RENAME33}.bam
ln -s ${BAMS}/${SAMPLE34} ./${RENAME34}.bam

ls -1 *.bam > bams.list
cp ~/HRJ_monocytes/AS_CHiC/WASP/mono_34reps_allReads/bams/samples.list ./

# OR could've done this:
#cat bams.list | sed 's/.bam//g' > samples.list

### Make other directories. These will be the same for each WASP run.
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

