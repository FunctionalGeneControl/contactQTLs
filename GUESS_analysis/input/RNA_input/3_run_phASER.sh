#PBS -S /bin/bash
#PBS -N run_phASER
#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=8:mem=96gb
#PBS -o ~/HRJ_monocytes/eqtls/WASP/mono_34reps/OU
#PBS -e ~/HRJ_monocytes/eqtls/WASP/mono_34reps/ER
#PBS -J 1-34

###########  CONDA ENVS ###########
###      phASER requires Py2    ###
###        Use "phASER"         ###
###  baseQTL input requires Py3 ###
###		   use "refBias"		###
###################################


### This script was originally used in the BaseQTL pipeline. Now using on WASP output bam files to check reference bias.

source activate phASER
TESTNAME=mono_34reps

### Use phASER to count haplotypes at heterozygous SNPs in RNA-seq data. 
### Note that this is not restricted to peaks perse.
### Using the de-duplicated output bam files from WASP.
### Try with and without hg38 blacklist region; see how many snps we lose.
### Will provide the default base quality of 10.

BASEDIR=~/HRJ_monocytes/eqtls/WASP
INPUT=${BASEDIR}/${TESTNAME}/unique
OUTPUT=${BASEDIR}/${TESTNAME}/phASER_out
VCF=~/HRJ_monocytes/genotyping/imputed/hg38_filtered/all.rsq3.maf05.hg38.dose.finalRenamed.newIDs.vcf.gz
PHASER=~/bin/phaser
BLACKLIST=~/HRJ_monocytes/AS_ATAC/BaseQTL/phASER_in

mkdir -p ${OUTPUT}

cd ${INPUT}
#ls *.keep.merge.unique.sort.bam > bams.list
#cat bams.list | sed 's/.keep.merge.unique.sort.bam//g' > samples.list
INFILE=$(head -n $PBS_ARRAY_INDEX samples.list | tail -n 1)

# Bam files need to be sorted and indexed, this is already done at the end of WASP.

# Run phASER with only the HLA blacklist
cd ${BASEDIR}/${TESTNAME}

# Run phASER with just the HLA blacklist
python ${PHASER}/phaser/phaser.py \
	--vcf ${VCF} \
       	--bam ${INPUT}/${INFILE}.keep.merge.unique.sort.bam \
       	--mapq 0 --paired_end 1 --baseq 10 --threads 4 --unphased_vars 1 \
	--blacklist ${BLACKLIST}/hg38_hla.bed \
       	--sample ${INFILE} \
	--o ${OUTPUT}/${INFILE}_noBlacklist

# Run phASER with haplo blacklist
python ${PHASER}/phaser/phaser.py \
        --vcf ${VCF} \
        --bam ${INPUT}/${INFILE}.keep.merge.unique.sort.bam \
        --mapq 0 --paired_end 1 --baseq 10 --threads 4 --unphased_vars 1 \
        --blacklist ${BLACKLIST}/hg38_hla.bed \
	--haplo_count_blacklist ${BLACKLIST}/hg38_haplo_count_blacklist.bed \
        --sample ${INFILE} \
        --o ${OUTPUT}/${INFILE}_withBlacklist


conda deactivate

