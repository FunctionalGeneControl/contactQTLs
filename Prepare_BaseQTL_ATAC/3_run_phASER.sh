#PBS -S /bin/bash
#PBS -N run_phASER
#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=8:mem=96gb
#PBS -o ~/HRJ_monocytes/AS_ATAC/BaseQTL/OU
#PBS -e ~/HRJ_monocytes/AS_ATAC/BaseQTL/ER
#PBS -J 1-34

###########  CONDA ENVS ###########
###      phASER requires Py2    ###
###        Use "phASER"         ###
###  baseQTL input requires Py3 ###
###       Tried "baseqtl"       ###
###   Made by IT, do not edit!! ###
###	But, doesnt work	###
###	So now using "refBias"	###
###################################

source activate phASER

### Use phASER to count haplotypes at heterozygous SNPs in ATAC data. 

BASEDIR=~/HRJ_monocytes/AS_ATAC/BaseQTL
INPUT=${BASEDIR}/bams_dedup
OUTPUT=${BASEDIR}/phASER_out
VCF=~/HRJ_monocytes/genotyping/imputed/hg38_filtered/all.rsq3.maf05.hg38.dose.finalRenamed.newIDs.vcf.gz
PHASER=~/bin/phaser

cd ${BASEDIR}/bams_dedup
INFILE=$(head -n $PBS_ARRAY_INDEX samples.list | tail -n 1)

# Bam files need to be indexed
cd ${INPUT}
samtools sort ${INFILE}.bam -o ${INFILE}.sorted.bam
samtools index ${INFILE}.sorted.bam

# Run phASER with only the HLA blacklist
cd ${BASEDIR}

python ${PHASER}/phaser/phaser.py \
	--vcf ${VCF} \
       	--bam ${INPUT}/${INFILE}.sorted.bam \
       	--mapq 0 --paired_end 1 --baseq 10 --threads 4 --unphased_vars 1 \
	--blacklist ${BASEDIR}/phASER_in/hg38_hla.bed \
       	--sample ${INFILE} \
	--o ${OUTPUT}/${INFILE}_noBlacklist

# Run phASER with haplo blacklist
python ${PHASER}/phaser/phaser.py \
        --vcf ${VCF} \
        --bam ${INPUT}/${INFILE}.sorted.bam \
        --mapq 0 --paired_end 1 --baseq 10 --threads 4 --unphased_vars 1 \
        --blacklist ${BASEDIR}/phASER_in/hg38_hla.bed \
	--haplo_count_blacklist ${BASEDIR}/phASER_in/hg38_haplo_count_blacklist.bed \
        --sample ${INFILE} \
        --o ${OUTPUT}/${INFILE}_withBlacklist


conda deactivate

