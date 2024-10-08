#PBS -S /bin/bash
#PBS -N run_phASER
#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=8:mem=64gb
#PBS -o ~/HRJ_monocytes/AS_CHiC/WASP/mono_34reps_allReads_using_eGene_filters/OU
#PBS -e ~/HRJ_monocytes/AS_CHiC/WASP/mono_34reps_allReads_using_eGene_filters/ER
#PBS -J 1-34

###########  CONDA ENVS ###########
###      phASER requires Py2    ###
###        Use "phASER"         ###
###  baseQTL input requires Py3 ###
###      	use "refBias"		###
###################################

source activate phASER

### Use phASER to count haplotypes at heterozygous SNPs in CHiC data. 
### Using de-duplicated files from hicup, that have been filtered to interactions with eGenes and split by the eGene in question. Refer back to the bed file to see which genes they were.
### HiCUP removes multimapping reads, so not providing a mapping quality.
### Will provide the default base quality of 10.

BASEDIR=~/HRJ_monocytes/AS_CHiC/WASP/mono_34reps_allReads_using_eGene_filters
INPUT=${BASEDIR}/bams_per_gene
OUTPUT=${BASEDIR}/phASER_out
mkdir -p ${OUTPUT}
VCF=~/HRJ_monocytes/genotyping/imputed/hg38_filtered/all.rsq3.maf05.hg38.dose.finalRenamed.newIDs.vcf.gz
PHASER=~/bin/phaser

cd ${BASEDIR}/bams
INFILE=$(head -n $PBS_ARRAY_INDEX samples.list | tail -n 1)

# Bam files need to be indexed
cd ${INPUT}
for Gene in {1..4}; do
	samtools sort ${INFILE}.Gene${Gene}.bam -o ${INFILE}.Gene${Gene}.sorted.bam
	samtools index ${INFILE}.Gene${Gene}.sorted.bam;
done

# Run phASER with only the HLA blacklist
cd ${BASEDIR}

for Gene in {1..4}; do 
python ${PHASER}/phaser/phaser.py \
	--vcf ${VCF} \
       	--bam ${INPUT}/${INFILE}.Gene${Gene}.sorted.bam \
       	--mapq 0 --paired_end 1 --baseq 10 --threads 4 --unphased_vars 1 \
	--blacklist ~/HRJ_monocytes/AS_CHiC/BaseQTL/phASER_in/hg38_hla.bed \
       	--sample ${INFILE} \
	--o ${OUTPUT}/${INFILE}_Gene${Gene}_noBlacklist;
done


conda deactivate

