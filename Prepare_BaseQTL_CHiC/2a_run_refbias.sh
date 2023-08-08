#PBS -S /bin/bash
#PBS -N run_refbias
#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=32:mem=124gb
#PBS -o ~/HRJ_monocytes/AS_CHiC/BaseQTL/refbias/OU
#PBS -e ~/HRJ_monocytes/AS_CHiC/BaseQTL/refbias/ER
#PBS -J 1-34

#### Prepare the input files for BaseQTL
DIR=~/HRJ_monocytes/AS_CHiC/BaseQTL
cd ${DIR}/bams
INFILE=$(head -n $PBS_ARRAY_INDEX samples.list | tail -n 1) # these are the filt.bam files. if using the deduplicated bam files use bams_dedup 

source activate refBias # Trying again to make a refbias env.
## Uses python 3.9.7, pysam 0.18.0, numpy 1.22.3, pandas 1.4.1, pytables 3.7.0, scipy 1.8.0

cd ${DIR}
mysnps=./refbias/all_snps # in this folder we just need the snps.txt files. I am running against all SNPs MAF5.
mybam=./bams_dedup/${INFILE}.bam

### 1. Intersect SNPs and create fastq files for remapping
python ~/HRJ_monocytes/BaseQTL/baseqtl_pipeline/input/refbias/intersecting_snps_v2.py --is_paired_end --base_qual 10 --output_dir refbias/step1 --snp_dir ${mysnps} ${mybam}


# Don't worry that the output tells us that none of the reads have pairs. Checked against WASP and RefBias is correctly writing the fastq files.

### 2. Re-run HiCUP and sort by coordinate
PATH=$PATH:~/spivakov/bin/HiCUP_0.8.2

cd ${DIR}
module load bowtie2/2.1.0

hicup_mapper \
--bowtie2 /apps/bowtie2/2.1.0/bowtie2 \
--format Illumina_1.5 \
--index ~/spivakov/MDR/GRCh38/bowtie2/Homo_sapiens.GRCh38 \
--outdir ${DIR}/refbias/step2 \
--threads 8 \
--zip \
	./refbias/step1/${INFILE}.remap.fq1.gz ./refbias/step1/${INFILE}.remap.fq2.gz

samtools sort -o ./refbias/step2/${INFILE}.remap.fq1_2.pair.sort.bam ./refbias/step2/${INFILE}.remap.fq1_2.pair.bam
samtools index ./refbias/step2/${INFILE}.remap.fq1_2.pair.sort.bam


### 3. Run the next refbias script
python ~/HRJ_monocytes/BaseQTL/baseqtl_pipeline/input/refbias/intersecting_snps_post_remap_v2.py --is_paired_end --base_qual 10 --is_sorted --output_dir refbias/step3 --snp_dir ${mysnps} ${DIR}/refbias/step2/${INFILE}.remap.fq1_2.pair.sort.bam

### The AI script then needs to be run on all of the samples at once

conda deactivate
