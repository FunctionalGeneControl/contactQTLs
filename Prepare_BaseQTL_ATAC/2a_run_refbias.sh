#PBS -S /bin/bash
#PBS -N run_refbias
#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=32:mem=124gb
#PBS -o ~/HRJ_monocytes/AS_ATAC/BaseQTL/refbias/OU
#PBS -e ~/HRJ_monocytes/AS_ATAC/BaseQTL/refbias/ER
# Running again for S025RE-03

#### Prepare the input files for BaseQTL
DIR=~/HRJ_monocytes/AS_ATAC/BaseQTL
cd ${DIR}/bams_dedup
#INFILE=$(head -n $PBS_ARRAY_INDEX samples.list | tail -n 1) # these are the cleaned, deduplicated ATAC bam files.

source activate refBias # Trying again to make a refbias env.
## Uses python 3.9.7, pysam 0.18.0, numpy 1.22.3, pandas 1.4.1, pytables 3.7.0, scipy 1.8.0

cd ${DIR}
mysnps=./refbias/all_snps # in this folder we just need the snps.txt files. I am running against all SNPs MAF5.
mybam=./bams_dedup/${INFILE}.bam

### 1. Intersect SNPs and create fastq files for remapping
python ~/HRJ_monocytes/BaseQTL/baseqtl_pipeline/input/refbias/intersecting_snps_v2_pairs.py --is_paired_end --base_qual 10 --output_dir refbias/step1 --snp_dir ${mysnps} ${mybam} # this py script was repaired to call read pairs correctly


### 2. Re-run alignment and sort by coordinate
cd ${DIR}
module load bowtie2/2.2.9 # the same version as for initial alignment

bowtie2 --very-sensitive -p 4 -X 2000 \
-x /rds/general/project/lms-spivakov-analysis/live/MDR/GRCh38/bowtie2/Homo_sapiens.GRCh38 \
-1 ${DIR}/refbias/step1/${INFILE}.remap.fq1.gz \
-2 ${DIR}/refbias/step1/${INFILE}.remap.fq2.gz \
-S ${DIR}/refbias/step2/${INFILE}.remap.fq1_2.sam

# This bit was originally all piped, but it kept breaking so now I've separated each step
cd ${DIR}/refbias/step2
samtools view -h -b -o ${INFILE}.remap.fq1_2.bam ${INFILE}.remap.fq1_2.sam
samtools sort -o ${INFILE}.remap.fq1_2.sort.bam ${INFILE}.remap.fq1_2.bam
samtools index ${INFILE}.remap.fq1_2.sort.bam
rm ${INFILE}.remap.fq1_2.sam


### 3. Run the next refbias script
cd ${DIR}
INFILE=S025RE-03
python ~/HRJ_monocytes/BaseQTL/baseqtl_pipeline/input/refbias/intersecting_snps_post_remap_v2.py --is_paired_end --base_qual 10 --is_sorted --output_dir refbias/step3 --snp_dir ${mysnps} ${DIR}/refbias/step2/${INFILE}.remap.fq1_2.sort.bam

### 4. Run the AI script
cd ${DIR}
python ~/HRJ_monocytes/BaseQTL/baseqtl_pipeline/input/refbias/AI.py --initial-AI ./refbias/step1/${INFILE}.initial.AI.txt --post-remap_AI ./refbias/step3/${INFILE}.remap.fq1_2.sort.post_remapping_AI.txt --output_file ./refbias/AI/${INFILE}_finalAI.txt # check the name of post remap AI here 


conda deactivate
