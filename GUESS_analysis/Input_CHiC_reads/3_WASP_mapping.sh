#PBS -S /bin/bash
#PBS -N WASP
#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=32:mem=124gb
#PBS -o ~/HRJ_monocytes/AS_CHiC/WASP/mono_34reps_allReads_using_eGene_filters/OU
#PBS -e ~/HRJ_monocytes/AS_CHiC/WASP/mono_34reps_allReads_using_eGene_filters/ER
#PBS -J 1-34

##### This script finds SNPs intersecting with the bam reads. It flips the alleles and tries to remap the reads. If reads map elsewhere, they are discarded.

module load anaconda3/personal
#source activate WASP

##########################################################

#### NOT RUN: Copied over files from mono_34reps_allReads

##########################################################

TESTNAME=mono_34reps_allReads_using_eGene_filters

DIR="~/HRJ_monocytes/AS_CHiC/WASP"
SNPS="~/HRJ_monocytes/genotyping/imputed/hg38_filtered"
PATH=$PATH:~/bin/WASP

## Step 1 is to make the hd5 SNP files
## Make HDF5 SNPs files - this is done already for all SNPs at MAF >= 0.05
## Note that sample names were also changed to the final version (e.g. S025NM-01)

## Step 2 is to map the fastq files, this I have already done with hicup, not removing interim files. 

## Step 3: Find intersecting SNPs - find reads that have mapping biases
module load samtools/1.3.1

cd ${DIR}/${TESTNAME}/bams
INFILE=$(head -n $PBS_ARRAY_INDEX redo.list | tail -n 1)


cd ${INFILE}

python ~/bin/WASP/mapping/find_intersecting_snps.py \
        --is_paired_end \
        --output_dir ${DIR}/${TESTNAME}/find_intersecting_snps \
        --snp_tab ${DIR}/input/snp_tab_maf05_finalRenamed.h5 \
        --snp_index ${DIR}/input/snp_index_maf05_finalRenamed.h5 \
        --haplotype ${DIR}/input/haplotypes_maf05_finalRenamed.h5 \
        ${INFILE}.bam

conda deactivate

## Step 4: Now need to rerun hicup - just mapping, so we do not need to use HiCUP combinations.
## Check which version of HiCUP we are using. Now using the latest version as of October 2021.

PATH=$PATH:/rds/general/project/lms-spivakov-analysis/live/bin/HiCUP_0.8.2

cd ${DIR}/${TESTNAME}
source activate Renv1
module load bowtie2/2.1.0

hicup_mapper \
--bowtie2 /apps/bowtie2/2.1.0/bowtie2 \
--format Illumina_1.5 \
--index /rds/general/project/lms-spivakov-analysis/live/MDR/GRCh38/bowtie2/Homo_sapiens.GRCh38 \
--outdir ${DIR}/${TESTNAME}/map2 \
--threads 8 \
--zip \
./find_intersecting_snps/${INFILE}.remap.fq1.gz ./find_intersecting_snps/${INFILE}.remap.fq2.gz

conda deactivate

# Sort the bam file

## Hopefully HiCUP zipped it, so don't need to do this: samtools view -S -b ./map2/${INFILE}.remap.fq1_2.pair.sam > ./map2/${INFILE}.remap.fq1_2.pair.bam
#samtools sort -o ./map2/${INFILE}.remap.fq1_2.pair.sort.bam ./map2/${INFILE}.remap.fq1_2.pair.bam
#samtools index ./map2/${INFILE}.remap.fq1_2.pair.sort.bam

##Use filter_remapped_reads.py to filter out reads where one or more of the allelic versions of the reads fail to map back to the same location as the original read
### seeing if I should install full new version of WASP. NO, it is fine. Only the filter_remapped_reads.py script was updated (due to my CIGAR problem).
source activate WASP
PATH=$PATH:~/home/bin/WASP
cd ${DIR}/${TESTNAME}

python ~/bin/WASP/mapping/filter_remapped_reads.py \
./find_intersecting_snps/${INFILE}.to.remap.bam \
./map2/${INFILE}.remap.fq1_2.pair.sort.bam \
./filter_remapped_reads/${INFILE}.remap.keep.bam


# merge the reserved, non-re-mapped reads with the re-mapped reads to keep
samtools merge ./merge/${INFILE}.keep.merge.bam \
	       ./filter_remapped_reads/${INFILE}.remap.keep.bam  \
	       ./find_intersecting_snps/${INFILE}.keep.bam

samtools sort -o  ./merge/${INFILE}.keep.merge.sort.bam \
                  ./merge/${INFILE}.keep.merge.bam

samtools index ./merge/${INFILE}.keep.merge.sort.bam


## Filter duplicate reads. rmdup.py which performs unbiased removal of duplicate reads. The script discards duplicate reads at random (independent of their score).

python ~/bin/WASP/mapping/rmdup_pe.py \
./merge/${INFILE}.keep.merge.sort.bam \
./unique/${INFILE}.keep.merge.unique.bam

samtools sort -o  ./unique/${INFILE}.keep.merge.unique.sort.bam \
                  ./unique/${INFILE}.keep.merge.unique.bam

samtools index ./unique/${INFILE}.keep.merge.unique.sort.bam


conda deactivate
