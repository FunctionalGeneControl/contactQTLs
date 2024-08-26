#PBS -S /bin/bash
#PBS -N WASP
#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=32:mem=124gb
#PBS -o ~/HRJ_monocytes/eqtls/WASP/mono_34reps/OU
#PBS -e ~/HRJ_monocytes/eqtls/WASP/mono_34reps/ER
#PBS -J 1-34

##### This script finds SNPs intersecting with the bam reads. It flips the alleles and tries to remap the reads. If reads map elsewhere, they are discarded.

module load anaconda3/personal
source activate WASP

##########################################################

TESTNAME=mono_34reps

DIR="~/HRJ_monocytes/eqtls/WASP"
CHIC="~/HRJ_monocytes/AS_CHiC/WASP/input"
SNPS="~/HRJ_monocytes/genotyping/imputed/hg38_filtered"
PATH=$PATH:~/bin/WASP

## Step 1 is to make the hd5 SNP files - we can use the ones already made in the CHiC folder
## Make HDF5 SNPs files - this is done already for all SNPs at MAF >= 0.05
## Note that sample names were also changed to the final version (e.g. S025NM-01)

## Step 2 is to map the fastq files, this I have already done with STAR in the BaseQTL pipeline.

## Step 3: Find intersecting SNPs - find reads that have mapping biases

cd ${DIR}/${TESTNAME}/bams
INFILE=$(head -n $PBS_ARRAY_INDEX samples.list | tail -n 1)
#
## Note that even though RNA-seq, it is paired end.
python ~/bin/WASP/mapping/find_intersecting_snps.py \
	--is_paired_end \
        --output_dir ${DIR}/${TESTNAME}/find_intersecting_snps \
        --snp_tab ${CHIC}/snp_tab_maf05_finalRenamed.h5 \
        --snp_index ${CHIC}/snp_index_maf05_finalRenamed.h5 \
        --haplotype ${CHIC}/haplotypes_maf05_finalRenamed.h5 \
        ${INFILE}.bam

conda deactivate

### Step 4: Now need to rerun mapping
source activate refBias

cd ${DIR}/${TESTNAME}/map2
READ1=../find_intersecting_snps/${INFILE}.remap.fq1.gz
READ2=../find_intersecting_snps/${INFILE}.remap.fq2.gz

~/anaconda3/envs/refBias/bin/STAR \
	--runThreadN 16 \
	--genomeDir ~/HRJ_monocytes/eqtls/BaseQTL/STAR_idx/indices \
        --readFilesIn $READ1 $READ2 \
        --readFilesCommand zcat \
        --outSAMtype BAM SortedByCoordinate \
        --outFileNamePrefix ${INFILE}.remap

samtools index -b ${INFILE}.remap.bam ${INFILE}.remap.idx

conda deactivate

## Use filter_remapped_reads.py to filter out reads where one or more of the allelic versions of the reads fail to map back to the same location as the original read
source activate WASP
PATH=$PATH:~/bin/WASP
cd ${DIR}/${TESTNAME}
python ~/bin/WASP/mapping/filter_remapped_reads.py \
	./find_intersecting_snps/${INFILE}.to.remap.bam \
	./map2/${INFILE}.remapAligned.sortedByCoord.out.bam \
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










