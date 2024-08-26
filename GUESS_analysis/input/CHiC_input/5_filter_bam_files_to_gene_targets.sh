#PBS -S /bin/bash
#PBS -N process_bams
#PBS -l walltime=12:00:00
#PBS -l select=1:ncpus=4:mem=12gb
#PBS -o ~/HRJ_monocytes/AS_CHiC/WASP/mono_34reps_allReads_using_eGene_filters/OU
#PBS -e ~/HRJ_monocytes/AS_CHiC/WASP/mono_34reps_allReads_using_eGene_filters/ER
#PBS -J 1-34

###### How to get AI per eQTL/eGene
###### Following the same strategy as used for BaseQTL input for CHi-C

source activate refBias
DIR=~/HRJ_monocytes/AS_CHiC/WASP/mono_34reps_allReads_using_eGene_filters
bams=${DIR}/unique

cd ${DIR}/bams
INFILE=$(head -n $PBS_ARRAY_INDEX samples.list | tail -n 1)

### 1. Convert bam files (WASP filtered) to bedpe

## Make bedpe from the WASP filtered bam files
cd ${DIR}/unique
bedtools bamtobed -i ${INFILE}.keep.merge.unique.bam -bedpe > ${DIR}/bedpe/${INFILE}.bedpe


### 2. Intersect the bedpe with the one detailing interactions between eQTLs and eGenes (EDIT THIS BEDPE TO SAY GENE 1, GENE 2 etc)
cd ${DIR}/bedpe
eQTL_to_egene=~/HRJ_monocytes/AS_CHiC/scripts/mono_34reps_analysis_using_eGene_filters/bedpe_for_chic_egenes_with_Gene_number.txt
## Run pairtopair, requiring 5% overlap in each case. Ignore strand.
bedtools pairtopair -a ${INFILE}.bedpe -b ${eQTL_to_egene} -f 0.05 -is > ${INFILE}_eQTLs_to_eGenes.bedpe

### 3. Split the output bedbe by gene number
cd ${DIR}/bedpe
# In awk you have to set the variable for the name
awk -v myname=$INFILE 'OFS="\t" {print $0 > myname "_eQTLs_to_eGenes." $18 ".bedpe"}' ${INFILE}_eQTLs_to_eGenes.bedpe 

### 4. Make the resultant bam files (Up to 4 per sample)
#mkdir -p ${DIR}/bams_per_gene

for Gene in {1..4}; do
	echo $INFILE
	echo $Gene
	cd ${DIR}/bedpe
	awk 'OFS="\t"{print $7}' ${INFILE}_eQTLs_to_eGenes.Gene_${Gene}.bedpe | sort -u > ${DIR}/bams_per_gene/${INFILE}.Gene${Gene}.ids
	cd ${DIR}/bams_per_gene
	samtools view -H ${bams}/${INFILE}.keep.merge.unique.bam -o ${INFILE}.header
	samtools view ${bams}/${INFILE}.keep.merge.unique.bam | fgrep -w -f ${INFILE}.Gene${Gene}.ids | cat ${INFILE}.header - | samtools view -Sb - > ${INFILE}.Gene${Gene}.bam
	rm ${INFILE}.header;
done


### 5. Use phASER to count alleles in these files - see following script
