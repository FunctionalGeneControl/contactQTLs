#PBS -S /bin/bash
#PBS -N process_bams
#PBS -l walltime=12:00:00
#PBS -l select=1:ncpus=32:mem=24gb
#PBS -o ~/HRJ_monocytes/AS_CHiC/BaseQTL/OU
#PBS -e ~/HRJ_monocytes/AS_CHiC/BaseQTL/ER
#PBS -J 1-34

###### How to get AI per eQTL/eGene

source activate refBias
DIR=~/HRJ_monocytes/AS_CHiC/BaseQTL
bams=${DIR}/bams_dedup

cd $bams
INFILE=$(head -n $PBS_ARRAY_INDEX samples.list | tail -n 1)

### 1. Convert bam files (the ones used as input to refbias) to bedpe

## Get rid of X and Y and MT and make bedpe
cd $bams
bedtools bamtobed -i ${INFILE}.bam -bedpe > ${DIR}/bedpe/${INFILE}.bedpe.temp
cd ${DIR}/bedpe
awk 'OFS="\t"{if (($1 != "X") && ($1 != "Y") && ($1 != "MT") && ($4 != "X") && ($4 != "Y") && ($4 != "MT")) print $0}' ${INFILE}.bedpe.temp > ${INFILE}.bedpe

rm ${INFILE}.bedpe.temp

### 2. Intersect the bedpe with the one detailing interactions between eQTLs and eGenes (EDIT THIS BEDPE TO SAY GENE 1, GENE 2 etc)
cd ${DIR}/bedpe
eQTL_to_egene=${DIR}/scripts/bedpe_for_chic_egenes_with_Gene_number_CORRECTED.txt # now unique numbers per gene/DpnII/SNP combo
## Run pairtopair, requiring 5% overlap in each case. Ignore strand.
bedtools pairtopair -a ${INFILE}.bedpe -b ${eQTL_to_egene} -f 0.05 -is > ${INFILE}_eQTLs_to_eGenes.bedpe

### 3. Split the output bedbe by gene number
cd ${DIR}/bedpe
# In awk you have to set the variable for the name
awk -v myname=$INFILE 'OFS="\t" {print $0 > myname "_eQTLs_to_eGenes." $18 ".bedpe"}' ${INFILE}_eQTLs_to_eGenes.bedpe 

### 4. Make the resultant bam files (Up to 33 per sample - orignally was just 4)
for Gene in {1..33}; do
	echo $INFILE
	echo $Gene
	cd ${DIR}/bedpe
	awk 'OFS="\t"{print $7}' ${INFILE}_eQTLs_to_eGenes.Gene_${Gene}.bedpe | sort -u > ${DIR}/bams_per_gene/${INFILE}.Gene${Gene}.ids
	cd ${DIR}/bams_per_gene
	samtools view -H ${bams}/${INFILE}.bam -o ${INFILE}.header
	samtools view ${bams}/${INFILE}.bam | fgrep -w -f ${INFILE}.Gene${Gene}.ids | cat ${INFILE}.header - | samtools view -Sb - > ${INFILE}.Gene${Gene}.bam
	rm ${INFILE}.header;
done


### 5. Use phASER to count alleles in these files - see following script
