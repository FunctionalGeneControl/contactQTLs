#PBS -S /bin/bash
#PBS -N bedtools
#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=8:mem=64gb
#PBS -o ~/HRJ_monocytes/AS_ATAC/BaseQTL/snps
#PBS -e ~/HRJ_monocytes/AS_ATAC/BaseQTL/snps

##########################################################

### This script was run for the WASP pipeline. Can just use the output here as well.
### I have made a new snp file that should have the correct hg38 strand for the alleles and IDs

# To see how the input files were made, see the script assign_snps_peaks.R

DIR="~/HRJ_monocytes/AS_ATAC/BaseQTL"
SNPS=${DIR}/snps/all_snp_pos_rsq3_maf05.bed
PEAKS=${DIR}/snps/HMMRATAC_peaks_mono34reps_withBlacklist.bed

##########################################################

source activate samtools_env

cd ${DIR}/snps

### Assign snps to the closest ATAC-seq peak

#a is snps
#b is peaks
#-d report the distance
#-t all incase of ties, report all

# Sort the files!
# Sorted the HMMRATAC.bed file on the command line, and removed X and Y.
# Sorted the snps.bed file also on the command line.

bedtools closest -d -t all -a ${SNPS} -b ${PEAKS} > closest_peaks.bed

# Require that a snp is within 5kb of a peak in order for it to be associated
awk 'OFS="\t"{if ($9 <= 5000) print $0}' closest_peaks.bed > closest_peaks_5kb.bed 

# This file can be used when collating allele counts.

