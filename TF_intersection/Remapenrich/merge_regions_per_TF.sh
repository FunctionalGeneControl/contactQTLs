#PBS -S /bin/bash
#PBS -N merge_TF
#PBS -l walltime=2:00:00
#PBS -l select=1:ncpus=8:mem=64gb
#PBS -o ~/HRJ_monocytes/remapenrich/OU
#PBS -e ~/HRJ_monocytes/remapenrich/ER

##############################################################################
### Merge regions per transcription factor in ATAC predictions or ChIP
##############################################################################

cd ~/HRJ_monocytes/remapenrich/Rev/input_by_TF
source activate samtools_env

ls *.bed > allfiles.txt

for file in $(cat allfiles.txt)
do
sort -k1,1 -k2,2n $file > ${file}.sorted.bed
bedtools merge -i ${file}.sorted.bed > ${file}.sorted.bed.merged.bed
done

## remove all files size zero
find . -size 0 -print -delete

### check then the script remapenrich_GUESS_BASE_monoFourSets_exploration.R for how these were used

##############################################################################
### Now, want to run remapenrich on just the ATAC predicted sites.
##############################################################################

### take the file mono_remap_input_fourSets_into2categories_merged.txt and split by ATAC or ChIP.
cd ~/HRJ_monocytes/remapenrich/Rev/
grep ATAC mono_remap_input_fourSets_into2categories_merged.txt > mono_remap_input_fourSets_into2categories_merged_ATAC_category.txt

##############################################################################
