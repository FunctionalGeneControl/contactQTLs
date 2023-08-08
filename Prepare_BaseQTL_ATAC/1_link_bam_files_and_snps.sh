#PBS -S /bin/bash
#PBS -N mod_bams
#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=2:mem=12gb
#PBS -o ~/HRJ_monocytes/AS_ATAC/BaseQTL/OU
#PBS -e ~/HRJ_monocytes/AS_ATAC/BaseQTL/ER

### We have decided to run the refbias analysis on the de-duplicated bam files.

BAMDIR=~/HRJ_monocytes/AS_ATAC/BaseQTL/bams_dedup
mkdir -p $BAMDIR
cd $BAMDIR
cp ~/HRJ_monocytes/ATAC-seq/clean/non_merge.list ./
cp ~/HRJ_monocytes/ATAC-seq/clean/merge.list ./
cat non_merge.list merge.list | sort -u > samples.list

for mysample in $(cat non_merge.list); do
	ln -s ~/HRJ_monocytes/ATAC-seq/clean/${mysample}/${mysample}*.filtered.unique.bam ./${mysample}.bam; 
done

for mysample in $(cat merge.list); do
        ln -s ~/HRJ_monocytes/ATAC-seq/clean_merged/${mysample}*.merged.filtered.unique.bam ./${mysample}.bam;
done

cd ~/HRJ_monocytes/AS_ATAC/BaseQTL
mkdir -p allele_counts
mkdir -p AI
mkdir -p ER
mkdir -p phASER_in
mkdir -p phASER_out
mkdir -p refbias

cd phASER_in
cp ~/HRJ_monocytes/AS_CHiC/BaseQTL/phASER_in/* ./

cd ~/HRJ_monocytes/AS_ATAC/BaseQTL/refbias
cp -r ~/HRJ_monocytes/AS_CHiC/BaseQTL/refbias/all_snps ./

mkdir -p step1
mkdir -p step2
mkdir -p step3
mkdir -p AI







