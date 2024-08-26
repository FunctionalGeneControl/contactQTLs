#PBS -S /bin/bash
#PBS -N TOBIAS_34r
#PBS -l walltime=48:00:00
#PBS -l select=1:ncpus=20:mem=240gb
#PBS -o ~/HRJ_monocytes/ATAC-seq/TOBIAS_snakemake/OU
#PBS -e ~/HRJ_monocytes/ATAC-seq/TOBIAS_snakemake/ER

module load anaconda3/personal
source activate TOBIAS_ENV

# Run the snakemake, which processes the data from the original bam files
cd ~/HRJ_monocytes/ATAC-seq/TOBIAS_snakemake
snakemake --configfile config/mono_34reps_cleaned_config.yaml --cores 20 --use-conda --keep-going

conda deactivate
