#PBS -S /bin/bash
#PBS -N get_pileups
#PBS -l walltime=1:00:00
#PBS -l select=1:ncpus=2:mem=32gb
#PBS -o ~/HRJ_monocytes/AS_CHiC/WASP/mono_34reps_allReads_using_eGene_filters/OU
#PBS -e ~/HRJ_monocytes/AS_CHiC/WASP/mono_34reps_allReads_using_eGene_filters/ER

# Get allele pileups using R script

source activate DT_DPLYR
cd ~/HRJ_monocytes/AS_CHiC/scripts/mono_34reps_analysis_using_eGene_filters

Rscript 6_collate_allele_pileups_paper.R

conda deactivate
