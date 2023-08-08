#PBS -S /bin/bash
#PBS -N get_bedpe
#PBS -l walltime=1:00:00
#PBS -l select=1:ncpus=2:mem=32gb
#PBS -o ~/HRJ_monocytes/AS_CHiC/BaseQTL/OU
#PBS -e ~/HRJ_monocytes/AS_CHiC/BaseQTL/ER


source activate DT_DPLYR
cd ~/HRJ_monocytes/AS_CHiC/BaseQTL/scripts

Rscript 3_get_eQTL_eGene_bedpe.R

conda deactivate
