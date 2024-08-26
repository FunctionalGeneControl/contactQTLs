#PBS -S /bin/bash
#PBS -N featureCounts
#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=32:mem=64gb
#PBS -o /rds/general/project/lms-spivakov-analysis/live/HRJ_monocytes/eqtls/WASP/mono_34reps/OU
#PBS -e /rds/general/project/lms-spivakov-analysis/live/HRJ_monocytes/eqtls/WASP/mono_34reps/ER

source activate DT_DPLYR

cd /rds/general/project/lms-spivakov-analysis/live/HRJ_monocytes/eqtls/WASP/scripts
Rscript 4_get_total_reads.R

conda deactivate

