#PBS -S /bin/bash
#PBS -N Probes
#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=2:mem=32gb
#PBS -o ~/ephemeral/
#PBS -e ~/ephemeral/

echo $PBS_O_WORKDIR
cd $PBS_O_WORKDIR


perl ~/eQTL_CHiC_design/design_all_grch38_primers_mboi.pl


