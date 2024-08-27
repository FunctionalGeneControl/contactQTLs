#!/bin/bash
#$ -S /bin/sh
#$ -N splitSNP
# execute the job from the current working directory:
#$ -cwd
#$ -l cpu=8
# Specify locations for STOUT and STERROR:
#$ -e /data/cmn_vamal/ChIP/ER
#$ -o /data/cmn_vamal/ChIP/OU

# Tell it where to find conda: (note "source" is substituted with a dot)
. ~/miniconda3/etc/profile.d/conda.sh

source activate py2

PATH=$PATH:/home/AD/hrayjones/bin/splitSNP

cd /data/cmn_vamal/ChIP/chip-ap/mono_CTCF/08_results

python /home/AD/hrayjones/bin/splitSNP/splitSNP.py mono_CTCF_chip_rep1.bam \
       mono_CTCF_chip_rep1_rs7146599 chr14:24058363:G:A

python /home/AD/hrayjones/bin/splitSNP/splitSNP.py mono_CTCF_ctrl_rep1.bam \
       mono_CTCF_ctrl_rep1_rs7146599 chr14:24058363:G:A

# the G allele should have more CTCF
# The A allele should have less CTCF


conda deactivate
