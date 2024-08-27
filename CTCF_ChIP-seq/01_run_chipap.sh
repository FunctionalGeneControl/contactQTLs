#!/bin/bash
#$ -S /bin/bash
#$ -N chip-ap
# execute the job from the current working directory:
#$ -cwd
#$ -l cpu=8
# Specify locations for STOUT and STERROR:
#$ -e /data/cmn_vamal/ChIP/ER
#$ -o /data/cmn_vamal/ChIP/OU

###################################################

# If using /bin/sh above, add the following:
# Tell it where to find conda: (note "source" is substituted with a dot)
#. ~/miniconda3/etc/profile.d/conda.sh

# If using /bin/bash above, do the following:
source ~/.bashrc
source ~/miniconda3/etc/profile.d/conda.sh # because the pipeline otherwise cannot activate the env


####################################################
# For chip-ap README, see:
# https://github.com/JSuryatenggara/ChIP-AP/tree/main

# Note that chip-ap includes adapter trimming and quality trimming of raw data. For now trying with the raw data - make sure that the adapters were adequately removed.

PATH=$PATH:~/bin/chip-ap/chipap_installation/chipap_scripts

conda activate chip-ap

cd /data/cmn_vamal/ChIP
RAW=/data/cmn_vamal/ChIP/monocytes_raw ### note: raw, not trimmed!
python ~/bin/chip-ap/chipap_installation/chipap_scripts/chipap.py --mode paired \
	--chipR1 ${RAW}/M18_CTCF_1.fq.gz ${RAW}/M28_CTCF_1.fq.gz ${RAW}/M33_CTCF_1.fq.gz \
	--chipR2 ${RAW}/M18_CTCF_2.fq.gz ${RAW}/M28_CTCF_2.fq.gz ${RAW}/M33_CTCF_2.fq.gz \
	--ctrlR1 ${RAW}/M18_input_1.fq.gz ${RAW}/M28_input_1.fq.gz ${RAW}/M33_input_1.fq.gz \
	--ctrlR2 ${RAW}/M18_input_2.fq.gz ${RAW}/M28_input_2.fq.gz ${RAW}/M33_input_2.fq.gz \
	--genome /home/AD/hrayjones/bin/chip-ap/chipap_installation/chipap_scripts/genomes \
	--output ./chip-ap \
	--setname mono_CTCF \
	--homer_motif both \
	--run

conda deactivate
