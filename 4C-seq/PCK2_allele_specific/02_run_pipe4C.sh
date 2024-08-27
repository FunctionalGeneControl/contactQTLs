#$ -S /bin/bash
#$ -N pipe4C_PCK2
## execute the job from the current working directory:
#$ -cwd
#$ -l cpu=8
## Specify locations for STOUT and STERROR:
#$ -e /data/cmn_vamal/4C/PCK2/pipe4C/ER
#$ -o /data/cmn_vamal/4C/PCK2/pipe4C/OU
## Tell it where to find conda:
source ~/miniconda3/etc/profile.d/conda.sh

BASE=/data/cmn_vamal/4C/PCK2/pipe4C
RAW=/data/cmn_vamal/4C/PCK2/raw
conda activate pipe4C

## For this locus it was very important to exclude "fix" and "patch" type chromosomes from the bowtie2 indices, because there is a fix covering 700kb around my primer. 

cd ~/bin/pipe4C


#Rscript pipe4C.R \
#	-v ${BASE}/scripts/VPinfo.txt \
#	-f ${RAW} \
#	-o ${BASE}/results \
#	--cores 8 \
#	--wig \
#	--plot \
#	--genomePlot

Rscript pipe4C.R \
	-v ${BASE}/scripts/VPinfo_splitReads.txt \
	-f ${RAW} \
        -o ${BASE}/results_splitReads \
        --cores 8 \
        --wig \
        --plot \
        --genomePlot


## Allow alnStart = FALSE (edited this in conf.yml)
#Rscript pipe4C.R \
#        -v ${BASE}/scripts/VPinfo.txt \
#        -f ${RAW} \
#        -o ${BASE}/results_alnF \
#        --cores 8 \
#        --wig \
#        --plot \
#        --genomePlot 

# trying with trimming
#Rscript pipe4C.R \
#        -v ${BASE}/scripts/VPinfo.txt \
#        -f ${RAW} \
#        -o ${BASE}/results_trimming75 \
#        --cores 8 \
#        --wig \
#        --plot \
#        --genomePlot \
#	--trimLength 75

conda deactivate
