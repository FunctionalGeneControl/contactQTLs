#!/bin/bash
#$ -S /bin/sh
#$ -N countsFile
# execute the job from the current working directory:
#$ -cwd
#$ -l cpu=8
# Specify locations for STOUT and STERROR:
#$ -e /data/cmn_vamal/4C/THBS1/4Cker/ER
#$ -o /data/cmn_vamal/4C/THBS1/4Cker/OU
# Tell it where to find conda: (note "source" is substituted with a dot)
. ~/miniconda3/etc/profile.d/conda.sh

### Activate the env with all required tools installed (including bedtools and UCSC oligoMatch)
#conda activate 4Cker
#conda activate deeptools # i need to convert bam to bedgraph
conda activate hicup # using bedtools

### Actually this PATH has oligoMatch
PATH=$PATH:/home/AD/hrayjones/bin
DIR=/data/cmn_vamal/4C/THBS1
RAW=${DIR}/raw
RGENOME=${DIR}/4Cker/reducedGenome

### Creating a counts file from mapped data (SAM file output from bowtie2)

# now trying with the files that were mapped using pipe4C
# with 
cd ${DIR}/4Cker/mapped_with_pipe4C/bam
for file in *.bam; 
do
	coverageBed -a ${DIR}/4Cker/design/dpnII_AseI.rmap -b $file > ../counts/$(echo $file | sed 's/.bam/.counts/g');
done

## this is giving counts per RF. Look at column 1,2,3 and 5
## then we can use these counts to run 4Cker. We remove the adjacent fragments to the bait first. 


conda deactivate
