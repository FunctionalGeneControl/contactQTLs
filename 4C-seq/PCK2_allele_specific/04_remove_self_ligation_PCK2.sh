#!/bin/bash
#$ -S /bin/sh
#$ -N selfLig
# execute the job from the current working directory:
#$ -cwd
#$ -l cpu=8
# Specify locations for STOUT and STERROR:
#$ -e /data/cmn_vamal/4C/PCK2/4Cker/ER
#$ -o /data/cmn_vamal/4C/PCK2/4Cker/OU
# Tell it where to find conda: (note "source" is substituted with a dot)
. ~/miniconda3/etc/profile.d/conda.sh

### Activate the env with all required tools installed (including bedtools and UCSC oligoMatch)
conda activate 4Cker

### Actually this PATH has oligoMatch
PATH=$PATH:/home/AD/hrayjones/bin
DIR=/data/cmn_vamal/4C/PCK2
RAW=${DIR}/raw
RGENOME=${DIR}/4Cker/reducedGenome

### This is a manual check for the fragments that we want to remove that are above or below our chosen fragment
### You have to use the Watson strand. i.e. if you primer was for the reverse strand you have to use the reverse complement.
### This was for primer ACCAGACTTTTCCACCAGAG 



### UPDATE: NOW DOING FOR THE BEDGRAPHS MADE WITH PIPE4C OUTPUT AND THE BEDTOOLS COVERAGE TOOL. extract the relevant regions and make into bedgraph.
### for the SNPs region files
cd ${DIR}/4Cker/design
grep -E -A 1 -B 1 'chr14.24058271.24058531' dpnII_hpych4v.rmap
cd ${DIR}/4Cker/mapped_with_pipe4C/counts

ls *.counts > splitReads.counts
for file in $(cat splitReads.counts);
do grep -v -E 'chr14.24058261.24058270|chr14.24058271.24058531|chr14.24058532.24058571' $file | awk 'OFS="\t"{print $1,$2,$3,$5}' > $(echo $file | sed 's/.counts/_rm_self_und.bedGraph/g');
done

conda deactivate
