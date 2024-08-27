#!/bin/bash
#$ -S /bin/sh
#$ -N selfLig
# execute the job from the current working directory:
#$ -cwd
#$ -l cpu=8
# Specify locations for STOUT and STERROR:
#$ -e /data/cmn_vamal/4C/THBS1/4Cker/ER
#$ -o /data/cmn_vamal/4C/THBS1/4Cker/OU
# Tell it where to find conda: (note "source" is substituted with a dot)
. ~/miniconda3/etc/profile.d/conda.sh

### Activate the env with all required tools installed (including bedtools and UCSC oligoMatch)
conda activate 4Cker

### Actually this PATH has oligoMatch
PATH=$PATH:/home/AD/hrayjones/bin
DIR=/data/cmn_vamal/4C/THBS1
RAW=${DIR}/raw
RGENOME=${DIR}/4Cker/reducedGenome

### This is a manual check for the fragments that we want to remove that are above or below our chosen fragment
### You have to use the Watson strand. i.e. if you primer was for the reverse strand you have to use the reverse complement.
### This was for primer TGAGTGATTGCAAATGGAAA

### UPDATE: NOW DOING FOR THE BEDGRAPHS MADE WITH PIPE4C OUTPUT AND THE BEDTOOLS COVERAGE TOOL. extract the relevant regions and make into bedgraph.
### for the SNPs region files
cd ${DIR}/4Cker/design
grep -E -A 1 -B 1 'chr15.39315234' dpnII_AseI.rmap
cd ${DIR}/4Cker/mapped_with_pipe4C/counts


### for the split alleles
for file in $(cat splitReads.counts);
do grep -v -E 'chr15.39315195.39315233|chr15.39315234.39315576|chr15.39315577.39315651' $file | awk 'OFS="\t"{print $1,$2,$3,$5}' > $(echo $file | sed 's/.counts/_rm_self_und.bedGraph/g');
done


#### Earlier notes, now fixed by mapping first with pipe4C
## there's a problem here; the first fragment is never there. I think because it is <122 bp, so how will my reads ever map to it?
## the problem with this tool is no 3' trimming and no identification of another cut site in the reads.
## pipe4C on the other hand deals with this so we get more mapped reads.

conda deactivate
