#!/bin/bash
#$ -S /bin/sh
#$ -N splitReads
# execute the job from the current working directory:
#$ -cwd
#$ -l cpu=8
# Specify locations for STOUT and STERROR:
#$ -e /data/cmn_vamal/4C/THBS1/ER
#$ -o /data/cmn_vamal/4C/THBS1/OU
# Tell it where to find conda: (note "source" is substituted with a dot)
. ~/miniconda3/etc/profile.d/conda.sh

# seqkit splits reads based on sequence (it has many other features but this one is of interest to us). 
# Please cite: W Shen, S Le, Y Li*, F Hu*. SeqKit: a cross-platform and ultrafast toolkit for FASTA/Q file manipulation. 
# PLOS ONE. doi:10.1371/journal.pone.0163962. 
# I look in the first 30 nt for the primer seq (TGAGTGATTGCAAATGGAAA) followed by the SNP rs2033937 (T or C). Then until the end of the DpnII frag (GATC).
# I first did a count using seqkit -C.
# T vs C allele (rs2033937)
# M4: 717680 vs 699158 reads
# M9: 1009737 vs 990701 reads
# M13: 777090 vs 777629 reads.

cd /data/cmn_vamal/4C/THBS1/raw
~/bin/seqkit grep -R 1:30 -s -p TGAGTGATTGCAAATGGAAATTGAGATC -m 0 M4_THBS1_1.fq.gz -o M4_THBS1_1_T.fq.gz
~/bin/seqkit grep -R 1:30 -s -p TGAGTGATTGCAAATGGAAACTGAGATC -m 0 M4_THBS1_1.fq.gz -o M4_THBS1_1_C.fq.gz
~/bin/seqkit grep -R 1:30 -s -p TGAGTGATTGCAAATGGAAATTGAGATC -m 0 M9_THBS1_1.fq.gz -o M9_THBS1_1_T.fq.gz
~/bin/seqkit grep -R 1:30 -s -p TGAGTGATTGCAAATGGAAACTGAGATC -m 0 M9_THBS1_1.fq.gz -o M9_THBS1_1_C.fq.gz
~/bin/seqkit grep -R 1:30 -s -p TGAGTGATTGCAAATGGAAATTGAGATC -m 0 M13_THBS1_1.fq.gz -o M13_THBS1_1_T.fq.gz
~/bin/seqkit grep -R 1:30 -s -p TGAGTGATTGCAAATGGAAACTGAGATC -m 0 M13_THBS1_1.fq.gz -o M13_THBS1_1_C.fq.gz
