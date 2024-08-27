#!/bin/bash
#$ -S /bin/sh
#$ -N splitReads
# execute the job from the current working directory:
#$ -cwd
#$ -l cpu=8
# Specify locations for STOUT and STERROR:
#$ -e /data/cmn_vamal/4C/PCK2/ER
#$ -o /data/cmn_vamal/4C/PCK2/OU
# Tell it where to find conda: (note "source" is substituted with a dot)
. ~/miniconda3/etc/profile.d/conda.sh

# seqkit splits reads based on sequence (it has many other features but this one is of interest to us). 
# Please cite: W Shen, S Le, Y Li*, F Hu*. SeqKit: a cross-platform and ultrafast toolkit for FASTA/Q file manipulation. 
# PLOS ONE. doi:10.1371/journal.pone.0163962. 
# I look in the first ~95 nt for the primer seq (ACCAGACTTTTCCACCAGAG) followed by the SNP rs7146599 (G or A). Then until the end of the DpnII frag (GATC).
# I first did a count using seqkit -C.
# G vs A allele (rs7146599)
# M4: G = 1213647, A = 1394331
# M9: G = 1002059, A = 1204763
# M13: 

cd /data/cmn_vamal/4C/PCK2/raw
#~/bin/seqkit grep -R 1:95 -s -p ACCAGACTTTTCCACCAGAGGGCAGGAGCAAGCTGTCTTAGGAATAGCAGCCTTCCTGGCCGAAGGGAGGGAAGCCAGCTCTAGGGAAGGATC -m 0 M4_PCK2_1.fq.gz -o M4_PCK2_1_G.fq.gz
#
#~/bin/seqkit grep -R 1:95 -s -p ACCAGACTTTTCCACCAGAGAGCAGGAGCAAGCTGTCTTAGGAATAGCAGCCTTCCTGGCCGAAGGGAGGGAAGCCAGCTCTAGGGAAGGATC -m 0 M4_PCK2_1.fq.gz -o M4_PCK2_1_A.fq.gz
#
#~/bin/seqkit grep -R 1:95 -s -p ACCAGACTTTTCCACCAGAGGGCAGGAGCAAGCTGTCTTAGGAATAGCAGCCTTCCTGGCCGAAGGGAGGGAAGCCAGCTCTAGGGAAGGATC -m 0 M9_PCK2_1.fq.gz -o M9_PCK2_1_G.fq.gz
#
#~/bin/seqkit grep -R 1:95 -s -p ACCAGACTTTTCCACCAGAGAGCAGGAGCAAGCTGTCTTAGGAATAGCAGCCTTCCTGGCCGAAGGGAGGGAAGCCAGCTCTAGGGAAGGATC -m 0 M9_PCK2_1.fq.gz -o M9_PCK2_1_A.fq.gz
#

~/bin/seqkit grep -R 1:95 -s -p ACCAGACTTTTCCACCAGAGGGCAGGAGCAAGCTGTCTTAGGAATAGCAGCCTTCCTGGCCGAAGGGAGGGAAGCCAGCTCTAGGGAAGGATC -m 0 M13_PCK2_1.fq.gz -o M13_PCK2_1_G.fq.gz

~/bin/seqkit grep -R 1:95 -s -p ACCAGACTTTTCCACCAGAGAGCAGGAGCAAGCTGTCTTAGGAATAGCAGCCTTCCTGGCCGAAGGGAGGGAAGCCAGCTCTAGGGAAGGATC -m 0 M13_PCK2_1.fq.gz -o M13_PCK2_1_A.fq.gz

