#PBS -S /bin/bash
#PBS -N get_AI_ATAC
#PBS -l walltime=8:00:00
#PBS -l select=1:ncpus=32:mem=124gb
#PBS -o ~/HRJ_monocytes/AS_ATAC/BaseQTL/refbias/OU
#PBS -e ~/HRJ_monocytes/AS_ATAC/BaseQTL/refbias/ER

DIR=~/HRJ_monocytes/AS_ATAC/BaseQTL

source activate refBias # Trying again to make a refbias env.
## Uses python 3.9.7, pysam 0.18.0, numpy 1.22.3, pandas 1.4.1, pytables 3.7.0, scipy 1.8.0

module load bowtie2/2.2.9 # the same version as for initial alignment

cd ${DIR}
python ~/HRJ_monocytes/BaseQTL/baseqtl_pipeline/input/refbias/AI.py \
	--initial-AI ./refbias/step1/*.initial.AI.txt \
	--post-remap_AI ./refbias/step3/*.remap.fq1_2.sort.post_remapping_AI.txt \
	--output_file ./AI/AI.txt 


conda deactivate
