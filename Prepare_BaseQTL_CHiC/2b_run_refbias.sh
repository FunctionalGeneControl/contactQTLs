#PBS -S /bin/bash
#PBS -N get_AI_CHiC
#PBS -l walltime=8:00:00
#PBS -l select=1:ncpus=32:mem=124gb
#PBS -o ~/HRJ_monocytes/AS_CHiC/BaseQTL/refbias/OU
#PBS -e ~/HRJ_monocytes/AS_CHiC/BaseQTL/refbias/ER

DIR=~/HRJ_monocytes/AS_CHiC/BaseQTL

source activate refBias # Trying again to make a refbias env.
## Uses python 3.9.7, pysam 0.18.0, numpy 1.22.3, pandas 1.4.1, pytables 3.7.0, scipy 1.8.0
cd $DIR
module load bowtie2/2.1.0

### 4. Run the AI script
python ~/HRJ_monocytes/BaseQTL/baseqtl_pipeline/input/refbias/AI.py \
	--initial-AI ./refbias/step1/*.initial.AI.txt \
	--post-remap_AI ./refbias/step3/*.remap.fq1_2.pair.sort.post_remapping_AI.txt \
	--output_file ./AI/AI.txt 


conda deactivate
