#PBS -S /bin/bash
#PBS -N chicago_5K
#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=32:mem=62gb
#PBS -o ~/HRJ_monocytes/CHiC/chicago/OU
#PBS -e ~/HRJ_monocytes/CHiC/chicago/ER
#PBS -J 1-3

##########################################################

#### CHANGE THIS
chinput=~/HRJ_monocytes/CHiC/chicago/combinations_input/5K-bins/combined_chinputs
cd ${chinput}
ls merge* | rev | cut -c9- | rev > samples.list
mysample=$(head -n $PBS_ARRAY_INDEX samples.list | tail -n 1)

# CHECK THE CHICAGO DESIGNDIR AND DESTINATION FOLDERS
DESIGNDIR=~/spivakov/Design/Human_eQTL_CHiC_DpnII_hg38_5K_bins_maxL2mb
BAITMAP=~/spivakov/Design/Human_eQTL_CHiC_DpnII_hg38_5K_bins_maxL2mb/eCHiC_grch38_5K.baitmap
RMAP=~/spivakov/Design/Human_eQTL_CHiC_DpnII_hg38_5K_bins_maxL2mb/human_DpnII_5K.rmap
CHINPUT_DIR=~/HRJ_monocytes/CHiC/chicago/combinations_input/5K-bins/combined_chinputs
CHICAGO_OUT=~/HRJ_monocytes/CHiC/chicago/combinations_output/5K-bins


module load bedtools
module load anaconda3/personal
source activate chicago_env

PATH=$PATH:~/bin/chicago_aug21/chicago/chicagoTools/

# Make the chinput files
bam2chicago.sh ${mysample}.captured.bam ${BAITMAP} ${RMAP} ${CHINPUT_DIR}/${mysample} nodelete

# Run Chicago
cd ${CHICAGO_OUT}

Rscript ~/bin/chicago_aug21/chicago/chicagoTools/runChicago.R \
 --design-dir ${DESIGNDIR} \
 --settings-file ${DESIGNDIR}/hindIII_orig_weights.txt \
 --en-feat-list ~/ext_datasets/blueprint_projected_segmentations/hg38/blueprint_monocytes_feat_list.txt \
 ${CHINPUT_DIR}/${mysample}.chinput \
 ${mysample}

conda deactivate











