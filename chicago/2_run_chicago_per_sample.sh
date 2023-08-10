#PBS -S /bin/bash
#PBS -N chicago
#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=32:mem=62gb
#PBS -o ~/HRJ_monocytes/CHiC/chicago/OU
#PBS -e ~/HRJ_monocytes/CHiC/chicago/ER
#PBS -J 1-34

##########################################################

#### CHANGE THIS
hicup=~/HRJ_monocytes/CHiC/aligned/aug21_chicago_input
cd ${hicup}
ls *.bam | rev | cut -c14- | rev > samples.list
mysample=$(head -n $PBS_ARRAY_INDEX samples.list | tail -n 1)

# CHECK THE CHICAGO DESIGNDIR AND DESTINATION FOLDERS
DESIGNDIR=~/spivakov/Design/Human_eQTL_CHiC_DpnII_hg38_binsize4000_maxL100K
BAITMAP=~/spivakov/Design/Human_eQTL_CHiC_DpnII_hg38_binsize4000_maxL100K/eCHiC_grch38.baitmap
RMAP=~/spivakov/Design/Human_eQTL_CHiC_DpnII_hg38_binsize4000_maxL100K/hg38_dpnII.rmap
CHINPUT_DIR=~/HRJ_monocytes/CHiC/chicago/input/dpnII/aug21_chicago_input
CHICAGO_OUT=~/HRJ_monocytes/CHiC/chicago/output/dpnII/aug21_chicago_output


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
 --settings-file ${DESIGNDIR}/mono_dpnII_fraglevel_5reps.settings.txt \
 --en-feat-list ~/ext_datasets/blueprint_projected_segmentations/hg38/blueprint_monocytes_feat_list.txt \
 ${CHINPUT_DIR}/${mysample}/${mysample}.chinput \
 ${mysample}

conda deactivate











