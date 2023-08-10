#PBS -S /bin/bash
#PBS -N chicago
#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=1:mem=62gb
#PBS -o ~/HRJ_monocytes/CHiC/chicago/OU
#PBS -e ~/HRJ_monocytes/CHiC/chicago/ER
#PBS -J 1-34

##########################################################

#### CHANGE THIS
chicago=~/HRJ_monocytes/CHiC/chicago
cd ${chicago}/combinations_input/bam-files
#ls * | rev | cut -c20- | rev > samples.list

mysample=$(head -n $PBS_ARRAY_INDEX samples.list | tail -n 1)
module load bedtools
module load samtools/1.3.1
module load anaconda3/personal
source activate chicago_env
PATH=$PATH:~/bin/chicago_aug21/chicago/chicagoTools/

# CHECK THE CHICAGO DESIGNDIR AND DESTINATION FOLDERS
DESIGN=~/spivakov/Design

## Make the dpnii-level chinput files
cd ${chicago}/combinations_input/DpnII
bam2chicago_V02.sh -b ${chicago}/combinations_input/bam-files/${mysample}.hicup.captured.bam \
  -t ${DESIGN}/Human_eQTL_CHiC_DpnII_hg38/eCHiC_grch38.baitmap \
  -r ${DESIGN}/Human_eQTL_CHiC_DpnII_hg38/hg38_dpnII.rmap \
  -o ${mysample}_DpnII \
  --nodelete \
  --combinations

# Make the 5K chinput files
cd ${chicago}/combinations_input/5K-bins
bam2chicago_V02.sh -b ${chicago}/combinations_input/bam-files/${mysample}.hicup.captured.bam \
  -t ${DESIGN}/Human_eQTL_CHiC_DpnII_hg38_5K_bins/eCHiC_grch38_5K.baitmap \
  -r ${DESIGN}/Human_eQTL_CHiC_DpnII_hg38_5K_bins/human_DpnII_5K.rmap \
  -o ${mysample}_5K \
  --nodelete \
  --combinations

# Run Chicago

# Using binsize 4k and maxL 100k
cd ${chicago}/combinations_output/DpnII
Rscript ~/bin/chicago_aug21/chicago/chicagoTools/runChicago.R \
 --design-dir ${DESIGN}/Human_eQTL_CHiC_DpnII_hg38_binsize4000_maxL100K \
 --settings-file ${DESIGN}/Human_eQTL_CHiC_DpnII_hg38_binsize4000_maxL100K/mono_dpnII_fraglevel_5reps.settings.txt \
 --en-feat-list ~/ext_datasets/blueprint_projected_segmentations/hg38/blueprint_monocytes_feat_list.txt \
 ${chicago}/combinations_input/DpnII/${mysample}_DpnII/${mysample}_DpnII.chinput \
 ${mysample}_DpnII

# Using standard 5kb settings
cd ${chicago}/combinations_output/5K-bins
Rscript ~/bin/chicago_aug21/chicago/chicagoTools/runChicago.R \
 --design-dir ${DESIGN}/Human_eQTL_CHiC_DpnII_hg38_5K_bins \
 --settings-file ${DESIGN}/Human_eQTL_CHiC_DpnII_hg38_5K_bins/mono_dpnII_5Kbins_5reps.settings.txt \
 --en-feat-list ~/ext_datasets/blueprint_projected_segmentations/hg38/blueprint_monocytes_feat_list.txt \
 ${chicago}/combinations_input/5K-bins/${mysample}_5K/${mysample}_5K.chinput \
 ${mysample}_5K

conda deactivate











