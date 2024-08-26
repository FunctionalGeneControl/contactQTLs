#PBS -S /bin/bash
#PBS -N maxatac_predict
#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=100:mem=800gb
#PBS -o ~/HRJ_monocytes/ATAC-seq/maxatac/OU
#PBS -e ~/HRJ_monocytes/ATAC-seq/maxatac/ER
#PBS -J 1-128

DIR=~/HRJ_monocytes/ATAC-seq/maxatac

### RE conda environment:
### Followed advice re pyBigWig and numpy: "pip uninstall pyBigWig" then "conda install pyBigWig"

source activate maxATAC
module load cuda

### Merge all the normalised bigwigs into one file and then run predict. This was run for all 34 samples.

maxatac average --bigwigs */*_IS_slop20_RP20M_minmax01.bw --prefix mono_normalized_average --output ./average 

###########
cd ~/opt/maxatac/data/models
TF=$(head -n $PBS_ARRAY_INDEX TF.list | tail -n 1)
###########

## Run against all TFs in the catalog.
cd ${DIR}/output
mkdir -p $TF
for chr in {1..22}
do
	echo $TF
	echo $chr
	maxatac predict --sequence ~/opt/maxatac/data/hg38/hg38.2bit -tf $TF --signal average/mono_normalized_average.bw -o ./${TF}/${TF}.averaged.chr${chr} --chromosomes chr${chr} --batch_size 2000
done

conda deactivate
