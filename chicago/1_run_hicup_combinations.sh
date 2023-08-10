#PBS -S /bin/bash
#PBS -N make_Scripts
#PBS -l walltime=1:00:00
#PBS -l select=1:ncpus=1:mem=12gb
#PBS -o ~/HRJ_monocytes/CHiC/hicup_combinations/OU
#PBS -e ~/HRJ_monocytes/CHiC/hicup_combinations/ER
#PBS -J 1-16

##########################################################

############ This script makes the permut fastq scripts from a predefined stub. All fastq files were put in their own folder. Then runs the script. After this, run hicup combinations in the same folders.
### HiCUP combinations version 0.7.2

cd ~/HRJ_monocytes/CHiC/hicup_combinations
mysample=$(head -n $PBS_ARRAY_INDEX combine_samples2.list | tail -n 1)

#### CHANGE THIS
cd ${mysample}


### Make the config file
echo "${mysample}_R1.fastq.gz" > ${mysample}_A
echo "${mysample}_R2.fastq.gz" > ${mysample}_B
cat ~/HRJ_monocytes/CHiC/scripts/wrappers/hicup_combinations_config_stub ${mysample}_A ${mysample}_B > ${mysample}.config


### Make the HiCUP combinations sh script
echo "cd ~/HRJ_monocytes/CHiC/hicup_combinations/${mysample}" > ${mysample}_1
echo "hicup -c ${mysample}.config" > ${mysample}_2
echo "conda deactivate" > ${mysample}_3
cat ~/HRJ_monocytes/CHiC/scripts/wrappers/hicup_combinations_sh_stub ${mysample}_1 ${mysample}_2 ${mysample}_3 > ${mysample}_run_hicup_combinations.sh


### Run HiCUP combinations
qsub ${mysample}_run_hicup_combinations.sh


### Clean up
rm ${mysample}_A ${mysample}_B ${mysample}_1 ${mysample}_2 ${mysample}_3

# Afterwards, will need to get capture efficiency










