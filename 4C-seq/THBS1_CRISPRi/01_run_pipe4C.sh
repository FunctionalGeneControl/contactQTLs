#$ -S /bin/bash
#$ -N pipe4C_THBS1
## execute the job from the current working directory:
#$ -cwd
#$ -l cpu=8
## Specify locations for STOUT and STERROR:
#$ -e /data/cmn_vamal/4C/THBS1/pipe4C/ER
#$ -o /data/cmn_vamal/4C/THBS1/pipe4C/OU
## Tell it where to find conda:
source ~/miniconda3/etc/profile.d/conda.sh

BASE=/data/cmn_vamal/4C/THBS1/pipe4C
RAW=/data/cmn_vamal/4C/THBS1/raw
conda activate pipe4C

## Not using a trimlength (previously had set in the conf.yml file trimLength : 75. 
cd ~/bin/pipe4C

Rscript pipe4C.R \
	-v ${BASE}/scripts/VPinfo.txt \
	-f ${RAW} \
	-o ${BASE}/results_with_bigwigs_all \
	--cores 8 \
	--bigwig \
	--plot \
	--genomePlot \
	--tsv \
	--bins

# trying with trimming 75bp
#Rscript pipe4C.R \
#        -v ${BASE}/scripts/VPinfo.txt \
#        -f ${RAW} \
#        -o ${BASE}/results_trimming75 \
#	--cores 8 \
#        --wig \
#        --plot \
#        --genomePlot \
#        --trimLength 75


conda deactivate
