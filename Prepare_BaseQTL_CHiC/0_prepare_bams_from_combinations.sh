#PBS -S /bin/bash
#PBS -N mod_bams
#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=2:mem=12gb
#PBS -o ~/HRJ_monocytes/AS_CHiC/BaseQTL/OU
#PBS -e ~/HRJ_monocytes/AS_CHiC/BaseQTL/ER
#PBS -J 1-34

### We have decided to run the refbias analysis on the de-duplicated bam files.
### Have to modify the readnames in these files so that pairs can be found.

BAMDIR=~/HRJ_monocytes/CHiC/chicago/combinations_input/bam-files
cd $BAMDIR

mysample=$(head -n $PBS_ARRAY_INDEX samples.list | tail -n 1)

source activate samtools1.3.1 # Older version of samtools does not throw header error

mkdir -p ${mysample}

samtools view ${mysample}.hicup.captured.bam | cut -f 1 > ${mysample}/${mysample}_qnames
bamcolno=`awk -F: '{print NF; exit}' ${mysample}/${mysample}_qnames`
if [ $bamcolno -gt 6 ]; then # makes sure that you have more than 6 fields
       echo "Removing excess fields in qnames made by HiCUP Combinations... "
       awk -F: '{for(i=0;++i<=NF-7;)printf $i":";print $(NF-6)}' ${mysample}/${mysample}_qnames > ${mysample}/${mysample}_newnames1
       sed 'N;s/\n/\t/' ${mysample}/${mysample}_newnames1 | cat -n | awk 'OFS="\t"{print $2"_"$1"\n"$3"_"$1}' > ${mysample}/${mysample}_newnames # Gives each pair an index, because some pairs will now have the same readname but represent different ditag permutations
       samtools view ${mysample}.hicup.captured.bam | cut -f 2- > ${mysample}/${mysample}_torename
       samtools view -H ${mysample}.hicup.captured.bam | grep -v "HiCUP Deduplicator" > ${mysample}/${mysample}_header # Remove HiCUP deduplicator lines, sometimes causing an error with later versions of samtools
       paste ${mysample}/${mysample}_newnames ${mysample}/${mysample}_torename > ${mysample}/${mysample}_temp
       cat ${mysample}/${mysample}_header ${mysample}/${mysample}_temp > ${mysample}/${mysample}_2.sam
       samtools view -h -b ${mysample}/${mysample}_2.sam -o ${mysample}/${mysample}.bam
       rm ${mysample}/${mysample}_qnames ${mysample}/${mysample}_newnames1 ${mysample}/${mysample}_newnames ${mysample}/${mysample}_header ${mysample}/${mysample}_torename ${mysample}/${mysample}_temp ${mysample}/${mysample}_2.sam
else
       echo "No excess fields detected in bam readnames. Are you sure that you used HiCUP Combinations?"
fi

### Link to the BaseQTL directory
cd ~/HRJ_monocytes/AS_CHiC/BaseQTL/bams_dedup
ln -s ${BAMDIR}/${mysample}/${mysample}.bam ./

conda deactivate
