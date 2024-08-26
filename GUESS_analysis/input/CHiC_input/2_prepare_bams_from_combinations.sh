#PBS -S /bin/bash
#PBS -N mod_bams
#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=2:mem=12gb
#PBS -o ~/HRJ_monocytes/AS_CHiC/WASP/mono_34reps_allReads_using_eGene_filters/OU
#PBS -e ~/HRJ_monocytes/AS_CHiC/WASP/mono_34reps_allReads_using_eGene_filters/ER
#PBS -J 1-34


BAMDIR=~/HRJ_monocytes/AS_CHiC/WASP/mono_34reps_allReads_using_eGene_filters/bams
cd $BAMDIR

mysample=$(head -n $PBS_ARRAY_INDEX samples.list | tail -n 1)

module load bedtools
module load samtools/1.3.1

mkdir -p ${mysample}

samtools view ${mysample}.filt.captured.bam | cut -f 1 > ${mysample}/${mysample}_qnames
bamcolno=`awk -F: '{print NF; exit}' ${mysample}/${mysample}_qnames`
if [ $bamcolno -gt 6 ]; then # makes sure that you have more than 6 fields
       echo "Removing excess fields in qnames made by HiCUP Combinations... "
       awk -F: '{for(i=0;++i<=NF-7;)printf $i":";print $(NF-6)}' ${mysample}/${mysample}_qnames > ${mysample}/${mysample}_newnames1
       sed 'N;s/\n/\t/' ${mysample}/${mysample}_newnames1 | cat -n | awk 'OFS="\t"{print $2"_"$1"\n"$3"_"$1}' > ${mysample}/${mysample}_newnames # Gives each pair an index, because some pairs will now have the same readname but represent different ditag permutations
       samtools view ${mysample}.filt.captured.bam | cut -f 2- > ${mysample}/${mysample}_torename
       samtools view -H ${mysample}.filt.captured.bam > ${mysample}/${mysample}_header
       paste ${mysample}/${mysample}_newnames ${mysample}/${mysample}_torename > ${mysample}/${mysample}_temp
       cat ${mysample}/${mysample}_header ${mysample}/${mysample}_temp > ${mysample}/${mysample}_2.sam
       samtools view -h -b ${mysample}/${mysample}_2.sam -o ${mysample}/${mysample}.bam
       rm ${mysample}/${mysample}_qnames ${mysample}/${mysample}_newnames1 ${mysample}/${mysample}_newnames ${mysample}/${mysample}_header ${mysample}/${mysample}_torename ${mysample}/${mysample}_temp ${mysample}/${mysample}_2.sam
else
       echo "No excess fields detected in bam readnames. Are you sure that you used HiCUP Combinations?"
fi

