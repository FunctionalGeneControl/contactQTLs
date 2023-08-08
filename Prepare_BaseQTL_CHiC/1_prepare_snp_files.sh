### Prepare SNP files to run in refBias pipeline

### I am using the VCF file with the alleles updated such that they match the hg38 strand.
cd ~/HRJ_monocytes/genotyping/imputed/hg38_filtered
cp all.rsq3.maf05.hg38.dose.finalRenamed.newIDs.vcf ~/HRJ_monocytes/AS_CHiC/BaseQTL/refbias/all_snps/

### Make a sites file, split by chromosome.
cd ~/HRJ_monocytes/AS_CHiC/BaseQTL/refbias/all_snps
tail -n +49 all.rsq3.maf05.hg38.dose.finalRenamed.newIDs.vcf | cut -f 1-2,4-5 > allsnps.txt # get rid of vcf header
# split by chromosome
awk 'OFS="\t" {print $0 > "chr" $1 ".maf05.hg38strand.txt"}' allsnps.txt

# now remove the chromosome column
for chr in {1..22}; do \
	awk 'OFS="\t"{print $2,$3,$4}' chr${chr}.maf05.hg38strand.txt > temp; \
	mv temp chr${chr}.maf05.hg38strand.txt; \
done

rm allsnps.txt all.rsq3.maf05.hg38.dose.finalRenamed.newIDs.vcf
