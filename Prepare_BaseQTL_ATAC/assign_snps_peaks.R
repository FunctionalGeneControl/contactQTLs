######## Get an input file  with all genotyped SNP positions.
######## Assign these snps to ATAC-seq peaks, in order to do the allele counting.
### Chr, pos, allele, readgroup
### Check how many ambiguous we get, by the way... it might be worse than with CHiC

library(data.table)
library(dplyr)
library(tidyr)

setwd("~/HRJ_monocytes/AS_ATAC/BaseQTL/snps")
## Note - not using a bim file because it lists minor and major alleles, instead of ref and alt
genosnps <- fread("./all.rsq3.maf05.hg38.dose.finalRenamed.newIDs.sites.vcf")
genosnps[, c("QUAL", "FILTER", "INFO") := NULL]
genosnps2 <- as.data.table(pivot_longer(genosnps, c("REF", "ALT"), names_to = "Genotype", values_to = "Allele"))
genosnps3 <- unique(genosnps2[, .(`#CHROM`, POS, Allele, Genotype)])

# Not using biostar any more
#write.table(genosnps3, file = "all_snp_pos_rsq3_maf05_for_biostar.txt", sep = "\t", 
           # quote = FALSE, row.names = FALSE, col.names = FALSE)

######## Assign snps to ATAC-seq peaks 
### Will use the gappedPeak file from HMMRATAC. Note that gappedPeak files are zero based.
peaks <- fread("~/HRJ_monocytes/ATAC-seq/analyses/HMMRATAC/consensus/output/ATAC_merged_34reps_HMMRATAC_withBlacklist_peaks.gappedPeak")
names(peaks) = c("Chrom", "start", "end", "name", "score", "strand", "thickStart", "thickEnd", "itemRgb", 
                 "blockCount", "blockSizes", "blockStarts", "signalValue", "pvalue", "qvalue")
peaks[, len := end-start]
peaks[len > 50000] # there is one that is 80kb! Keep it in for now.

### We don't want peaks to overlap, necessarily. Can we assign snps to the nearest peak? and discard, if it is >5kb.
### Lets do this using bedtools closest. 
### For that, make bed files of 1) snps and 2) ATAC-seq peaks.
options(scipen = 999) # Disable scientific notation
snpsbed <- copy(genosnps2)
snpsbed[, start := as.numeric(POS-1)]
snpsbed2 <- unique(snpsbed[, .(`#CHROM`, start, POS, ID)])
setorder(snpsbed2, `#CHROM`, start)
write.table(snpsbed2, file = "all_snp_pos_rsq3_maf05.bed", sep = "\t", 
            quote = FALSE, row.names = FALSE, col.names = FALSE)

peaks
peaks2 <- unique(peaks[, .(Chrom, start, end, name)])
setorder(peaks2, Chrom, start)
write.table(peaks2, file = "HMMRATAC_peaks_mono34reps_withBlacklist.bed", sep = "\t", 
            quote = FALSE, row.names = FALSE, col.names = FALSE)

# Now run the script 4_match_SNPs_with_peaks.sh
