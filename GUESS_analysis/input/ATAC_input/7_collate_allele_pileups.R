### Load in the allelic counts from phASER
### Note that the allelic counts file seem to ignore the blacklisting option and blacklisting is only in the haplo counts files
### We can blacklist manually if we want to
### Per sample, load the AS counts and the per-peak counts.
### SNPs have been assigned to nearest peak. This is fine for within allele, but for total counts needed to generate bed files.

library(data.table)
library(dplyr)
library(tidyr)

# Will need to intersect the snps with the peaks found in this bed file (see script 5_match_SNPs_with_peaks.sh; this bed file used "newIDs"):
setwd("~/spivakov/HRJ_monocytes/AS_ATAC/WASP/mono_34reps_noMultiMapped")
alleles <- fread("./snps/closest_peaks_5kb.bed")
names(alleles) = c("Chr", "hg38Proxy_zero", "hg38Proxy_pos", "hg19SNP_ID", "PeakChr", "Peak_start", "Peak_end", "Peak_name", "Dist2Peak")

##############################################################################

# Make individual files and one large matrix per analysis

###### Testing ######
#file <- "S025NM-01_noBlacklist.allelic_counts.txt"
#mycounts <- fread("./S025NM-01_noBlacklist.allelic_counts.txt")
#snps2peaks <- alleles
#####################

get_counts <- function(counts_list, snps2peaks) { 
  all_counts <- data.table()
  for(file in counts_list) {
    mycounts <- fread(file, sep = "\t")
    names(mycounts)[1] = "Chr"
    names(mycounts)[2] = "hg38Proxy_pos"
    id <- "_noBlacklist.allelic_counts.txt"
    myname <- sub(id, "", file)
    
    ######## SNPs-to-peaks
    setkey(snps2peaks, Chr, hg38Proxy_pos)
    
    with_peaks <- mycounts[snps2peaks, on = c("Chr", "hg38Proxy_pos"), nomatch = NULL]
    with_peaks[, c("totalCount", "hg38Proxy_zero", "variantID", "Peak_name", "Dist2Peak") := NULL]
    with_peaks[, feature1 := paste(Peak_start, Peak_end, sep = "-")]
    with_peaks[, feature := paste(PeakChr, feature1, sep = ":")]
    
    with_peaks2 <- unique(with_peaks[, .(Chr, hg19SNP_ID, feature, hg38Proxy_pos, refAllele, altAllele, refCount, altCount)])
    with_peaks2[, Sample := myname]
    
    all_counts <- rbind(all_counts, with_peaks2, fill= TRUE)
  }
  return(all_counts)
}


setwd("./phASER_out")
counts <- list.files("./", pattern = "\\_noBlacklist.allelic_counts.txt$")
my_result <- get_counts(counts, alleles)

# Should be in the format: feature, SNP, hg38Pos, allele, counts

###### To test with Dplyr

get_matrix <- function(counts) {
  matrix1 <- as.data.table(counts %>% pivot_longer (
    cols = c("refAllele", "altAllele"),
    names_to = "Allele_type",
    values_to = "Allele",
  ))

  matrix1[Allele_type == "refAllele", Counts := refCount]
  matrix1[Allele_type == "altAllele", Counts := altCount]
  matrix2 <- unique(matrix1[, .(feature, hg19SNP_ID, hg38Proxy_pos, Allele, Counts, Sample)])

  matrix3 <- as.data.table(matrix2 %>% pivot_wider (
    id_cols = c("feature", "hg19SNP_ID", "hg38Proxy_pos", "Allele"),
    names_from = Sample,
    values_from = Counts,
    values_fill = 0
  ))

  setnames(matrix3, c("hg19SNP_ID", "hg38Proxy_pos", "Allele"), c("SNP", "hg38SNP_pos", "allele")) ## matches the last AS_ATAC format for Elena
  setorder(matrix3, feature)
  return(matrix3)
}

####
atac_matrix <- get_matrix(my_result)


### How does this compare with the reads that we had before?
old_matrix <- fread("~/HRJ_monocytes/AS_ATAC/WASP/mono_34reps/for_elena/mono34reps_ATAC_allele_counts_noMiss_hg38pos_newIDs_noMulti.txt")
old_matrix[SNP == "1:752566:G:A"]
atac_matrix[SNP == "1:752566:G:A"]
old_matrix[SNP == "1:100503564:T:C"]
atac_matrix[SNP == "1:100503564:T:C"]

######## Write files.

setwd("../allele_counts")
write.table(atac_matrix, file = "mono_34reps_ATAC_allele_counts_phASER.txt",
            quote = FALSE, row.names = FALSE, sep = "\t")


############# phASER does not handle homozygotes so need to think about how to represent them in the table.
############# 

#allele_counts <- fread("./allele_counts/mono_34reps_ATAC_allele_counts_phASER.txt")

################################################################################################
############### Now to get the TOTAL counts per DpnII fragment for those interactions of interest.
################################################################################################

# To get total counts, can use the bed files.
# Need to intersect with peaks/snps, using the file "alleles"
# Writing a script to do this per sample.
# Here, combine those counts into one table.

setwd("../peak_counts")

### Write for all samples
counts_list <- list.files("./", pattern = "\\_total_counts.txt$") 

###### Testing ######
#mycounts <- fread("./S025NM-01_total_counts.txt")
#####################

#### The counts were obtained from the bed files using parallel scripts.

collate_total_counts <- function(counts_list) { 
  all_counts <- data.table()
  for(file in counts_list) {
    mycounts <- fread(file)
    all_counts <- rbind(all_counts, mycounts, fill= TRUE)
  }
  return(all_counts)
}

all_totals <- collate_total_counts(counts_list = counts_list)

####### Make the large matrix

totals_mat <- as.data.table(all_totals %>% pivot_wider (
  id_cols = c("feature", "SNP", "hg38SNP_pos"), 
  names_from = Sample,
  values_from = n,
  values_fill = 0
))

setorder(totals_mat, feature)


######## Write file
getwd()
setwd("~/spivakov/HRJ_monocytes/AS_ATAC/WASP/mono_34reps_noMultiMapped/allele_counts/")
write.table(totals_mat, file = "mono_34reps_ATAC_total_counts.txt", 
            quote = FALSE, row.names = FALSE, sep = "\t")

### Compare against previous - make a plot.

### TODO: check that we have avoided the issues last time with SNP IDs
### Make a file where we collate everything

### The issues last time were:
# 1. Multi-allelic snps - remove them from the counts analysis (in ATAC only, I think - make sure I used the "no multi" vcf for everything)
# 2. Strand issues for the alleles. Make sure that the strand in the SNP ID matches the strand for the allele.

### 1. remove mulli-allelic snps

## Using a txt file written for duplicate SNPs from the VCF with the new IDs (correct hg38 strand alleles)
multi <- fread("~/HRJ_monocytes/genotyping/imputed/hg38_filtered/multialleleic_snps.txt",
              header = FALSE)

atac_matrix <- fread("mono_34reps_ATAC_allele_counts_phASER.txt")
no_multi_alleles <- atac_matrix[!SNP %in% multi$V3]

### Make sure we have exactly two alleles for each snp
to_test <- unique(no_multi_alleles[, .(SNP, allele)])
nrow(to_test) # 2018838
length(unique(to_test[, SNP])) # 1009419 # Good!

### 2.  Make sure the alleles do not have any strand flips. the SNP IDs should match the hg38 strand.
genosnps <- fread("~/HRJ_monocytes/genotyping/imputed/hg38_filtered/all.rsq3.maf05.hg38.dose.finalRenamed.newIDs.sites.vcf")
# This file has the correct for the SNP IDs
genos <- genosnps[, .(ID, REF, ALT, POS)] ###### change from here to add hg38 positions. Use REF and ALT in the SNP name for Elena's analysis.

# Might not need the genos file. Making sure that REF and ALT match the SNP ID here. IF so, we are fine already
smaller <- unique(no_multi_alleles[, .(SNP, hg38SNP_pos, allele)])
x<- rep(c("REF", "ALT"), times = nrow(smaller)/2)
smaller[, "allele_type" := x]

smaller_wide <- as.data.table(pivot_wider(smaller, id_cols = c("SNP", "hg38SNP_pos"),
                            names_from = "allele_type", values_from = "allele"))
smaller_wide[, c("Chr", "hg19pos", "ID_REF", "ID_ALT") := tstrsplit(SNP, split = ":")]
smaller_wide[REF == ID_REF]
smaller_wide[REF != ID_REF]
smaller_wide[ALT == ID_ALT]
smaller_wide[ALT != ID_ALT] # OK, we are good!

# Make sure that the alleles match the IDs, just in case.
setnames(smaller_wide, c("REF", "ALT"), c("oldREF", "oldALT"))
with_genos <- smaller_wide[genos, on = c(SNP = "ID", hg38SNP_pos = "POS"), nomatch = NULL]
with_genos[REF != oldREF]
with_genos[REF != ID_REF]
with_genos[ALT != oldALT]
with_genos[ALT != ID_ALT]

## Check the couple that Elena flagged
no_multi_alleles[SNP %like% "10:51566879"]
no_multi_alleles[SNP %like% "10:51568378"]

no_multi_alleles[SNP %like% "22:50468907:C:G"]
no_multi_alleles[SNP %like% "10:20320888:T:A"] # should not exist any more

#### Save the file with no multiallelic snps.
fwrite(no_multi_alleles, file = "./mono_34reps_ATAC_allele_counts_WASP_phASER_noMulti_hg38Strand.txt",
       sep = "\t", quote = FALSE, row.names = FALSE)

### Also remove multialleleic in the total counts
total_counts <- fread("./mono_34reps_ATAC_total_counts.txt")
no_multi_totals <- total_counts[!SNP %in% multi$V3]

### Double check the IDs are correct
no_multi_totals[!SNP %in% genos$ID]

#### Save
fwrite(no_multi_totals, file = "./mono_34reps_ATAC_total_counts_WASP_noMulti_hg38Strand.txt",
       sep = "\t", quote = FALSE, row.names = FALSE)

############# Will do same for CHiC

#### Double check them
alleles <- fread("./mono_34reps_ATAC_allele_counts_WASP_phASER_noMulti_hg38Strand.txt")
totals <- fread("./mono_34reps_ATAC_total_counts_WASP_noMulti_hg38Strand.txt")
####

### Compare against previous - make a plot
old <- fread("~/spivakov/HRJ_monocytes/AS_ATAC/BaseQTL/for_elena/mono_34reps_ATAC_total_counts_noMulti_hg38Strand.txt")
new <- copy(totals)

length(unique(old$SNP))
length(unique(new$SNP))

#old[SNP == "10:101760700:G:A"]
#new[SNP == "10:101760700:G:A"]

old[, SNP_feat := paste(feature, SNP, sep = "_")]
new[, SNP_feat := paste(feature, SNP, sep = "_")]

# Avg. no reads across samples
mySamples <- names(old[, 4:37])

old[, avg_reads := rowMeans(.SD), .SDcols = c(mySamples)]
new[, avg_reads := rowMeans(.SD), .SDcols = c(mySamples)]

mySNPfeats <- unique(old$SNP_feat)
mySNPfeatsSampled <- sample(mySNPfeats, 500)

#print(mySNPfeatsSampled)

old_sampled <- old[SNP_feat %in% mySNPfeatsSampled]
new_sampled <- new[SNP_feat %in% mySNPfeatsSampled]

setnames(old_sampled, "avg_reads", "avg_reads_noCorrection")
setnames(new_sampled, "avg_reads", "avg_reads_WASPCorrection")

old_sampled_plot <- unique(old_sampled[, .(SNP_feat, avg_reads_noCorrection)])
new_sampled_plot <- unique(new_sampled[, .(SNP_feat, avg_reads_WASPCorrection)])

to_plot <- old_sampled_plot[new_sampled_plot, on = "SNP_feat", nomatch = NULL]

print(to_plot)

library(ggplot2)
pdf(file = "./ATAC_WASP_comparison_with_no_correction.pdf")
p <- ggplot(to_plot, aes(x = avg_reads_noCorrection, y = avg_reads_WASPCorrection)) + geom_point() +
        ggtitle("ATAC-seq WASP correction, 500 samples")
print(p)
dev.off()


