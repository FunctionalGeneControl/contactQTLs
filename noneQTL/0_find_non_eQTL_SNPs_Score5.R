#### Look for control interactions between non eQTLs and gene promoters at chicago score 5.
#### For this, can filter the baseQTL results directly.
#### We can only keep results for which ALL the interactions occurred at a score of 5.

library(data.table)
library(dplyr)
setwd("~/HRJ_monocytes/AS_CHiC/BaseQTL/scripts/non_eqtls/")

# Load DpnII-level interactions at score 5
ints <- fread("~/HRJ_monocytes/CHiC/chicago/combinations_output/DpnII/combined_34reps_DpnII_binsize1500_maxL150K_mergeWeights/data/combined_34reps_DpnII_binsize1500_maxL150K_mergeWeights_washU_text.txt", 
              sep = "\t")
intsb <- washu2bedpe(ints)

# Load the input for non-eQTLs
non1 <- fread("../non_eqtls_non_gwas_otherEnd.bedpe")
non2 <- fread("../non_eqtls_gwas_otherEnd.bedpe")
bothnon <- rbind(non1, non2)
names(bothnon) = c("chr1", "start1", "end1", "chr2", "start2", "end2", "feature")

## Ask if the interactions are in score5
# for this do an antijoin
# put the lower fragment on the left, each time
bothnon[start1 < start2, `:=` (leftchr = chr1, leftstart = start1, leftend = end1, rightchr = chr2, rightstart = start2, rightend = end2)]
bothnon[start1 > start2, `:=` (leftchr = chr2, leftstart = start2, leftend = end2, rightchr = chr1, rightstart = start1, rightend = end1)]
bothnon[, leftchrom := paste0("chr", leftchr)]
bothnon[, rightchrom := paste0("chr", rightchr)]
bothnon[, c("leftchr", "rightchr") := NULL]
setnames(bothnon, c("leftchrom", "rightchrom"), c("leftchr", "rightchr"))

intsb[start1 < start2, `:=` (leftchr = chr1, leftstart = start1, leftend = end1, rightchr = chr2, rightstart = start2, rightend = end2)]
intsb[start1 > start2, `:=` (leftchr = chr2, leftstart = start2, leftend = end2, rightchr = chr1, rightstart = start1, rightend = end1)]

score3_match <- bothnon[, .(leftchr, leftstart, leftend, rightchr, rightstart, rightend, feature)]
score5_match <- intsb[, .(leftchr, leftstart, leftend, rightchr, rightstart, rightend)]
score5_match[, `:=` (leftstart = as.numeric(leftstart), leftend = as.numeric(leftend), rightstart = as.numeric(rightstart), rightend = as.numeric(rightend))]

# now do the antijoin
not_score5 <- score3_match[!score5_match, on = c("leftchr", "leftstart", "leftend", "rightstart", "rightend")]

keep <- bothnon[!feature %in% not_score5$feature]

keep[, SNP := tstrsplit(feature, split = "_", keep = 1)]
length(unique(keep$SNP)) # 14,779 SNPs.

# write file and use to filter the BaseQTL results.
fwrite(keep, file = "../non_eqtls_score5_keep.txt", sep = "\t", quote = F, row.names = F, col.names = T)




