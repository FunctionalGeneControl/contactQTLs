#### Compare the non-eQTL (GWAS/non GWAS) BaseQTL results with the eQTL BaseQTL results
## 1. Purge the non-eQTL data of DpnII efects by looking for SNPs in LD and then comparing against DpnII cut sites
## 2. Filter for Rhat and AI
## 3. Generate graphs of effects

library(Biostrings)
library(data.table)
library(stringr)
library(dplyr)
library(tidyr)
library(FSA)
library(ggplot2)
library(EnvStats)


#nongwas <- fread("~/HRJ_monocytes/AS_CHiC/BaseQTL/scripts/non_eqtls_non_gwas_otherEnd.snps") # do not use the ones not saying otherEnd!
#gwas <- fread("~/HRJ_monocytes/AS_CHiC/BaseQTL/scripts/non_eqtls_gwas_otherEnd.snps")
## combine them for this.
#both <- rbind(nongwas, gwas)

# for LD needs to be files with lists of SNP IDs, separated by chromosome.

#for(chr in 1:22) {
#  my_chr <- both[Chr == chr]
#  towrite <- unique(my_chr[, .(hg19ID)])
#  myname <- paste0("GWAS_and_nonGWAS_noneQTLs_otherEnd_IDs_chr", chr, ".txt")
#  fwrite(towrite, file = paste0("~/HRJ_monocytes/AS_CHiC/BaseQTL/scripts/non_eqtls/LD/", myname), sep = "\t", 
#         quote = F, row.names = F, col.names = F)
#}
#list.files("~/HRJ_monocytes/AS_CHiC/BaseQTL/scripts/non_eqtls/LD/")

# now run script get_LD_noneQTLs.sh

######## begin from here!! 21/03/23

### Read in SNPs and get hg38 positions
#setwd("~/HRJ_monocytes/AS_CHiC/BaseQTL/scripts/non_eqtls/LD")
#my_ld <- list.files("./", pattern = "*.maf5.ld$")
#all_snps <- data.table()
#for(file in my_ld) {
#  mychr <- fread(file)
#  all_snps <- rbind(all_snps, mychr)
#}

#nongwas[, type := "nongwas"]
#gwas[, type := "gwas"]
#both_type <- rbind(nongwas, gwas)

#all_snps_pos <- both_type[all_snps, on = c(Chr = "CHR_A", hg19ID = "SNP_A"), nomatch = NULL]
#all_snps_pos[Chr != CHR_B] # none

# Now need to ask, if any of SNP_B affect DpnII sites and are in same frag as SNP A, then we remove SNP A.
# copying and adapting code from /rds/general/project/lms-spivakov-analysis/live/HRJ_monocytes/findmotifs/scripts/get_dpnII_new_sites_CHi-C_SNPs_allProxies.R
# to run the below, qsub find_dpnii_motifs.sh

#myseq <- fread("/rds/general/project/lms-spivakov-analysis/live/HRJ_monocytes/findmotifs/snp_fasta/ALL_genotyped_SNPs_REF_40bp_fasta_with_locations.txt")
#setkey(myseq, SNP, Chr) # Going to join to the proxies, i.e. "SNP_B"
#sets_seqs1 <- all_snps_pos[myseq, on = c(SNP_B = "SNP", "Chr"), nomatch = NULL]
#sets_seqs1 <- sets_seqs1[CHR_B == Chr]
#sets_seqs <- unique(sets_seqs1[, .(hg19ID, SNP_B, Chr, hg38SNP_loc, REF, ALT, Sequence, seq_start, seq_stop)])


## Change the genotypes if the SNPs are in the same SNP set. 
## Now requires that the alleles are listed as REF and ALT
mutate_seqs_nt <- function(myseq) {
  final_refalt <- data.table()
  set_list <- unique(as.list(myseq$hg19ID))
  for(SET in set_list) { # Split by the SNP set (these are the SNPs we care about btw!)
    my_set <- myseq[hg19ID == SET]
    snps <- unique(my_set[, .(SNP_B, REF, ALT, Chr, hg38SNP_loc, hg19ID)]) ### hg19ID is the SNP to filter, SNP_B is the one to check.
    snps[, hg38SNP_loc2 := hg38SNP_loc]
    locations <- unique(my_set[, .(Sequence, hg38SNP_loc, Chr, seq_start, seq_stop)])
    locations[, mod_seq_start := seq_start + 1]
    setkey(locations, Chr, mod_seq_start, seq_stop)
    
    ### The following overlap will get all the snps in LD within the 40bp fasta region
    snp_overlap <- foverlaps(snps, locations, by.x = c("Chr", "hg38SNP_loc", "hg38SNP_loc2"), nomatch = NULL, type = "within")
    snp_overlap[, seqPos := i.hg38SNP_loc - seq_start] # used to be seq_start + 1 but I had the positions one off.
    # Used i.hg38SNP_loc bc this is the position of the SNP in question (whereas hg38Proxy_pos is just the middle SNP.)
    seqtable <- unique(snp_overlap[, .(Sequence, SNP_B, REF, ALT, seqPos, hg19ID)])
    
    
    sequences <- as.list(unique(seqtable$Sequence)) # Get all the sequences for this set
    refalt <- data.table()
    for(s in sequences) {
      mytable <- seqtable[Sequence == s]
      mysequence_REF <- DNAString(unique(mytable$Sequence))
      i <- 0
      repeat{                             # Start
        i <- i + 1                        # Update running index
        if(i > nrow(mytable)) {           # Break condition
          break
        }
        mysequence_REF <- replaceLetterAt(mysequence_REF, mytable[i, seqPos], mytable[i, REF]) # Change the base at each position to REF
      }
      mysequence_ALT <- DNAString(unique(mytable$Sequence))
      i <- 0
      repeat{                             # Start
        i <- i + 1                        # Update running index
        if(i > nrow(mytable)) {           # Break condition
          break
        }
        mysequence_ALT <- replaceLetterAt(mysequence_ALT, mytable[i, seqPos], mytable[i, ALT]) # Change the base at each position to ALT
      }
      mysequence_id <- paste(c(mytable$SNP_B), collapse = "_") # per sequence, get all the snp IDs pasted
      
      mysequence_REF_char <- as.character(mysequence_REF)
      mysequence_ALT_char <- as.character(mysequence_ALT)
      
      current_set <- as.character(unique(mytable$hg19ID))
      
      # get a column with all the snps in that particular sequence, and keep the set ID
      seqs <- data.table(REFSeq = mysequence_REF_char, ALTSeq = mysequence_ALT_char, id = mysequence_id, set = current_set)
      refalt <- rbind(refalt, seqs, fill = TRUE)
    }
    #return(refalt)
    final_refalt <- unique(rbind(final_refalt, refalt, fill = TRUE))
  }
  return(final_refalt)
}


#### Run the mutate function 
#eqs_mutated <- mutate_seqs_nt(sets_seqs)


#### write the results
#setwd("/rds/general/project/lms-spivakov-analysis/live/HRJ_monocytes/AS_CHiC/BaseQTL/scripts/non_eqtls/dpnII_sites")
#fwrite(eqs_mutated, file = "gwas_and_nongwas_noneQTLs_otherEnd_LD.REFALTseqs.txt", 
#       row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

### then carry on and look at dpnii efects as in script ~/HRJ_monocytes/AS_CHiC/BaseQTL/scripts/non_eqtls/compare_eQTL_noneQTL.R - delete that file to avoid confusion!
setwd("~/HRJ_monocytes/AS_CHiC/BaseQTL/scripts/non_eqtls/dpnII_sites")
new_seqs <- fread("gwas_and_nongwas_noneQTLs_otherEnd_LD.REFALTseqs.txt")
setnames(new_seqs, c("id", "set"), c("hg19Proxy_ID", "hg19SNP_ID"))

## In fact here, look in both the REF and the ALT.
new_seqs[, REFseq1 := substr(REFSeq, 18, 21)]
new_seqs[, REFseq2 := substr(REFSeq, 19, 22)]
new_seqs[, REFseq3 := substr(REFSeq, 20, 23)]
new_seqs[, REFseq4 := substr(REFSeq, 21, 24)]
### Any of these four seqs, if containing GATC, have the ref SNP in the cut site for the reference allele
new_seqs[REFseq1 == "GATC" | REFseq2 == "GATC" | REFseq3 == "GATC" | REFseq4 == "GATC", Dpn_site_REF := TRUE]
new_seqs[Dpn_site_REF == TRUE]

new_seqs[, ALTseq1 := substr(ALTSeq, 18, 21)]
new_seqs[, ALTseq2 := substr(ALTSeq, 19, 22)]
new_seqs[, ALTseq3 := substr(ALTSeq, 20, 23)]
new_seqs[, ALTseq4 := substr(ALTSeq, 21, 24)]
### Any of these four seqs, if containing GATC, have the ref SNP in the cut site for the alternative allele
new_seqs[ALTseq1 == "GATC" | ALTseq2 == "GATC" | ALTseq3 == "GATC" | ALTseq4 == "GATC", Dpn_site_ALT := TRUE]
new_seqs[Dpn_site_ALT == TRUE] # For these cases, will need to remove the SNPs in LD within the same original DpnII fragment.

fwrite(new_seqs, file = "gwas_and_nongwas_noneQTLs_otherEnd_LD.DpnII.REFALTseqs.txt", 
       sep = "\t", quote = F, row.names = F, col.names = T)


# Now continue with analysis only looking at those that do not intersect DpnII cut sites.
# We need to also remove SNPs in LD that were in the same original DpnII fragment as the one that creates a cut site.
new_seqs <- fread("gwas_and_nongwas_noneQTLs_otherEnd_LD.DpnII.REFALTseqs.txt")
# have to separate the new_seqs id column by underscore.
new_seqs_sep <- as.data.table(new_seqs %>% separate_rows(hg19Proxy_ID, sep = "_"))
length(unique(new_seqs_sep$hg19SNP_ID)) # 73,741
# put them into DpnII fragments.
dpnII <- fread("/rds/general/project/lms-spivakov-analysis/live/Design/Human_eQTL_CHiC_DpnII_hg38/hg38_dpnII.rmap")
names(dpnII) = c("Chr", "baitStart", "baitEnd", "baitID")
str(dpnII)

# get hg38 info from prev file and DpnII from rmap
# firstly the hg38pos for the SNP of interest
nongwas <- fread("~/HRJ_monocytes/AS_CHiC/BaseQTL/scripts/non_eqtls_non_gwas_otherEnd.snps")
gwas <- fread("~/HRJ_monocytes/AS_CHiC/BaseQTL/scripts/non_eqtls_gwas_otherEnd.snps")
# combine them for this.
both <- rbind(nongwas, gwas)

new_seqs_hg381 <- both[new_seqs_sep, on = c(hg19ID = "hg19SNP_ID")]
setnames(new_seqs_hg381, "hg38pos", "hg38SNP_pos")

# next the hg38pos for the proxy snp
# sadly have to go back to geno snps file
myseq <- fread("/rds/general/project/lms-spivakov-analysis/live/HRJ_monocytes/findmotifs/snp_fasta/ALL_genotyped_SNPs_REF_40bp_fasta_with_locations.txt")
myseq_small <- myseq[, .(SNP, hg38SNP_loc)]
setnames(myseq_small, c("SNP", "hg38SNP_loc"), c("hg19Proxy_ID", "hg38Proxy_pos"))
setkey(myseq_small, hg19Proxy_ID)

new_seqs_hg38 <- myseq_small[new_seqs_hg381, on = "hg19Proxy_ID", nomatch = NULL]

# intersect with DpnII once on the SNP and once on the proxy
new_seqs_hg38[, hg38SNP_pos2 := hg38SNP_pos]
new_seqs_hg38[, hg38Proxy_pos2 := hg38Proxy_pos]

dpnII[, Chr := as.numeric(Chr)]
setkey(dpnII, Chr, baitStart, baitEnd)
new_seqs_hg38_dpnII_1 <- foverlaps(new_seqs_hg38, dpnII, by.x = c("Chr", "hg38Proxy_pos", "hg38Proxy_pos2"), nomatch = NULL)
new_seqs_hg38_dpnII_1 # we keep everything
setnames(new_seqs_hg38_dpnII_1, c("baitStart", "baitEnd", "baitID"), c("Proxy_DpnStart", "Proxy_DpnEnd", "Proxy_baitID"))
new_seqs_hg38_dpnII <- foverlaps(new_seqs_hg38_dpnII_1, dpnII, by.x = c("Chr", "hg38SNP_pos", "hg38SNP_pos2"), nomatch = NULL)
new_seqs_hg38_dpnII # we keep everything
setnames(new_seqs_hg38_dpnII, c("baitStart", "baitEnd", "baitID"), c("SNP_DpnStart", "SNP_DpnEnd", "SNP_baitID"))


# Mark if the SNP is in the same DpnII site as the proxy.
new_seqs_hg38_dpnII[, sameFrag := SNP_baitID == Proxy_baitID]

### write the file of DpnII effects
setwd("~/HRJ_monocytes/AS_CHiC/BaseQTL/scripts/non_eqtls/dpnII_sites")
fwrite(new_seqs_hg38_dpnII, file = "gwas_and_nongwas_noneQTLs_otherEnd_LD.DpnII.marked.REFALTseqs.txt", 
       sep = "\t", quote = F, col.names = T, row.names = F)


# scenarios where we need to remove the snp:
to_remove <- new_seqs_hg38_dpnII[sameFrag == T & !is.na(Dpn_site_REF) | sameFrag == T & !is.na(Dpn_site_ALT) | SNP_baitID == Proxy_baitID - 1 & !is.na(Dpn_site_REF)]


new_seqs_hg38_dpnII_PURGED <- new_seqs_hg38_dpnII[!hg19ID %in% to_remove$hg19ID]
length(unique(new_seqs_hg38_dpnII_PURGED$hg19ID)) # 1029
length(unique(new_seqs_hg38_dpnII$hg19ID)) # 1058

new_seqs_hg38_dpnII_PURGED[, is_non_gwas := feature %like% "non_gwas"]
length(unique(new_seqs_hg38_dpnII_PURGED[is_non_gwas == T, hg19ID])) # 987
length(unique(new_seqs_hg38_dpnII_PURGED[is_non_gwas == F, hg19ID])) # 42

# Write this file and then carry on
fwrite(new_seqs_hg38_dpnII_PURGED, file = "gwas_and_nongwas_noneQTLs_otherEnd_LD.DpnII.marked.REFALTseqs.PURGED.txt", 
       sep = "\t", quote = F, row.names = F, col.names = T)

###################### from here 9/5/23
### Looking at baseQTL tesults
setwd("~/HRJ_monocytes/AS_CHiC/BaseQTL/scripts/non_eqtls/dpnII_sites/")
#setwd("~/Documents/analysis/BaseQTL_results/CHiC/round2_refBias/noneQTLs")
new_seqs_hg38_dpnII_PURGED <- fread("./gwas_and_nongwas_noneQTLs_otherEnd_LD.DpnII.marked.REFALTseqs.PURGED.txt")
new_seqs_hg38_dpnII <- fread("./gwas_and_nongwas_noneQTLs_otherEnd_LD.DpnII.marked.REFALTseqs.txt")
#res <- fread("~/HRJ_monocytes/BaseQTL/elena_results/CHi-C/round2_refBias/all_stan.summary.txt")
#res <- fread("~/Documents/analysis/BaseQTL_results/CHiC/round2_refBias/all_stan.summary.txt")
res <- fread("~/HRJ_monocytes/BaseQTL/elena_results/CHi-C/round2_refBias/all_stan.summary.txt")
sig <- res[Signif.99 == "yes" & Rhat < 1.01 & AI_estimate >= 0.4 | Signif.99 == "yes" & Rhat < 1.01 & model == "NB"]

check <- new_seqs_hg38_dpnII[feature %like% "3:180693504:C:T_non_gwas_chicago_score3_otherEnd"] # The SNP itself creates a DpnII site
new_seqs_hg38_dpnII_PURGED[feature %like% "3:180693504:C:T_non_gwas_chicago_score3_otherEnd"] # yep, this one is now removed


new_seqs_hg38_dpnII_PURGED[feature %like% "10:13260857:G:A_non_gwas_chicago_score3_otherEnd"]
new_seqs_hg38_dpnII_PURGED[feature %like% "10:13261773:C:T_non_gwas_chicago_score3_otherEnd"]
new_seqs_hg38_dpnII_PURGED[feature %like% "10:13261559:T:C_non_gwas_chicago_score3_otherEnd"]
# yep, these ones are ok

# generally it seems like if a locus has more than one snp, it is more likely to be real!

##### Overall compare the results for eQTL versus non eQTL (non gwas and gwas)
### have to first of all purge the results for dpnII effects in eQTL res and in gwas/nongwas res.
eqtl_purge <- fread("~/HRJ_monocytes/findmotifs/find_new_dpnII/snps_went_into_loopingQTL_analysis.DpnII.marked.REFALTseqs.PURGED.DISTAL.txt")
# This purge file was updated with proper filtering.
#eqtl_purge <- fread("~/Documents/analysis/BaseQTL_results/CHiC/round2_refBias/snps_went_into_loopingQTL_analysis.DpnII.marked.REFALTseqs.PURGED.DISTAL.txt")
# the SNPs were removed based on:
# [all_dpn_ALT %like% "TRUE" | all_dpn_REF %like% "TRUE" | dpn_REF_in_baitRight == TRUE]

## Filter on Rhat and AI estimates ... not sure if necessary but want to remove technical artefacts
res2 <- res[Rhat < 1.01 & AI_estimate >= 0.4 | Rhat < 1.01 & is.na(AI_estimate)]
res2[, feat := sub("-[^-]+$", "", int_id)] # remove SNP name after final "-". Because some genes have "-" in the name.
res2[, "hg19Proxy_ID" := tstrsplit(feat, split = "_", keep = 1)]

## remove DpnII and proximal effects: for either eQTLs or non eQTLs
resd <- res2[hg19Proxy_ID %in% new_seqs_hg38_dpnII_PURGED$hg19ID | hg19Proxy_ID %in% eqtl_purge$hg19Proxy_ID]

## Further feature processing
resd[, feature := sub("^[^_]*_", "", feat)]
# ^       beginning of string
# ^[^_]*  any character at the beginning except for "_" repeated zero or more times;
# ^[^_]*_ the pattern in point 2 above followed by an underscore

# binary type: eQTL or not?
resd[feature %like% "gwas", type.2 := "non_eQTL"]
resd[substr(feature, 1, 4) == "ENSG", type.2 := "eQTL"]
resd[is.na(type.2)]

# tertiary type: eQTL, gwas, or non gwas?
resd[feature %like% "non_gwas", type.3 := "non_gwas"]
resd[substr(feature, 1, 4) == "gwas", type.3 := "gwas"]
resd[substr(feature, 1, 4) == "ENSG", type.3 := "eQTL"]
resd[is.na(type.3)]

length(unique(resd[type.2 == "non_eQTL", hg19Proxy_ID])) # 55,591
length(unique(resd[type.2 == "eQTL", hg19Proxy_ID])) # 9,362
length(unique(resd[type.3 == "non_gwas", hg19Proxy_ID])) # 54,257
length(unique(resd[type.3 == "gwas", hg19Proxy_ID])) # 1,334
length(unique(resd[type.3 == "eQTL", hg19Proxy_ID])) # 9,362

### Normal approximation
resd[, inverse_norm_approx := 1-pnorm(abs(log_mean_aFC), mean = 0, sd = sd, lower.tail=FALSE) ]
fwrite(resd, file = "~/HRJ_monocytes/BaseQTL/findings_round2_dpnIICorrection/BaseQTL_CHiC_categories.DISTAL.txt", 
       sep = "\t", quote = F, row.names = F, col.names = T)
###

##### Determine if the read counts are similar or not.
## the counts col in resd only refers to total counts, so get the allele counts too.
hets <- fread("~/HRJ_monocytes/AS_CHiC/BaseQTL/for_elena/mono_34reps_CHiC_allTests_allele_counts_noMulti_hg38Strand.txt")
tots <- fread("~/HRJ_monocytes/AS_CHiC/BaseQTL/for_elena/mono_34reps_CHiC_allTests_total_counts_noMulti_hg38Strand.txt")
#hets <- fread("~/Documents/analysis/BaseQTL_input/CHiC/mono_34reps_CHiC_allTests_allele_counts_noMulti_hg38Strand.txt")
#tots <- fread("~/Documents/analysis/BaseQTL_input/CHiC/mono_34reps_CHiC_allTests_total_counts_noMulti_hg38Strand.txt")

hets[, hetCounts := rowSums(.SD), .SDcols = c(5:38)] # this is per allele
hets_bothAlleles <- as.data.table(hets %>% group_by(feature, SNP, hg38SNP_pos) %>% 
                                    summarise(totalAlleleCounts = sum(hetCounts))) # this is sum of both alleles

tots[, totCounts := rowSums(.SD), .SDcols = c(4:37)] # this is for both alleles
tots_allCounts <- tots[, .(feature, SNP, hg38SNP_pos, totCounts)]

# have to merge them on local machine, as openOnDemand cannot cope.
CountsCompare <- merge.data.table(tots_allCounts, hets_bothAlleles, by = c("feature", "SNP", "hg38SNP_pos"), all = T)
CountsCompare[is.na(CountsCompare)] <- 0

# Add the effects from BaseQTL
resd_small <- resd[, .(feat, inverse_norm_approx)]

### how do the different categories look?
# tertiary type: eQTL, gwas, or non gwas?
CountsCompare[, feat := feature]
CountsCompare[, int_id := paste(feature, SNP, sep = "-")]
CountsCompare[, feature := sub("^[^_]*_", "", feature)]
CountsCompare[feature %like% "non_gwas", type.3 := "non_gwas"]
CountsCompare[substr(feature, 1, 4) == "gwas", type.3 := "gwas"]
CountsCompare[substr(feature, 1, 4) == "ENSG", type.3 := "eQTL"]
CountsCompare[is.na(type.3)]
CountsCompare_resd <- CountsCompare[resd_small, on = c("feat"), nomatch = NULL]

# total counts
quantile(CountsCompare_resd[type.3 == "eQTL", totalAlleleCounts])
quantile(CountsCompare_resd[type.3 == "gwas", totalAlleleCounts])
quantile(CountsCompare_resd[type.3 == "non_gwas", totalAlleleCounts])

quantile(CountsCompare_resd[type.3 == "eQTL", totCounts])
quantile(CountsCompare_resd[type.3 == "gwas", totCounts])
quantile(CountsCompare_resd[type.3 == "non_gwas", totCounts])

CountsCompare_resd[, read_count := totalAlleleCounts + totCounts]
quantile(CountsCompare_resd[type.3 == "eQTL", read_count])
quantile(CountsCompare_resd[type.3 == "gwas", read_count])
quantile(CountsCompare_resd[type.3 == "non_gwas", read_count])

# gwas and non gwas are a bit lower.
# save this file

CountsCompare_resd_save <- CountsCompare_resd[, .(feat, SNP, hg38SNP_pos, type.3, read_count, inverse_norm_approx)]
fwrite(CountsCompare_resd_save, file = "./CountsCompare_BaseQTL_normAprrox_CHiC.DISTAL.txt", sep = "\t", 
       quote = F, row.names = F, col.names = T)


########################## Can be run from here on ##########################
############## UPDATE 12th JUNE 2023: look at non eQTL at score 5 ###########
#### How can we sample from the same distribution of readcounts?
library(data.table)
library(ggplot2)
library(dplyr)

setwd("/rds/general/project/lms-spivakov-analysis/live/HRJ_monocytes/BaseQTL/findings_round2_dpnIICorrection/")
#### First try a stratified random sampling method to normalise for read counts. Then, check median "normal approximation" from BaseQTL.
old <- fread("./CountsCompare_BaseQTL_normAprrox_CHiC.txt") # just to check
CountsCompare_resd <- fread("./CountsCompare_BaseQTL_normAprrox_CHiC.DISTAL.txt") # we have removed eQTL distal effects, but allowed for AI to be "NA"
length(unique(CountsCompare_resd[type.3 != "eQTL", SNP])) # 58,046

### Now, if the non eQTL did not have all its interactions at score5, remove it
score5 <- fread("~/HRJ_monocytes/AS_CHiC/BaseQTL/scripts/non_eqtls_score5_keep.txt")

keep_eQTL <- CountsCompare_resd[type.3 == "eQTL"]
keep_noneQTL <- CountsCompare_resd[feat %in% score5$feature]

CountsCompare_resd <- rbind(keep_eQTL, keep_noneQTL)

length(unique(keep_noneQTL$SNP))
length(unique(keep_eQTL$SNP))
# 11,582 non eQTL
# 8,582 eQTL

# 
fwrite(CountsCompare_resd, file = "./CountsCompare_BaseQTL_normAprrox_CHiC.DISTAL_noneQTLscore5.txt", sep = "\t", quote = F, row.names = F, col.names = T)

###### check for significant results

CountsCompare_resd <- fread("./CountsCompare_BaseQTL_normAprrox_CHiC.DISTAL_noneQTLscore5.txt") # we have removed eQTL distal effects, but allowed for AI to be "NA"
res <- fread("../elena_results/CHi-C/round2_refBias/all_stan.summary.txt")
noneq <- res[int_id %like% "gwas"]
noneq[, feature := tstrsplit(int_id, split = "-", keep = 1)]
noneq_sig <- noneq[Signif.99 == "yes" & AI_estimate >= 0.4 & Rhat < 1.01]

CountsCompare_resd[feat == "3:180693504:C:T_non_gwas_chicago_score3_otherEnd"]
CountsCompare_resd[feat == "10:13261773:C:T_non_gwas_chicago_score3_otherEnd"] # this one
CountsCompare_resd[feat == "10:13260857:G:A_non_gwas_chicago_score3_otherEnd"]
CountsCompare_resd[feat == "10:13261559:T:C_non_gwas_chicago_score3_otherEnd"] # this one

noneq_sig_score5 <- noneq_sig[feature %in% CountsCompare_resd$feat] # two loci
fwrite(noneq_sig_score5, file = "BaseQTL_CHiC_noneQTLs_signif_score5.txt", sep = "\t", quote = F, row.names = F, col.names = T)
noneq_sig_score5 # # They are in the same DpnII fragment and interacting with RNU6-6P. but it is a very short interaction! 1kb...


########################## ########################## ##########################
CountsCompare_resd <- fread("./CountsCompare_BaseQTL_normAprrox_CHiC.DISTAL_noneQTLscore5.txt") # we have removed eQTL distal effects, but allowed for AI to be "NA"
length(unique(CountsCompare_resd[type.3 != "eQTL", SNP])) # 11,582

eqtl <- CountsCompare_resd[type.3 == "eQTL"]
gwas <- CountsCompare_resd[type.3 == "gwas"]
nongwas <- CountsCompare_resd[type.3 == "non_gwas"]

max_count <- max(c(max(eqtl$read_count), max(gwas$read_count), max(nongwas$read_count)))


######## CODING WITH MIKHAIL
quantile(CountsCompare_resd[type.3 == "eQTL", read_count])
nrow(CountsCompare_resd[type.3!="eQTL" & read_count>112])
nrow(CountsCompare_resd[type.3=="eQTL" & read_count>112])

me <- vector("numeric", 10000)
for(i in 1:10000) {
  q4=sample(CountsCompare_resd[type.3!="eQTL" & read_count>112]$inverse_norm_approx, 200) #2187
  q3=sample(CountsCompare_resd[type.3!="eQTL" & read_count>45 & read_count<=112]$inverse_norm_approx, 200)
  q2=sample(CountsCompare_resd[type.3!="eQTL" & read_count>18 & read_count<=45]$inverse_norm_approx, 200)
  q1=sample(CountsCompare_resd[type.3!="eQTL" & read_count>5 & read_count<=18]$inverse_norm_approx, 200)
  expected=c(q1,q2,q3,q4)
  me[i] <- mean(expected)
}
length(me[me<mean(observed)])/length(me)
# therefore, we can reject the null hypothesis.

quantile(expected)
quantile(CountsCompare_resd[type.3 == "eQTL", inverse_norm_approx])
observed <- sample(CountsCompare_resd[type.3 == "eQTL", inverse_norm_approx], 8748)
# you can do that mulitple times but might not make that much difference. Then:
# see the difference at high inverse norm approx:
mean(expected)
mean(observed)

quantile(expected, 0.99)
quantile(observed, 0.99)
# fiure out the qqplot
pdf(file = "qq_CHiCBaseQTL_score5.pdf", width = 9, height = 4.5)
qqplot(expected, observed)
abline(a = 0, b = 1)
dev.off()

################## Use package EnvStats
library(EnvStats)
pdf(file = "tukey_qq_CHiCBaseQTL_score5.pdf", width = 9, height = 4.5)
qqPlot(observed, expected, plot.type = "T", add.line = T, 
       xlab = "Mean of percentiles for eQTLs and non-eQTLs", 
       ylab = "Non-eQTL percentiles minus eQTL percentiles", 
       main = "Tukey mean-difference QQ plot for eQTLs and non-eQTLs")
dev.off()

###################

## Can we just build a ECDF plot? For both expected and observed?
library(ggplot2)
library(RColorBrewer)
mycols <- 
dt <- data.table(`Non eQTL` = expected, eQTL = observed)
dtplot <- melt.data.table(dt, measure.vars = c("Non eQTL", "eQTL"), value.name = "effect", 
                          variable.name = "Type of SNP")
pdf(file = "ECDF_CHiCBaseQTL_score5.pdf", width = 6, height = 3)
ggplot(dtplot, aes(effect, color = `Type of SNP`)) + stat_ecdf(geom = "point", size = 0.2, alpha = 0.4, pad = FALSE) +
  xlab("Sampled BaseQTL inverse normal approximation") +
  ylab("Percentile") + 
  theme_bw() +
  scale_color_brewer(palette = "Paired") +
  guides(colour = guide_legend(override.aes = list(size=2, alpha = 1)))
dev.off()

# are they significanrly different?
wilcox.test(effect ~ `Type of SNP`, data = dtplot)
# yes, significantly different from each other.
# p = 4.531 e-12

################### older stuff below






qqplot(-log10(1-expected), -log10(1-observed))
lm(-log10(1-expected) ~ -log10(1-observed))
expected2 <- -log10(1-expected)
observed2 <- -log10(1-observed)
lm(observed2 ~ expected2)
# error because some valyes of expected are infinite
any(is.infinite(expected2)) # because log10(0) = -inf
any(is.infinite(expected))
#Error in lm.fit(x, y, offset = offset, singular.ok = singular.ok, ...) : 
#  NA/NaN/Inf in 'x'

# so we can exclude these
both <- data.table("expected" = expected, "observed" = observed)

both2 <- both[expected < 1 & observed < 1]
lm((-log10(1-both2$observed)) ~ (-log10(1-both2$expected)))
# ????????????????????????????????????
qqplot(-log10(1-both2$expected), -log10(1-both2$observed))
model <- lm(-log10(1-observed) ~ -log10(1-expected), data = both2)
qqplot(-log10(1-both2$expected), -log10(1-both2$observed))
abline(model)
# does not work.
# try the package QQperm with function qqplot


# Combine the data frames and use the ipw column as weights to create a weighted sample
sampled_data <- as.data.table(CountsCompare_resd %>%
                                sample_n(size = 300, replace = TRUE, weight = ipw) %>%
                                select(feat, type.3, read_count, inverse_norm_approx)) # remove ipw column from sampled data
setindex(sampled_data, NULL)
quantile(CountsCompare_resd[type.3 == "eQTL", read_count])
quantile(CountsCompare_resd[type.3 == "gwas", read_count])
quantile(CountsCompare_resd[type.3 == "non_gwas", read_count])
# they look more similar.
## try with a few iterations.
# Set the number of iterations
n_iter <- 100

# Initialize an empty list to store the sampled data frames
sampled_data_list <- list()

# Repeat the sampling process n_iter times
for (i in 1:n_iter) {
  
  # Code for stratified random sampling or IPW (depending on the approach)
  sampled_data <- CountsCompare_resd %>%
    sample_n(size = 300, replace = FALSE, weight = ipw) %>%
    select(feat, type.3, read_count, inverse_norm_approx) # remove ipw column from sampled data
  
  # Store the sampled data frame in the list
  sampled_data_list[[i]] <- sampled_data
}

# Combine the sampled data frames into a single data frame
sampled_data_all <- do.call(rbind, sampled_data_list)

## Now, set up a list to store the effect sizes
effect_sizes <- list()

# Loop over the sampled data frames and calculate the effect size for each category
dt <- data.table()
for (i in 1:n_iter) {
  
  # Extract the category and measurement variable columns from the sampled data frame
  data_i <- as.data.table(sampled_data_list[[i]])
  
  setindex(data_i, NULL)
  eqtl_res <- median(data_i[type.3 == "eQTL", inverse_norm_approx], na.rm = T)
  gwas_res <- median(data_i[type.3 == "gwas", inverse_norm_approx], na.rm = T)
  nongwas_res <- median(data_i[type.3 == "non_gwas", inverse_norm_approx], na.rm = T)
  dti <- data.table(eQTL = eqtl_res, gwas = gwas_res, non_gwas = nongwas_res)
  
  # Store the effect size in the effect_sizes list
  dt <- rbind(dt, dti)
}

to_plot <- melt.data.table(dt, measure.vars = c("eQTL", "gwas", "non_gwas"), variable.name = "Category", value.name = "median_inverse_normal_approx")

# Plot a dotplot of the effect sizes
library(ggplot2)
ggplot(to_plot, aes(x = Category, y = median_inverse_normal_approx, col = Category)) +
  geom_point() +
  xlab("Category") +
  ylab("median_inverse_normal_approx")

# Plot boxplots of the effect sizes
ggplot(to_plot, aes(x = Category, y = median_inverse_normal_approx, col = Category)) +
  geom_boxplot() +
  xlab("Category") +
  ylab("median_inverse_normal_approx")

### After sampling based on readCount, the normal approximations still look different between categories.

### Earlier stuff doing a basic comparison of effects (not comparing for read counts)
### Looking in either 2 categories (non eQTL or eQTL) or 3 categories, as above.
### Also splitting by test type: within indiviudals and between individuals.
resd <- fread("/rds/general/project/lms-spivakov-analysis/live/HRJ_monocytes/BaseQTL/findings_round2_dpnIICorrection/BaseQTL_CHiC_categories.txt")
ggplot(resd, aes(x = type.2, y = inverse_norm_approx)) + geom_boxplot()+
  scale_fill_brewer(palette = "Paired") + theme_light() + 
  ylab("Inverse normal approximation, BaseQTL CHi-C") +
  xlab("Type of SNP")
wilcox.test(inverse_norm_approx ~ type.2, data = resd)
# very significant.

# split by the type of model.

ggplot(resd, aes(x = type.2, y = inverse_norm_approx, fill = model)) + geom_boxplot() +
  scale_fill_brewer(palette = "Paired") + theme_light() +
  ylab("Inverse normal approximation, BaseQTL CHi-C") +
  xlab("Type of SNP")

# split by the type of non-eQTL: GWAS or non GWAS
ggplot(resd, aes(x = type.3, y = inverse_norm_approx, fill = type.2)) + geom_boxplot() +
  scale_fill_brewer(palette = "Paired") + theme_light() + 
  ylab("Inverse normal approximation, BaseQTL CHi-C") +
  xlab("Type of SNP")

kruskal.test(inverse_norm_approx ~ type.3, data = resd)
dunnTest(inverse_norm_approx ~ type.3, data = resd, method = "bonferroni")
# all significantly different from each other.







