library(data.table)

########### The basis for this script is ~/HRJ_monocytes/eqtls/normalise_for_leo/normalise_rnaseq_for_leo.R

########### Make a file of eQTL/proxies to eGene promoters.
##### Have DpnII fragments at the SNP and TSS intersected with the 5Kb baitmap (see also method in "~/HRJ_monocytes/AS_CHiC/scripts/mono_34reps_analysis/get_CHiC_bed_Regions_testControl.R")

# SNPs_to_test are from here: ~/HRJ_monocytes/CHiC/chicago/analysis_scripts/contacting_noncontacting_and_TADs_34reps.R
all_info <- fread("~/HRJ_monocytes/eqtls/snps_went_into_loopingQTL_analysis_with_gene_targets.txt")

### Check that this file has all the TSS??
### Yes, looks fine to me.

rmap <- fread("/rds/general/project/lms-spivakov-analysis/live/Design/Human_eQTL_CHiC_DpnII_hg38_binsize4000_maxL100K/hg38_dpnII.rmap")
names(rmap) = c("Chr", "hg38Dpn_start", "hg38Dpn_end", "hg38DpnID")
# Get the DpnID for the SNPs
setkey(rmap, Chr, hg38Dpn_start, hg38Dpn_end)
rmap[, Chr := as.numeric(Chr)]
all_info_dpn <- rmap[all_info, on = c("Chr", "hg38Dpn_start", "hg38Dpn_end")]
setnames(all_info_dpn, c("hg38Dpn_start", "hg38Dpn_end", "hg38DpnID"), c("SNP_Dpnstart", "SNP_Dpnend", "SNP_DpnID"))


### But, we are not getting as many TSS into this file as we did before (using get_CHiC_bed_Regions_testControl.R))
#length(unique(all_info_dpn$TSS))
#[1] 4183
#> length(unique(eQTL_Genes_targets_egene$cap_TSS_start))
#[1] 8789

### Need to find out why not. 
### In this script I use "all info" - is this only the "CAPTURED" TSS?
### Yes! First of all, I used the file of all TSS that we "TRIED" to capture. Will add those TSS.
original <- fread("~/eCHiC/design/final_design/all_captured_genes_finalhg38.txt") # Has info on mWindow or eWindow. This file actually includes all the TSS that we *tried* to capture.
tried_egenes <- unique(original[captured_gene_type == "in_eGene_window", .(hg38SNP_ID, cap_Gene, cap_ENSG_ID, cap_TSS_start)]) # keep "captrued gene type" to keep track of extra TSS
setnames(tried_egenes, c("cap_Gene", "cap_ENSG_ID"), c("Gene", "ENSG_ID"))
### Check if all these TSS are in "all_info_dpn", if so then we can just join based on SNP/Gene
all_info_dpn[!TSS %in% tried_egenes$cap_TSS_start] # almost all, but some are not - so would need to merge (probably due to the windows)

##### Join based on these cols:
mycols = c("hg38SNP_ID", "Gene", "ENSG_ID")

both1 <- tried_egenes[all_info_dpn, on = c(mycols), allow.cartesian = TRUE] # this is the intersection of the SNP/Gene info but KEEP the ones from all_info_dpn
# Separate the types of TSS into two tables and then rbind

first <- both1[, .(hg38SNP_ID, hg38Proxy_ID, Gene, ENSG_ID, Chr, SNP_Dpnstart, SNP_Dpnend, SNP_DpnID, hg38Proxy_pos, hg38Proxy_end, REF, ALT, Strand, hg19Proxy_ID, hrc_hg38Proxy_ID, hg19SNP_ID, hg38SNP_pos, rep_Gene, TSS)]
second0 <- both1[, .(hg38SNP_ID, hg38Proxy_ID, Gene, ENSG_ID, Chr, SNP_Dpnstart, SNP_Dpnend, SNP_DpnID, hg38Proxy_pos, hg38Proxy_end, REF, ALT, Strand, hg19Proxy_ID, hrc_hg38Proxy_ID, hg19SNP_ID, hg38SNP_pos, rep_Gene, cap_TSS_start)]
second <- second0[!is.na(cap_TSS_start)]
setnames(second, "cap_TSS_start", "TSS")

all_tss <- unique(rbind(first, second))
length(unique(all_tss$TSS)) # 8749!
length(unique(all_tss$hg38Proxy_ID)) # 14647
length(unique(all_tss$hg38SNP_ID)) # 1369 perfect


# Now intersect the TSS with 5K baitmap
rmap5k <- fread("/rds/general/project/lms-spivakov-analysis/live/Design/Human_eQTL_CHiC_DpnII_hg38_5K_bins/human_DpnII_5K.rmap")
names(rmap5k) = c("Chr", "start5k", "end5k", "index5k")
rmap5k[, Chr := as.numeric(Chr)]
setkey(rmap5k, Chr, start5k, end5k)

all_tss[, TSS0 := TSS]
all_info_dpn_TSS5K <- foverlaps(all_tss, rmap5k, by.x = c("Chr", "TSS0", "TSS"), nomatch = NULL)
all_info_dpn_TSS5K[is.na(index5k)] # None
setnames(all_info_dpn_TSS5K, c("start5k", "end5k", "index5k"), c("hg38TSS_5Kstart", "hg38TSS_5Kend", "hg38TSS_5KDpnID"))
all_info_dpn_TSS5K[, TSS0 := NULL]

# Check how far the TSS are from the ends of the 5K region
to_check <- copy(all_info_dpn_TSS5K)
to_check[, TSS_5Kdist1 := TSS - hg38TSS_5Kstart]
to_check[, TSS_5Kdist2 := hg38TSS_5Kend - TSS]
to_check[TSS_5Kdist1 < TSS_5Kdist2, minDist := TSS_5Kdist1]
to_check[TSS_5Kdist1 > TSS_5Kdist2, minDist := TSS_5Kdist2]
to_check[TSS_5Kdist1 == TSS_5Kdist2, minDist := TSS_5Kdist2]

to_plot <- unique(to_check[, .(Gene, TSS, minDist)])
hist(to_plot$minDist, breaks = 20)
range(to_plot$minDist) # From 0 to 4045
unique(to_check[Gene == "SLC25A29", c("Gene", "TSS", "hg38TSS_5KDpnID")]) 

length(unique(all_info_dpn_TSS5K$hg19Proxy_ID)) #14,647
length(unique(all_info_dpn_TSS5K$hg19SNP_ID)) #1,369


### To be honest, it looks ok. But we could add a buffer around TSS such that if the 500bp region around them is captured.
### i.e. extend by 250bp each way of the TSS.
all_tss[, TSS0 := TSS-250]
all_tss[, TSS1 := TSS+250]
all_info_dpn_TSS5K <- foverlaps(all_tss, rmap5k, by.x = c("Chr", "TSS0", "TSS1"), nomatch = NULL)
all_info_dpn_TSS5K[is.na(index5k)] # None
setnames(all_info_dpn_TSS5K, c("start5k", "end5k", "index5k"), c("hg38TSS_5Kstart", "hg38TSS_5Kend", "hg38TSS_5KDpnID"))

## Compare with before
to_check
all_info_dpn_TSS5K # We go from 122,125 lines to 132,604 (probably worth doing??)
# To be honest, I'm still not sure why we have fewer counts, because before we only had 27,889 "interactions" to test for the eGenes.
prev <- fread("~/HRJ_monocytes/AS_CHiC/WASP/input/eGene_loci_to_filter_Oct21.txt")
check <- unique(all_info_dpn_TSS5K[, .(Chr, hg38TSS_5Kstart, hg38TSS_5Kend, SNP_Dpnstart, SNP_Dpnend)])
### Now we have 29,606 lines!



# Check how far the TSS are from the ends of the 5K region
to_check <- copy(all_info_dpn_TSS5K)
to_check[, TSS_5Kdist1 := abs(TSS - hg38TSS_5Kstart)]
to_check[, TSS_5Kdist2 := abs(hg38TSS_5Kend - TSS)]
to_check[TSS_5Kdist1 < TSS_5Kdist2, minDist := TSS_5Kdist1]
to_check[TSS_5Kdist1 > TSS_5Kdist2, minDist := TSS_5Kdist2]
to_check[TSS_5Kdist1 == TSS_5Kdist2, minDist := TSS_5Kdist2]
to_plot <- unique(to_check[, .(Gene, TSS, minDist)])
hist(to_plot$minDist, breaks = 20)
range(to_plot$minDist) # Still from 1 to 4045 (but now we have included the adjacent frag when necessary)
unique(to_check[Gene == "SLC25A29", c("Gene", "TSS", "hg38TSS_5KDpnID")]) # we have added a frag

######### CONTINUE

short <- all_info_dpn_TSS5K[, .(Chr, hg38TSS_5Kstart, hg38TSS_5Kend, hg38TSS_5KDpnID, SNP_Dpnstart, SNP_Dpnend, SNP_DpnID, hg19Proxy_ID, ENSG_ID, Gene)]
short[, id := paste(SNP_DpnID, hg38TSS_5KDpnID, hg19Proxy_ID, ENSG_ID, Gene, sep = "_")]

# Get a number for each gene target

prox_gene <- unique(short[, .(hg19Proxy_ID, ENSG_ID, Gene)])
# Get a unique number for each gene when there is more than one gene in the locus.
prox_gene[, Gene_number := seq_len(.N), by = "hg19Proxy_ID"]
#prox_gene2 <- as.data.table(prox_gene %>% group_by(Gene) %>% mutate(max_number = max(Gene_number))) This doesnt work for e.g. 10:13159362:G:A; you would end up with three different genes in the same bam
more_than_one <- prox_gene[Gene_number != 1]
unique(prox_gene$Gene_number) # There are up to four
## Check that OSCAR etc exist?
prox_gene[Gene == "OSCAR"]
prox_gene[hg19Proxy_ID == "19:54549026:T:C"]
# Looks ok
prox_gene[Gene %like% "MS4"]
prox_gene[hg19Proxy_ID == "11:59964992:G:A"]
# yep
# Check: Have three: 10:13159362:G:Ahas 4 11:134140586:A:T, 11:134108059:G:T 
prox_gene[hg19Proxy_ID == "10:13159362:G:A"]
prox_gene[hg19Proxy_ID == "11:134108059:G:T"]
prox_gene[hg19Proxy_ID == "11:134140586:A:T"]
# Check in original file
mik <- fread("~/eCHiC/design/source/cd14_eqtls_no_freq.txt")
mik[snp_id == "19:54553697:T:C"]
mik[snp_id == "10:13160035:G:A"]
# Well it seems ok, we can double check again later
prox_gene[, Gene_number := paste0("Gene_", Gene_number)]

short_with_id <- short[prox_gene, on = c("hg19Proxy_ID", "Gene", "ENSG_ID")]

short_with_id[, c("hg19Proxy_ID", "ENSG_ID", "Gene") := NULL]
short_un <- unique(short_with_id)
short_un[, Chr := as.character(Chr)]

short_un_snps <- short_un[, .(Chr, SNP_Dpnstart, SNP_Dpnend, id)]
short_un_tss <- short_un[, .(Chr, hg38TSS_5Kstart, hg38TSS_5Kend, id)]

## Write our "bedpe"
bedpe_to_write <- short_un[, .(Chr, SNP_Dpnstart, SNP_Dpnend, Chr, hg38TSS_5Kstart, hg38TSS_5Kend, id, Gene_number)]
fwrite(bedpe_to_write, "~/HRJ_monocytes/AS_CHiC/BaseQTL/scripts/bedpe_for_chic_egenes_with_Gene_number.txt", 
       quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

################## Reading in the results from pair2pair - see ~/HRJ_monocytes/eqtls/normalise_for_leo/normalise_rnaseq_for_leo.R
