library(data.table)

########### The basis for this script is ~/HRJ_monocytes/eqtls/normalise_for_leo/normalise_rnaseq_for_leo.R

########### Make a file of eQTL/proxies to eGene promoters.
##### Have DpnII fragments at the SNP and TSS intersected with the 5Kb baitmap (see also method in "~/HRJ_monocytes/AS_CHiC/scripts/mono_34reps_analysis/get_CHiC_bed_Regions_testControl.R")

# SNPs_to_test are from here: ~/HRJ_monocytes/CHiC/chicago/analysis_scripts/contacting_noncontacting_and_TADs_34reps.R
all_info <- fread("~/HRJ_monocytes/eqtls/snps_went_into_loopingQTL_analysis_with_gene_targets.txt")
rmap <- fread("~/spivakov/Design/Human_eQTL_CHiC_DpnII_hg38_binsize4000_maxL100K/hg38_dpnII.rmap")
names(rmap) = c("Chr", "hg38Dpn_start", "hg38Dpn_end", "hg38DpnID")
# Get the DpnID for the SNPs
setkey(rmap, Chr, hg38Dpn_start, hg38Dpn_end)
rmap[, Chr := as.numeric(Chr)]
all_info_dpn <- rmap[all_info, on = c("Chr", "hg38Dpn_start", "hg38Dpn_end")]
setnames(all_info_dpn, c("hg38Dpn_start", "hg38Dpn_end", "hg38DpnID"), c("SNP_Dpnstart", "SNP_Dpnend", "SNP_DpnID"))

# Now intersect the TSS with 5K baitmap
rmap5k <- fread("~/spivakov/Design/Human_eQTL_CHiC_DpnII_hg38_5K_bins/human_DpnII_5K.rmap")
names(rmap5k) = c("Chr", "start5k", "end5k", "index5k")
rmap5k[, Chr := as.numeric(Chr)]
setkey(rmap5k, Chr, start5k, end5k)
all_info_dpn[, TSS0 := TSS]
all_info_dpn_TSS5K <- foverlaps(all_info_dpn, rmap5k, by.x = c("Chr", "TSS0", "TSS"), nomatch = NULL)
all_info_dpn_TSS5K[is.na(index5k)] # None
setnames(all_info_dpn_TSS5K, c("start5k", "end5k", "index5k"), c("hg38TSS_5Kstart", "hg38TSS_5Kend", "hg38TSS_5KDpnID"))
all_info_dpn_TSS5K[, TSS0 := NULL]

short <- all_info_dpn_TSS5K[, .(Chr, hg38TSS_5Kstart, hg38TSS_5Kend, hg38TSS_5KDpnID, SNP_Dpnstart, SNP_Dpnend, SNP_DpnID, hg19Proxy_ID, ENSG_ID, Gene)]
short[, id := paste(SNP_DpnID, hg38TSS_5KDpnID, hg19Proxy_ID, ENSG_ID, Gene, sep = "_")]

# Get a number for each gene target

prox_gene <- unique(short[, .(hg19Proxy_ID, ENSG_ID, Gene)])
prox_gene[, Gene_number := seq_len(.N), by = hg19Proxy_ID]
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
# Looks OK
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
