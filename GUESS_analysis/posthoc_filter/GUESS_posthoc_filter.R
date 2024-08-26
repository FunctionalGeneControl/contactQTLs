library(data.table)
library(tidyr)
library(dplyr)
library(stringr)

### Filter and explore the triplet results from Leo and Alex.
### Now running for the new results (dec22) after bug fixes
### Jan 24: Running the new results from Leo after including WASP filtering.

#setwd("~/Documents/analysis/Data/Leo")
setwd("/rds/general/project/lms-spivakov-analysis/live/HRJ_monocytes/leo_triplets")
#leo_mppi <- fread("./results/aug22/PP_stats_triplet_mPPI_SMALL.txt") # in this file, mPPI is filtered to >0.99
#leo <- fread("./results/aug22/PP_stats_triplet_SMALL.txt")
#leo <- fread("./results/dec22/PP_stats_triplet_SMALL_Dec22.txt")
#leo <- fread("./results/Rev/R2GUESS_PP_stats_triplet_mPPI_FDR_0.05_RS_FINAL.txt") # RS is the correct one. CC is comaparisons
# this is the final revised table sent by Leo in July 2024
leo <- fread("./results/Rev/Revision_01072024/PP_stats_triplet_mPPI_FDR_0.05_RS_Revision_01072024.txt")
trip <- fread("./input/aug22/triplet_list_100622.txt")
#chic <- fread("./input/aug22/eGenes_CHiC_counts_rlogtransform.txt")
chic <- fread("~/HRJ_monocytes/eqtls/normalise_for_leo/submission_1stDecember2023/eGenes_CHiC_counts_rlogtransform.txt")


####################  FILTERING   ####################  

######## DpnII cut sites.

#### We need to remove those features (triplets) where SNPs could alter a DpnII cut site.
## IF THE SNP IS IN A KNOWN CUT SITE: We remove the features that involve fragments to the left and right of the cut site. Do not keep SNP sets if other chic features still remain.
## IF THE SNP CREATES A CUT SITE: We remove the feature involving the original fragment. Do not keep SNP sets if other chic features still remain.
## We have to remove all the SNP sets because that fragment could have been dragged along in the ditag of an affected fragment and counted by HiCUP combinations.

######## Remove features where the SNPs fell within - or altered - a DpnII cut site.
set_seqs <- fread("~/HRJ_monocytes/findmotifs/find_new_dpnII/LeoSNPsets_DpnII_REFALT.txt")
#set_seqs <- fread("./filtering/LeoSNPsets_DpnII_REFALT.txt")
# These were generated using the script in findmotifs: get_dpnII_new_sites_onlyLeoSig.R (runs from the command line in a couple of hrs but would be quicker if optimised)
# In this table, have mutated the sequence based on SNP locations and indicated if there is a DpnII cut site intersecting the SNP.
set_seqs_dpnII <- set_seqs[Dpn_site_REF == TRUE | Dpn_site_ALT == TRUE] # 42% - but, we only care if they also intersect a feature.
# need to get SNP hg38 locations.
snps_hg38 <- fread("~/HRJ_monocytes/genotyping/imputed/hg38_filtered/all.rsq3.maf05.hg38.dose.finalRenamed.bim")
#snps_hg38 <- fread("./filtering/all.rsq3.maf05.hg38.dose.finalRenamed.bim")
names(snps_hg38) <- c("Chr", "hg19ID", "what", "hg38loc", "REF", "ALT")
snps_hg38[, c("what", "REF", "ALT") := NULL]
set_seqs2 <- as.data.table(set_seqs %>% separate_rows(id, sep = "_")) # Remember that some of these will have had more than one SNP in the sequence. Need to come back to this later on.
set_seqs_hg38 <- snps_hg38[set_seqs2, on = c(hg19ID = "id"), nomatch = NULL]
set_seqs_hg38[, hg38loc2 := hg38loc]
fwrite(set_seqs_hg38, file = "~/spivakov/HRJ_monocytes/findmotifs/find_new_dpnII/LeoSNPsets_DpnII_REFALT_sepWhg38.txt", 
       sep = "\t", quote = F, row.names = F, col.names = T)
##################

# Get the chic locations (these are in hg38)
chic_locs <- chic[, c("feature", "hg38_location")]
chic_locs[, feature := str_replace_all(feature, "-", ".")]
leo[, chic_feat := tstrsplit(`Triplet name`, split = "/", keep = 1)]
#leo[, chic_feat := substr(chic_feat, start = 2, stop = 1000)] # only if it starts with "X"
leo[, DpnID := tstrsplit(chic_feat, split = "_", keep = 1)]

leo_locs <- chic_locs[leo, on = c(feature = "chic_feat")]
leo_locs[is.na(hg38_location)] # Good.

# Intersect the chic features with the snp locations and indicate if they are in the cut site.
leo_locs[, c("Chr", "Dpn_start", "Dpn_end") := tstrsplit(hg38_location, split = ":", type.convert = TRUE)]
setkey(leo_locs, Chr, Dpn_start, Dpn_end)

# set_seqs_hg38 was generated earlier. If not, read in again.
set_seqs_hg38 <- fread("~/spivakov/HRJ_monocytes/findmotifs/find_new_dpnII/LeoSNPsets_DpnII_REFALT_sepWhg38.txt")
setkey(set_seqs_hg38, Chr, hg38loc, hg38loc2)
leo_wsnps <- foverlaps(set_seqs_hg38, leo_locs, by.x = c("Chr", "hg38loc", "hg38loc2"), nomatch = NULL)
# Set the triplet ID to be the same - so that we can filter on it later.
leo_wsnps_match <- leo_wsnps[`Triplet ID` == set]

######### What do we do?
######### if any of the snps make a site, maybe we have to remove the whole set. Cannot guarantee that the frag was not brought along.
######### So what if we just remove those that are in the site and intersecting a CHiC fragment - how many do we lose?

leo_wsnps_match_dpnII <- leo_wsnps_match[Dpn_site_REF == TRUE | Dpn_site_ALT == TRUE] # 221 sets out of 8974 - Dec results
# So we can probably live with losing these?
# Need to think about: in cases of REF DpnII sites we should also check for SNPs in adjacent fragments.
# Therefore from the chic locations we also need -1 and ask if there is a Dpn site REF. (The SNPs are mapping to the RHS of the cut site).
rmap <- fread("~/spivakov/Design/Human_eQTL_CHiC_DpnII_hg38/hg38_dpnII.rmap")
#rmap <- fread("~/Documents/analysis/Data/Design/hg38_dpnII.rmap")
names(rmap) = c("Chr", "DpnStartLeft", "DpnEndLeft", "DpnIDLeft")
chic_locs_extended1 <- copy(chic_locs)
chic_locs_extended1[, DpnID := tstrsplit(feature, split = "_", keep = 1, type.convert = T)]
chic_locs_extended1[, DpnIDLeft := DpnID - 1]
setkey(rmap, DpnIDLeft)
chic_locs_extended <- rmap[chic_locs_extended1, on = "DpnIDLeft"]

## intersect again with triplets.
leo_locs_extended <- chic_locs_extended[leo, on = c(feature = "chic_feat")]
leo_locs_extended[is.na(hg38_location)] # Good.

# Intersect the chic features with the snp locations for 1) the fragment and 2) the fragment -1 and indicate if they are in the cut site.
leo_locs_extended[, c("Chr", "Dpn_start", "Dpn_end") := tstrsplit(hg38_location, split = ":", type.convert = TRUE)]
setkey(set_seqs_hg38, Chr, hg38loc, hg38loc2)

leo_wsnps1 <- foverlaps(leo_locs_extended, set_seqs_hg38, by.x = c("Chr", "Dpn_start", "Dpn_end"))
leo_wsnps1_match <- leo_wsnps1[`Triplet ID` == set] # any here with Dpn_site_REF = true or Dpn_site_ALT = true should be kicked out. 
leo_wsnps1_match_dpn <- leo_wsnps1_match[Dpn_site_REF == TRUE | Dpn_site_ALT == TRUE]
leo_wsnps2 <- foverlaps(leo_locs_extended, set_seqs_hg38, by.x = c("Chr", "DpnStartLeft", "DpnEndLeft"))
leo_wsnps2_match <- leo_wsnps2[`Triplet ID` == set] # any here with Dpn_site_REF = true should be kicked out.
leo_wsnps2_match_dpn <- leo_wsnps2_match[Dpn_site_REF == TRUE]

sets_purge1 <- unique(leo_wsnps1_match_dpn[, .(set)])
sets_purge2 <- unique(leo_wsnps2_match_dpn[, .(set)])
sets_purge <- rbind(sets_purge1, sets_purge2)
#fwrite(sets_purge, file = "./results/aug22/filtered/sets_with_DpnII_sites_to_purge.txt", sep = "\t", 
#       quote = F, row.names = F, col.names = T)

dir.create("./results/Rev/Revision_01072024/filtered")
fwrite(sets_purge, file = "./results/Rev/Revision_01072024/filtered/sets_with_DpnII_sites_to_purge.txt", sep = "\t", 
       quote = F, row.names = F, col.names = T)


############################# save files
leo_dpnFilt <- leo[!`Triplet ID` %in% sets_purge$set]
#fwrite(leo_dpnFilt, file = "./results/aug22/filtered/PP_stats_triplet_SMALL_DpnFilt.txt", sep = "\t", quote = F, row.names = F, col.names = T)
#fwrite(leo_dpnFilt, file = "./results/dec22/filtered/PP_stats_triplet_SMALL_Dec22_DpnFilt.txt", sep = "\t", quote = F, row.names = F, col.names = T)
fwrite(leo_dpnFilt, file = "./results/Rev/Revision_01072024/filtered/R2GUESS_PP_stats_triplet_mPPI_FDR_0.05_RS_FINAL_DpnFilt.txt", sep = "\t", quote = F, row.names = F, col.names = T)

#print(leo_dpnFilt)
range(leo_dpnFilt$`Post-proc. FDR mPPI`)
length(unique(leo_dpnFilt$`Post-proc. FDR mPPI SNP`))
length(unique(leo_dpnFilt$chic_feat))
length(unique(leo_dpnFilt$`Triplet ID`))
leo_dpnFilt[, gene := tstrsplit(`Triplet name`, split = "/", keep = 3)]
length(unique(leo_dpnFilt$gene)) # 419 genes
print(leo)
#########################################################

