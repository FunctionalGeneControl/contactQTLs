library(DESeq2)
library(data.table)
#library(readxl)
library(binhf)

## UPDATE 1st December 2023
## Now I am adding results from WASP correction.

## UPDATE 22nd April 2022
## Now I am using the CHiC total counts per gene target that I have already generated for Elena.
## I am also going to update ATAC with the latest counts.
## Just using the counts files themselves and normalising.

setwd("~/HRJ_monocytes/eqtls/normalise_for_leo")
## Had to run this using conda environment DEseq2
## Read in the count data. Rather than from Novogene, this was now analysed in the WASP pipeline and then using phASER/FeatureCounts.
## We don't use the phASER results here.
genesdt <- fread("~/HRJ_monocytes/eqtls/WASP/mono_34reps/total_counts/featureCounts_mono34reps_GTF.txt", sep = "\t", header = TRUE)
## We use the same novogene file for genes as before, so that a proper comparison can be made. This is ensembl 94.
gene_locs <- fread("~/novogene/X204SC20081070-Z01-F001/X204SC20081070-Z01-F001_Homo_sapiens_result/3.Quantification/Count/readcount_genename.csv")
gene_locs_format <- gene_locs[, c("gene_id", "gene_chr", "gene_start", "gene_end", "gene_name")]

genesdt_loc <- genesdt[gene_locs_format, on = "gene_id", nomatch = NULL]

## Restrict to the genes we care about
cd14 <- fread("~/eCHiC/design/source/cd14_eqtls_no_freq.txt")
genes_filt1 <- genesdt_loc[gene_id %in% cd14$ensembl_gid]

# Also should bring in the gene name, incase they match rather than the ENSG ID
genes_filt2 <- genesdt_loc[gene_name %in% cd14$hgnc]

# combine
genes_filt <- unique(rbind(genes_filt1, genes_filt2)) # 4329 (previously 4762, probably due to mismatch with Novogene)

# Check those that are brought in by gene name alone
genes_filt[!gene_id %in% cd14$ensembl_gid]
tester <- genesdt_loc[cd14, on = c(gene_id = "ensembl_gid"), nomatch = NULL]
tester[gene_name != hgnc]
# Make sure no duplicates on ENSG_ID
genes_filt[duplicated(gene_id)]
# Continue

#hist(genes_filt$S025NM, breaks = 100)
genes_filt[, gene_loc := paste(gene_chr, gene_start, gene_end, sep = ":")]
genes_filt[, id := paste(gene_id, gene_name, gene_loc, sep = "_")]
genes_filt[, c("gene_id", "gene_name", "gene_chr", "gene_start", "gene_end", "gene_loc") := NULL]
setkey(genes_filt, id)
genes <- as.matrix(genes_filt, rownames = TRUE)
## Make the samples table. 
samples <- data.frame(row.names = colnames(genes), group = rep("1", 34), ID = colnames(genes)) # I've set group to "1" because there is no condition
## Make the DEseq object
dds <- DESeqDataSetFromMatrix(countData = genes, colData = samples, design = ~ 1) # have to use design ~1 because there is only one treatment.

## Collapsing technical replicates - this was previosuly done for M31B S026T6 and M19 S02699 in BaseQTL. Here, only M31B S026T6 was processed through WASP.
# M31B S026T6 and M19 S02699 as shown in ID col...
#dds_col <- collapseReplicates(dds, dds$ID)

## Remove sample S0266F_M16 - this was already done.
#test <- as.data.table(colData(dds))
#test[16, ]
#dds_col <- dds_col[, -16]

rna_un_to_write <- data.table(assay(dds), keep.rownames = TRUE)
names(rna_un_to_write)[1] = "feature"

##### Save the non-transformed counts
write.table(rna_un_to_write, file = "./WASP_temp/genes_interest_untransformed_counts.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)


##### Transform
dds_col_rlog <- rlog(dds)
rna_tr_to_write <- data.table(assay(dds_col_rlog), keep.rownames = TRUE)
names(rna_tr_to_write)[1] = "feature"
write.table(rna_tr_to_write, file = "./WASP_temp/genes_interest_rlog_counts.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)


##### Checking
rlog_counts <- fread("./WASP_temp/genes_interest_rlog_counts.txt")
hist(rlog_counts$S025NM, breaks = 100)
hist(rna_un_to_write$S025NM, breaks = 100)


#### Note I also tried vst function, but the data did not look as normal.
#### Now to do: normalise the counts from capture HiC and ATAC.
#### This time we are not adding SNPs.

### 1. ATAC.
atac <- fread("~/HRJ_monocytes/AS_ATAC/WASP/mono_34reps_noMultiMapped/allele_counts/mono_34reps_ATAC_total_counts_WASP_noMulti_hg38Strand.txt")
hist(atac$`S025NM-01`)

# Filter to peaks NEARBY our SNPs of interest!! Will remove SNPs as did before, but not add back in afterwards.
all_info <- fread("~/HRJ_monocytes/eqtls/snps_went_into_loopingQTL_analysis_with_gene_targets.txt")
atac_interest <- atac[SNP %in% all_info$hg19Proxy_ID]
write.table(atac_interest, file = "./WASP_temp/ATAC_peaks_withSNPs_counts_Dec23.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)


# Make into a matrix
atac_interest[, c("hg38SNP_pos", "SNP") := NULL] # Now remove snps, later intersect with the SNPs of interest based on their hg38 positions
atac_feat <- unique(atac_interest)
hist(atac_feat$`S025NM-01`, breaks = 100)
write.table(atac_feat, file = "./WASP_temp/ATAC_peaks_counts_Dec23.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

setkey(atac_feat, feature)
at_mat <- as.matrix(atac_feat, rownames = TRUE)
## Make the samples table. Set M19 ID to the same as M31B, these are technical reps
samples <- data.frame(row.names = colnames(at_mat), group = rep("1", 34), ID = colnames(at_mat)) # I've set group to "1" because there is no condition
## Make the DEseq object 
dds <- DESeqDataSetFromMatrix(countData = at_mat, colData = samples, design = ~ 1) # have to use design ~1 because there is only one treatment.


##### Transform
dds_rlog <- rlog(dds)
atac_tr_to_write <- data.table(assay(dds_rlog), keep.rownames = TRUE)
hist(atac_tr_to_write$`S025NM-01`)

names(atac_tr_to_write)[1] = "feature"
write.table(atac_tr_to_write, file = "./WASP_temp/ATAC_peaks_rlog_counts_Dec23.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

#### Now check how it looks.
atac_rlog <- fread("./WASP_temp/ATAC_peaks_rlog_counts_Dec23.txt")
hist(atac_rlog$`S025NM-01`)

# Also tried with VST (ignore)
#atac_vst <- fread("./ATAC_peaks_vst_counts.txt")

### 2. CHiC

### Need to remove snps and define the counts by Dpn bait to Gene promoters, because a bait can have more than one SNP in it (and more than one gene target!)
chic <- fread("~/HRJ_monocytes/AS_CHiC/WASP/mono_34reps_allReads_using_eGene_filters/allele_counts/mono_34reps_CHiC_total_counts_from_bedpe.txt")
chic[, c("SNP") := NULL]

# Add in DpnII frags before writing.
### Add in DpnII frags
dpn <- fread("~/spivakov/Design/Human_eQTL_CHiC_DpnII_hg38/hg38_dpnII.rmap")
names(dpn) = c("Chrom", "start", "stop", "DpnID")
chic[, chr := tstrsplit(feature, split = ":", keep = 1)]
chic[, hg38SNP_pos0 := hg38SNP_pos]
setkey(dpn, Chrom, start, stop)
chic_dpn <- foverlaps(chic, dpn, by.x = c("chr", "hg38SNP_pos0", "hg38SNP_pos"), nomatch = NULL)
to_normalise <- copy(chic_dpn)
to_normalise[, c("ENSG_ID", "gene") := tstrsplit(feature, split = "_", keep = c(2, 3))] # keeps the "gwas" and "non gwas" ones correctly as well.
to_normalise[, Dpn_loc := paste(chr, start, stop, sep = ":")]
to_normalise[, c("hg38SNP_pos0", "hg38SNP_pos", "chr", "start", "stop", "feature") := NULL]
to_normalise[, feature := paste(DpnID, Dpn_loc, ENSG_ID, gene, sep = "_")]
to_normalise[, c("ENSG_ID", "gene", "DpnID", "Dpn_loc") := NULL]

to_normalise_unique <- unique(to_normalise)
to_normalise_unique[duplicated(feature)] # none

hist(to_normalise_unique$`S025NM-01`, breaks = 100)
qqnorm(to_normalise_unique$`S025NM-01`)
fwrite(to_normalise_unique, file = "./WASP_temp/chic_eGene_counts_untransformed_Dec23.txt", 
       sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

setkey(to_normalise_unique, feature)
chic_mat <- as.matrix(to_normalise_unique, rownames = TRUE)

## Make the samples table. 
samples <- data.frame(row.names = colnames(chic_mat), group = rep("1", 34), ID = colnames(chic_mat)) # I've set group to "1" because there is no condition
## Make the DEseq object
dds <- DESeqDataSetFromMatrix(countData = chic_mat, colData = samples, design = ~ 1) # have to use design ~1 because there is only one treatment.

##### Transform
dds_rlog <- rlog(dds)

##### Have a look
dds_rlogdt <- as.data.table(assay(dds_rlog))
hist(dds_rlogdt$`S025NM-01`)
hist(dds_rlogdt$`S025QG-02`)
#qqnorm(dds_rlogdt$`S025NM-01`)


dds_rlogdt <- as.data.table(assay(dds_rlog), keep.rownames = TRUE)


names(dds_rlogdt)[1] = "feature"

fwrite(dds_rlogdt, file = "./WASP_temp/chic_eGene_counts_rlog_Dec23.txt", 
       sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

###############################################################################

######### Make the files with just feautre names, and make the intersections.

all_info # This has snp-gene combos with hg38 pos
rna_un <- fread("./WASP_temp/genes_interest_untransformed_counts.txt")
rna_tr <- fread("./WASP_temp/genes_interest_rlog_counts.txt")
ata_un <- fread("./WASP_temp/ATAC_peaks_counts_Dec23.txt")
ata_tr <- fread("./WASP_temp/ATAC_peaks_rlog_counts_Dec23.txt")
chic_un <- fread("./WASP_temp/chic_eGene_counts_untransformed_Dec23.txt")
chic_tr <- fread("./WASP_temp/chic_eGene_counts_rlog_Dec23.txt")

dir.create("submission_1stDecember2023")
setwd("submission_1stDecember2023")

########## RNASEQ (I have normalised before adding SNPs)
# Fix sample names in rna seq - add correct numbers
mynames <- c(colnames(chic_un[, c(35, 1:34)])) 
names(rna_un) = mynames

rna_un[, c("ensg", "gene_name", "hg38_location") := tstrsplit(feature, split = "_")]
rna_un[, feature := paste(ensg, gene_name, sep = "_")]
rna_un[, c("ensg", "gene_name") := NULL]
setcolorder(rna_un, c(1, 36, 2:35))
rna_un2 <- rna_un[order(hg38_location)]

fwrite(rna_un2, file = "eGenes_RNA_counts_untransformed.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

names(rna_tr) = mynames

rna_tr[, c("ensg", "gene_name", "hg38_location") := tstrsplit(feature, split = "_")]
rna_tr[, feature := paste(ensg, gene_name, sep = "_")]
rna_tr[, c("ensg", "gene_name") := NULL]
setcolorder(rna_tr, c(1, 36, 2:35))
rna_tr2 <- rna_tr[order(hg38_location)]

fwrite(rna_tr2, file = "eGenes_RNA_counts_rlogtransform.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)


########## ATACSEQ (I have normalised before adding SNPs, because we are counting peaks. 
## Not run again - computationally heavy
## get peak names
#peaks <- fread("~/ATAC-seq/analyses/HMMRATAC/consensus/output/ATAC_merged_34reps_HMMRATAC_withBlacklist_peaks.bed")
#peaks[, chr := tstrsplit(V1, split = "r", keep = 2)]
#peaks[, feature1 := paste(chr, V2, sep = ":")]
#peaks[, feature2 := paste(feature1, V3, sep = "-")]
#peaks[, feature1 := NULL]
#setnames(peaks, "feature2", "feature")

## Need to make peak names unique
#temp <- peaks[V4 == "HighCoveragePeak_0"]
#temp[, rn := row.names(temp)]
#temp[, peakname := tstrsplit(V4, split = "_", keep = 1)]
#temp[, V4 := paste(peakname, rn, sep = "_")]
#temp[, c("rn", "peakname") := NULL]

#nodup <- peaks[V4 != "HighCoveragePeak_0"]
#nodup[duplicated(V4)] # good

#allpeaks <- rbind(temp, nodup)
#allpeaks_to_write <- allpeaks[order(feature)]

#allpeaks_to_write[, c("chr", "feature") := NULL]
#fwrite(allpeaks_to_write, file = "~/HRJ_monocytes/ATAC-seq/analyses/HMMRATAC/consensus/output/ATAC_merged_34reps_HMMRATAC_withBlacklist_peaks_uniqueNames.bed", 
#       sep = "\t", row.names = F, col.names = F, quote = F)

allpeaks_to_write <- fread("~/HRJ_monocytes/ATAC-seq/analyses/HMMRATAC/consensus/output/ATAC_merged_34reps_HMMRATAC_withBlacklist_peaks_uniqueNames.bed")

allpeaks_to_write[, chr_start := paste(V1, V2, sep = ":")]
allpeaks_to_write[, chr_start_end := paste(chr_start, V3, sep = "-")]
allpeaks_to_write[, c("ch", "hg38_location") := tstrsplit(chr_start_end, split = "r")]
allpeaks <- allpeaks_to_write[, .(V4, hg38_location)]
names(allpeaks)[1] = "feature"

# Add peak names to the reads
names(ata_un)[1] = "hg38_location"
ata_un_peaks <- ata_un[allpeaks, on = "hg38_location", nomatch = NULL]

ata_un_peaks_to_write <- ata_un_peaks[order(hg38_location)]
setcolorder(ata_un_peaks_to_write, c(36, 1:35))

fwrite(ata_un_peaks_to_write, file = "peaks_ATAC_counts_untransformed.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# Normalised 
names(ata_tr)[1] = "hg38_location"
ata_tr_peaks <- ata_tr[allpeaks, on = "hg38_location", nomatch = NULL]

ata_tr_peaks_to_write <- ata_tr_peaks[order(hg38_location)]
setcolorder(ata_tr_peaks_to_write, c(36, 1:35))

fwrite(ata_tr_peaks_to_write, file = "peaks_ATAC_counts_rlogtransform.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

########## CHiC (I have normalised before adding SNPs i.e. normalised counts for each snp-containing DpnII fragment to 5K-bin-containing eGene promoters)
#### Before I had quite a complicated "feature", like the other traits. But now I have these "controls". Need to decide if we actually want them for this analysis?
#### Doesn't make sense, since they are non eqtls and we are looking for important predictors. And we would need to include them in all traits.
#### Removing from higher up

chic_un[, c("DpnID", "Dpnloc", "ENSG_ID", "Gene") := tstrsplit(feature, split = "_")]
chic_un[, feature := paste(DpnID, ENSG_ID, Gene, sep = "_")]
chic_un[, c("DpnID", "ENSG_ID", "Gene") := NULL]
setnames(chic_un, "Dpnloc", "hg38_location")


setcolorder(chic_un, c(35:36, 1:34))
chic_un_to_write <- chic_un[order(hg38_location)]

fwrite(chic_un_to_write, file = "eGenes_CHiC_counts_untransformed.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

chic_tr[, c("DpnID", "Dpnloc", "ENSG_ID", "Gene") := tstrsplit(feature, split = "_")]
chic_tr[, feature := paste(DpnID, ENSG_ID, Gene, sep = "_")]
chic_tr[, c("DpnID", "ENSG_ID", "Gene") := NULL]
setnames(chic_tr, "Dpnloc", "hg38_location")
setcolorder(chic_tr, c(1, 36, 2:35))

fwrite(chic_tr, file = "eGenes_CHiC_counts_rlogtransform.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# Done!

################### Now just getting snp locations. gzip this file and upload.
#snps <- fread("~/HRJ_monocytes/genotyping/imputed/hg38_filtered/all.rsq3.maf05.hg38.dose.finalRenamed.newIDs.sites.vcf")
#snps
################### 

################### UPDATE 10th June 2022
################### Now Providing "triplets" i.e. each combination of features between the three traits that we want to test.
################### ATAC to CHi-C bait = 5kb from each end
################### CHi-C bait to RNA-seq = "other end" ID (gene to gene)
################### 1st December 2023 I re-ran this for the WASP data.

### Get the feature names for each trait.

setwd("~/HRJ_monocytes/eqtls/normalise_for_leo")

chic <- fread("./submission_1stDecember2023/eGenes_CHiC_counts_untransformed.txt")
atac <- fread("./submission_1stDecember2023/peaks_ATAC_counts_untransformed.txt")
rna <- fread("./submission_1stDecember2023/eGenes_RNA_counts_untransformed.txt")

# I have been inconsistent with the way that hg38 locations are field separated, so deal with them separately.
chic[, c("Chr", "start", "end") := tstrsplit(hg38_location, split = ":", type.convert = TRUE)]
chic[, start_minus5 := start - 5000]
chic[, end_plus5 := end + 5000]

atac[, c("Chr", "loc") := tstrsplit(hg38_location, split = ":", type.convert = TRUE)]
atac[, c("start", "end") := tstrsplit(loc, split = "-", type.convert = TRUE)]

# Note that ATAC counts were already filtered to be within 5kb of our eQTLs of interest.
# Here we just assign them to the correct DpnII fragments. Since the SNPs were within these fragments, we can go based on the location of the DpnII fragment.
atac_feat <- unique(atac[, .(feature, Chr, start, end)])
names(atac_feat)[1] = "atac_feature"
chic_feat <- unique(chic[, .(feature, Chr, start_minus5, end_plus5)])
names(chic_feat)[1] = "chic_feature"

setkey(atac_feat, Chr, start, end)

at_chic <- foverlaps(chic_feat, atac_feat, by.x = c("Chr", "start_minus5", "end_plus5"))
# there are a few with no atac peaks nearby, as we might expect

chic_na <- at_chic[is.na(atac_feature)]
length(unique(chic_na$chic_feature)) # was 5482 (now 5532 with WASP)

chic_nona <- at_chic[!is.na(atac_feature)]
length(unique(chic_nona$chic_feature)) # 5670 (now 5675 with WASP...)

## Now add the RNAseq. Based on ENSG ID (check)
at_chic[, c("DpnID", "ENSG_ID", "Gene") := tstrsplit(chic_feature, split = "_")]

rna[, c("ENSG_ID", "Gene") := tstrsplit(feature, split = "_")]
rna_feat <- unique(rna[, .(feature, ENSG_ID)])
names(rna_feat)[1] = "rna_feature"

at_chic_rna <- rna_feat[at_chic, on = "ENSG_ID"]
# There were only a couple that didnt match based on ENSG_ID.
# remove when they have neither atac or rna.
to_save1 <- at_chic_rna[!is.na(atac_feature) | !is.na(rna_feature)]

to_save <- unique(to_save1[, .(chic_feature, atac_feature, rna_feature)])

dir.create("./triplet_features_WASP")
fwrite(to_save, file = "./triplet_features_WASP/triplet_list_011223.txt", quote = FALSE, sep = "\t", 
       row.names = FALSE, col.names = TRUE)


