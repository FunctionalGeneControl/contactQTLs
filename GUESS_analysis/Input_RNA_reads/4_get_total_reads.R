### Need to use the same method as Novogene - they mapped the reads to the GRCh38 genome using STAR v2.6.1d allowing for 2 mismatches and assigned to Ensembl Release 94 genes. For Matrix eQTL, total gene counts were determined using FeatureCounts v1.5.0-p3 with default settings in the Novogene analysis pipeline.

### Run in DT_DPLYR environment, which has rsubread installed.
### We use a GTF file from Ensembl.

### "Note that, when counting at the meta-feature level, reads that overlap multiple features of the same meta-feature are always counted exactly once for that meta-feature, provided there is no overlap with any other meta-feature. For example, an exon-spanning read will be counted only once for the corresponding gene even if it overlaps with more than one exon."

library(Rsubread)
library(data.table)

setwd("~/spivakov/HRJ_monocytes/eqtls/WASP/mono_34reps/unique")

#### Testing
#setwd("~/spivakov/HRJ_monocytes/eqtls/WASP/mono_34reps/downsample_unique")
#bamFiles <- list.files(pattern = "\\.DS.bam$")
####

bamFiles <- list.files(pattern = "\\.sort.bam$")
#bamFiles

##### Run feature counts for the sample
##### Testing using a downsampled bam file.
#mySample <- file.path("~/spivakov/HRJ_monocytes/eqtls/WASP/mono_34reps/downsample_unique/S025NM-01.keep.merge.unique.sort.DS.bam")

#

###### Now trying the GTF formatted file that I have downloaded from ensembl FTP server. It is in hg38 and contains exons. Try with 1) no multi assignment and 2) multi assignment.


myGTF <- file.path("~/HRJ_monocytes/external_data/ensembl/Homo_sapiens.GRCh38.94.gtf.gz")

mycounts <- featureCounts(bamFiles,
                         annot.inbuilt = "hg38",
                         annot.ext = myGTF, 
			 isGTFAnnotationFile = TRUE, 
                         allowMultiOverlap = FALSE,
                         isPairedEnd = TRUE,
                         nthreads = 8)
## Note, the object also has: $targets - what were the bam files; $annotation - where were the transcripts, and $stat - how many reads were aligned per sample

saveRDS(mycounts, file = "../total_counts/featureCounts_mono34reps_GTF.Rds")
counts_dt <- as.data.table(mycounts$counts, keep.rownames = "GeneID")

### Modify the names of the counts table.
counts_novo <- fread("~/novogene/X204SC20081070-Z01-F001/X204SC20081070-Z01-F001_Homo_sapiens_result/3.Quantification/Count/readcount.txt", 
                     sep = "\t", header = TRUE)
str(counts_novo)

### Change the names to be the same as in novogene table: gene_id and then the sample names but without "-N_etc"
names2Change <- names(counts_dt[, 2:35])
namesChanged <- gsub("-.*", "", names2Change)
allNames <- c("gene_id", namesChanged)
#allNames
names(counts_dt) = allNames

print(head(counts_dt))

fwrite(counts_dt, file = "../total_counts/featureCounts_mono34reps_GTF.txt", sep = "\t", quote = F, row.names = F, col.names = T)
rm(counts_dt)
rm(mycounts)


######################################################################
#### Not using: counts allowing for multiple assignment
#mycountsmulti <- featureCounts(bamFiles,
#                         annot.inbuilt = "hg38",
#                         annot.ext = myGTF, 
#			 isGTFAnnotationFile = TRUE, 
#                         allowMultiOverlap = TRUE,
#                         isPairedEnd = TRUE,
#                         nthreads = 8)
#
#saveRDS(mycountsmulti, file = "../total_counts/featureCounts_mono34reps_multiFeatures_GTF.Rds")
#counts_dt <- as.data.table(mycountsmulti$counts, keep.rownames = "GeneID")
#fwrite(counts_dt, file = "../total_counts//featureCounts_mono34reps_multiFeatures_GTF.txt", sep = "\t", quote = F, row.names = F, col.names = T)
#rm(counts_dt)
#rm(mycountsmulti)
#


