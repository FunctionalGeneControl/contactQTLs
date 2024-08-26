###### Script to run on all GUESS and BaseQTL variants combined.
###### Will use the combined universe between BaseQTL and GUESS.
###### Using revision data 5th May 2024

library(data.table)
library(ReMapEnrich)
library(stringr)
#library(dplyr)
#library(tidyr)

###### Filter for DpnII and distal

setwd("~/HRJ_monocytes/remapenrich/Rev")

##########################################################################################
##########################################################################################
# READ IN FILES
#testGUESS <- fread("~/HRJ_monocytes/leo_triplets/results/mar23/filtered/final_biorxiv/final_GUESS_results_DpnII_distal_betas_AI_hg38_forRemap_CORRECTED.txt")
testGUESS <- fread("~/HRJ_monocytes/leo_triplets/results/Rev/filtered/final_GUESS_results_DpnII_distal_betas_AI_hg38_forRemap_CORRECTED.txt")
testGUESS[, position0 := position-1]
setnames(testGUESS, "Post-proc. FDR mPPI SNP single", "SNP")
backgroundGUESS <- fread("~/HRJ_monocytes/leo_triplets/results/mar23/filtered/GUESS_background_forRemap.txt")
setnames(backgroundGUESS, "id", "SNP")
backgroundGUESS[, set := NULL]
                
testBASE <- fread("~/HRJ_monocytes/BaseQTL/findings_round2_dpnIICorrection/Sig99_for_remap.txt")
names(testBASE) = c("SNP", "chrom", "position", "position0")
backgroundBASE <- fread("~/HRJ_monocytes/BaseQTL/findings_round2_dpnIICorrection/background_for_remap.txt")
names(backgroundBASE) = c("SNP", "chrom", "position", "position0")

test <- unique(rbind(testGUESS, testBASE))
background <- unique(rbind(backgroundGUESS, backgroundBASE))

##########################################################################################
##########################################################################################


############################## Run remapenrich ###########################################
query <- makeGRangesFromDataFrame(test, seqnames.field= "chrom", start.field= "position0", 
                                  end.field = "position", starts.in.df.are.0based = TRUE, ignore.strand = TRUE, keep.extra.columns = TRUE)
background_query <- makeGRangesFromDataFrame(background, seqnames.field= "chrom", start.field= "position0", 
                                             end.field = "position", starts.in.df.are.0based = TRUE, ignore.strand = TRUE, keep.extra.columns = TRUE)

# Using all cell types. Do not separate the cell:TF
remap <- fread("~/HRJ_monocytes/external_data/remap/remap2022_nr_macs2_hg38_v1_0.bed")
names(remap)[1:5] = c("seqnames", "start", "end", "id", "score") 
# remove cell type?
#remap[, TF := gsub("(.+?)(\\:.*)", "\\1", id)]
#remap2 <- remap[, .(V1, V2, V3, TF, V5)]
#names(remap2) = c("seqnames", "start", "end", "id", "score") # have to name the ID col for it to work.

## or mono. Here it would makde more sense to separate the cell:TF
#mono <- fread("~/HRJ_monocytes/external_data/remap/remap2022_nr_macs2_hg38_v1_0_monocyte_only.bed")
##names(mono)[1:5] = c("seqnames", "start", "end", "id", "score") # have to name the ID col for it to work.
## or remove the cell type?
#mono[, TF := gsub("(.+?)(\\:.*)", "\\1", V4)]
#mono2 <- mono[, .(V1, V2, V3, TF, V5)]
#names(mono2) = c("seqnames", "start", "end", "id", "score") # have to name the ID col for it to work.

#################
set.seed(1840) ##
#################


######## All cell types
catalog <- makeGRangesFromDataFrame(remap, keep.extra.columns = TRUE)

enrichment.df_shuffles100 <- enrichment(query, catalog, byChrom = TRUE, fractionCatalog=0, nCores = 1, 
                                        universe = background_query, shuffles = 100) 
enrichment.dt <- as.data.table(enrichment.df_shuffles100)
enrichment.dt[1:30, c("category", "q.value", "p.value", "effect.size")]
sig <- enrichment.dt[q.value < 1e-5]
fwrite(enrichment.dt, file = "./GUESS_BASE_vUniverseSNPs_remap_shuffles100.txt", sep = "\t", quote = F, row.names = F, col.names = T)

pdf(file = "./GUESS_BASE_vUniverseSNPs_remap_shuffles100_barplot.pdf")
par(mfrow = c(1, 1), mar = c(4, 10, 6, 1), oma = c(0.5, 10, 1, 0.5), mgp = c(2.2, 0.7, 0))
enrichmentBarPlot(enrichment.df_shuffles100, sigDisplayQuantile = 0.5, top = 27, aRisk = 0.00001)
dev.off()
pdf(file = "./GUESS_BASE_vUniverseSNPs_remap_shuffles100_dotplot.pdf")
par(mfrow = c(1, 1), mar = c(4, 10, 6, 1), oma = c(0.5, 10, 1, 0.5), mgp = c(2.2, 0.7, 0))
enrichmentDotPlot(enrichment.df_shuffles100, top = 7) # these were not all sig, only up to RELA
dev.off()
enrichmentVolcanoPlot(na.omit(enrichment.df_shuffles100), sigDisplayQuantile = 0.9, aRisk = 0.00001, pch=16)


# what are the cell lines?
#MV4-11 is macrophage like cell line
#TSU-1621MT is myeloid leukaemia - include it (myeloid cells are monocytes and granulocytes)
#RWPE-1 is prostate epithelium
#501-mel is malignant melanoma
#AML is myeloid - include it
#HL-60 is promyeoblast (myeloid) - include it
#CTV-1 is T cell but disputed
#HUVEC - no
#MCF−7,HAEC,MCF10A−Er−Src,endothelial,Karpas−299,Hep−G2

### Check the results.
### We can use q.significance > 5 as the baseline.
### This is at an FDR of 5%.
enrichment.dt <- fread("./GUESS_BASE_vUniverseSNPs_remap_shuffles100.txt")
#enrichment.dt[e.value < 0.05]
enrichment.dt[q.significance > 5]
enrichment.df <- as.data.frame(enrichment.dt)

enrichment.dt[category %like% "monocyte", Cell_type := "Monocyte"]
enrichment.dt[category %like% "THP-1", Cell_type := "Monocyte"]
enrichment.dt[category %like% "CD14", Cell_type := "Monocyte"]
enrichment.dt[category %like% "HL-60", Cell_type := "Myeloid"]
enrichment.dt[category %like% "SKH1", Cell_type := "Myeloid"]
enrichment.dt[category %like% "macrophage", Cell_type := "Macrophage"]
enrichment.dt[is.na(Cell_type), Cell_type := "Other"]
setnames(enrichment.dt, "category", "Category")
enrichment.dt[, category := tstrsplit(Category, split = ":", keep = 1)]
enrichment.dt

# proportion of cell types in remap versus in top TFs?
enrichment.dt[, sig := q.significance > 5]
sig_dt <- enrichment.dt[sig == T]
# todo.

# 1. AR:THP-1 -> monocyte cell line
# 2. SKI:HL-60 -> myeloid cell line
# 3. RXR:macrophage -> macrophage
# 4. GATA2:SKH1 -> myeloid cell line
# 5. STAT1:CD14 -> monocytes
# 6. DAXX:PC-3 -> prostate cancer cell line
# 7. SOX2:NPC,HNSC -> Nasopharyngeal carcinoma, human neural stem cells
# 8. FOSL1:MDA-MB-231,BT-549,HCT-116,MNNG-HOS -> other cancers
# 9. FOSL2:A-549,MDA-MB-231,NPC,LPS141,SK-N-SH -> other cancers
# 10. EP300:fibroblast,A-549,SK-N-SH -> other
# 11. YY1:SK-N-SH,RH4,HCT-116 -> other cancers
# 12. JUN:HAEC,endothelial,MCF10A-Er-Src,BT-549,HUES-8,786-O -> other cancers
# 13. TEAD4:A-549,MDA-MB-231,hMSC-TERT4,PC-9,SK-MEL-147,SK-N-SH,HCT-116,Ishikawa -> other cancers
# 14. RBPJ:GIC,GSC8-11,SCC,HKC,HCC1599,MDA-MB-157 -> other cancers
# 15. MED1:cardiomyocyte,SGBS,SUM159PT,U-87MG,hMSC-TERT4 -> other cancers
# 16. RELA:Detroit-562,SW480,786-O,HAEC,HUVEC-C,SGBS,AC16,FaDu,786-M1A,aortic-endothelial-cell,mammary-epithelial-cell,IMR-90 -> other cancer
# 17. KMT2B:AML -> myeloid
# 18. SREBP2:monocyte -> monocyte
# 19. MYB:THP-1 -> monocyte
# 20. NR3C1:BEAS-2B,Ishikawa,U2OS,IMR-90,MDA-MB-231,SUM159PT -> other
# 21. BRD4:LPS141,IMR-90,Hs-352-Sk,HCC1395,MPNST,RH4,SUM159PT,MDA-MB-231,HUVEC-C
# 22. STAT3:HCC1187,A139,MDA-MB-231,MDA-MB-468,MDA-MB-157,MCF-7,MCF10A-Er-Src,SUM159PT,NCI-H358,MCF-10A,T-47D,HCC1143,HCC70,A-137 -> other
# 23. JMJD1C:HL-60 -> myeloid
# 24. EGR1:macrophage -> macrophage
# 25. MED1:RH4,hMSC-TERT4 -> other
# 26. BRD4:LNCaP-C4-2,HCC1395 -> other
# 27. CEBPB:monocyte,MV4-11,THP-1,HL-60 -> monocyte

sig_dt[, overall_category := Category]
sig_dt[, c("Category", "Cell_type", "category") := NULL]

enrichment.df.sig <- as.data.frame(sig_dt)
row.names(enrichment.df.sig) <- sig_dt$overall_category

### I was getting some error about colors. Now just using the dt -> df that was first saved.
enrichment.dt <- fread("./GUESS_BASE_vUniverseSNPs_remap_shuffles100.txt")
enrichment.df <- as.data.frame(enrichment.dt)

pdf(file = "./GUESS_BASE_vUniverseSNPs_remap_shuffles100_dotplot.pdf", height = 7, width = 5)
par(mfrow = c(1, 1), mar = c(4, 10, 6, 1), oma = c(0.5, 10, 1, 0.5), mgp = c(2.2, 0.7, 0))
enrichmentDotPlot(enrichment.df, inches = 1/7) # these were all the sig ones.
dev.off()

#sig[, category := paste(category, Cell_type, sep = "_")]
pdf(file = "./GUESS_BASE_vUniverseSNPs_remap_shuffles100_barplot.pdf")
par(mfrow = c(1, 1), mar = c(4, 10, 6, 1), oma = c(0.5, 10, 1, 0.5), mgp = c(2.2, 0.7, 0))
enrichmentBarPlot(enrichment.df.sig, sigDisplayQuantile = 1, top = 27, aRisk = 0.00001)
dev.off()

### Then run remapenrich agains the following combined datasets:
# TOBIAS
# remap mono
# MaxATAC
# CTCF ChIP mono
### need to think about this - can we have overlaps?
### I think we keep the dataset that these TFs come from. In remap, the ChIP-seq datasets 
### are made into non-redundant peaks (merge similar targets)
### Keep everything separate and see how this looks
mono_remap <- fread("~/HRJ_monocytes/external_data/remap/remap2022_nr_macs2_hg38_v1_0_monocyte_only.bed")
mono_remap_bed1 <- mono_remap[, .(V4, V1, V2, V3)]
names(mono_remap_bed1) = c("category", "chrom", "start", "end")

tobias <- fread("~/HRJ_monocytes/ATAC-seq/TOBIAS_snakemake/clean_mono_34reps/overview/all_monocytes_bound_format.bed")
tobias[, category := paste(TF, "TOBIASmono", sep = ":")]
tobias_bed1 <- tobias[, .(category, chrom, start, end)]

# MaxATAC: combine TFs
#TFs <- fread("/rds/general/user/hrayjone/home/opt/maxatac/data/models/TF.list", header = F)
#TFs <- TFs[V1 != "TF.list"]
#TFs <- TFs[, V1]
#DIR <- file.path("~/HRJ_monocytes/ATAC-seq/maxatac/output")

#allTFs <- data.table()
#for(myTF in TFs) {
#  myTF_filename <- paste0(DIR, "/", myTF, "/", myTF, ".averaged.allChr_maxatac_predict_32bp.bed")
#  tf_dt <- fread(myTF_filename)
#  names(tf_dt) = c("chrom", "start", "end", "maxatac_score")
#  tf_dt[, TF := myTF]
#  allTFs <- rbind(allTFs, tf_dt)
#}
#fwrite(allTFs, file = "~/HRJ_monocytes/ATAC-seq/maxatac/combined_scores_cQTLs/allTFs_maxatac.txt", 
#       sep = "\t", quote = F, row.names = F, col.names = T)
allTFs <- fread("~/HRJ_monocytes/ATAC-seq/maxatac/combined_scores_cQTLs/allTFs_maxatac.txt")
allTFs[, category := paste(TF, "MaxATACmono", sep = ":")]
maxatac_bed1 <- allTFs[, .(category, chrom, start, end)]

# CTCF ChIP-seq
ctcf <- fread("~/HRJ_monocytes/ChIP/ChIP_AP/mono_CTCF_all_peaks_calculated_for_UCSC.narrowPeak", skip = 1)
setnames(ctcf, c("V1", "V2", "V3"), c("chrom", "start", "end"))
ctcf[, category := "CTCF:ChIPmono"]
ctcf_bed1 <- ctcf[, .(category, chrom, start, end)]

allMono <- rbind(mono_remap_bed1, tobias_bed1, maxatac_bed1, ctcf_bed1)
fwrite(allMono, file = "./mono_remap_input_fourSets.txt", sep = "\t", 
       quote = F, row.names = F, col.names = T)

# what is the format actually?
# we need chr, start, end, factor, score, strand. Will see if we can do without score and strand.
allMono <- fread("./mono_remap_input_fourSets.txt")
allMono_bed <- allMono[, .(chrom, start, end, category)]
names(allMono_bed)[1:4] = c("seqnames", "start", "end", "id") 

########## Run remapenrich
testGUESS <- fread("~/HRJ_monocytes/leo_triplets/results/Rev/filtered/final_GUESS_results_DpnII_distal_betas_AI_hg38_forRemap_CORRECTED.txt")
testGUESS[, position0 := position-1]
setnames(testGUESS, "Post-proc. FDR mPPI SNP single", "SNP")
backgroundGUESS <- fread("~/HRJ_monocytes/leo_triplets/results/mar23/filtered/GUESS_background_forRemap.txt")
setnames(backgroundGUESS, "id", "SNP")
backgroundGUESS[, set := NULL]

testBASE <- fread("~/HRJ_monocytes/BaseQTL/findings_round2_dpnIICorrection/Sig99_for_remap.txt")
names(testBASE) = c("SNP", "chrom", "position", "position0")
backgroundBASE <- fread("~/HRJ_monocytes/BaseQTL/findings_round2_dpnIICorrection/background_for_remap.txt")
names(backgroundBASE) = c("SNP", "chrom", "position", "position0")

test <- unique(rbind(testGUESS, testBASE))
background <- unique(rbind(backgroundGUESS, backgroundBASE))

##########################################################################################
##########################################################################################


############################## Run remapenrich ###########################################
query <- makeGRangesFromDataFrame(test, seqnames.field= "chrom", start.field= "position0", 
                                  end.field = "position", starts.in.df.are.0based = TRUE, ignore.strand = TRUE, keep.extra.columns = TRUE)
background_query <- makeGRangesFromDataFrame(background, seqnames.field= "chrom", start.field= "position0", 
                                             end.field = "position", starts.in.df.are.0based = TRUE, ignore.strand = TRUE, keep.extra.columns = TRUE)

#################
set.seed(1840) ##
#################


######## All cell types
catalog <- makeGRangesFromDataFrame(allMono_bed, keep.extra.columns = TRUE)

enrichment.df_shuffles100 <- enrichment(query, catalog, byChrom = TRUE, fractionCatalog=0, nCores = 1, 
                                        universe = background_query, shuffles = 100) 
enrichment.dt <- as.data.table(enrichment.df_shuffles100)
enrichment.dt[1:30, c("category", "p.value", "effect.size", "q.significance")]
sig <- enrichment.dt[q.value < 1e-5] # same as q.significance > 5. We have 46 enrichments!
enrichment.dt[category %like% "CTCF:ChIPmono"] # not enriched
fwrite(enrichment.dt, file = "./GUESS_BASE_vUniverseSNPs_monoChIP_shuffles100.txt", sep = "\t", quote = F, row.names = F, col.names = T)

pdf(file = "./GUESS_BASE_vUniverseSNPs_monoChIP_shuffles100_barplot.pdf", height = 9, width = 6)
par(mfrow = c(1, 1), mar = c(4, 10, 6, 1), oma = c(0.5, 10, 1, 0.5), mgp = c(2.2, 0.7, 0))
enrichmentBarPlot(enrichment.df_shuffles100, sigDisplayQuantile = 0.5, top = 46, aRisk = 0.00001)
dev.off()
pdf(file = "./GUESS_BASE_vUniverseSNPs_monoChIP_shuffles100_dotplot.pdf", height = 9, width = 6)
par(mfrow = c(1, 1), mar = c(4, 10, 6, 1), oma = c(0.5, 10, 1, 0.5), mgp = c(2.2, 0.7, 0))
enrichmentDotPlot(enrichment.df_shuffles100, top = 46, minCircleSize = 0.01) # these were not all sig, only up to RELA
dev.off()





