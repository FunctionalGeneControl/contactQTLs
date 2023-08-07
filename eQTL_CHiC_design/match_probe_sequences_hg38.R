# 01/07/2019

library(dplyr)
library(tidyr)
library(data.table)

## Read in filtered datasets

eqtls <- fread(file = "~/eCHiC/design/final_design/final_filtered_eqtls_illumina_july19_hg38.txt", 
                    sep = "\t", header = TRUE)


genes <- fread(file = "~/eCHiC/design/final_design/all_captured_genes_finalhg38.txt", 
                    sep = "\t", header = TRUE)


frags <- fread(file = "~/eCHiC/design/final_design/dpnIIhg38.bed", sep = "\t", header = FALSE)
names(frags) = c("Chrom", "DpnII_start", "DpnII_end", "DpnII_index")

## Get the probe sequences

setwd("~/eCHiC/design/probes")
fail <- fread(file = "failed_grch38_fragments.txt", sep = "\t", header = FALSE)
names(fail) = c("Chrom", "DpnII_start", "DpnII_end")
fail$Frag_length = fail$DpnII_end - fail$DpnII_start
fail2 <- fail[Frag_length < 2000]
median(fail2$Frag_length)
hist(fail2$Frag_length, breaks = 20)

success <- fread(file = "probe_grch38_positions.txt", sep = "\t", header = FALSE)
names(success) = c("Chrom", "Probe_start", "Probe_end", "Strand", "Sequence", "DpnII_start", "DpnII_end")
success$Frag_length = success$DpnII_end - success$DpnII_start
success$Chrom = paste("chr", success$Chrom, sep="")
success2 <- success[Frag_length < 2000]
median(success2$Frag_length)
hist(success2$Frag_length, breaks = 20)


all <- fread(file = "all_grch38_restriction_fragments.txt", sep = "\t", header = FALSE)
names(all) = c("Chrom", "DpnII_start", "DpnII_end")
all$Frag_length = all$DpnII_end - all$DpnII_start

all2 <- all[Frag_length < 2000]
median(all2$Frag_length)
hist(all2$Frag_length, breaks = 20)



# Match required DpnII fragments to probes

length(unique(frags$DpnII_index))

with_probes <- frags[success, on = c("Chrom", "DpnII_start", "DpnII_end"), nomatch = NULL]
with_probes
length(unique(with_probes$DpnII_index))
length(unique(with_probes$Sequence))


#############################################################################

# How many eQTLs can we capture?

eqtls2 <- unique(eqtls[, list(Chrom, hg38SNP_ID, DpnII_start, DpnII_end, DpnII_index)])
names(eqtls2) = c("Chrom", "Feature", "DpnII_start", "DpnII_end", "DpnII_index")

eqtls3 <- unique(eqtls[, list(Chrom, hg38SNP_ID, hg38Proxy_ID, DpnII_start, DpnII_end, DpnII_index)])

genes2 <- unique(genes[, list(Chrom, cap_ENSG_ID, gen_DpnII_start, gen_DpnII_end, gen_DpnII_index, cap_TSS_start)])
names(genes2) = c("Chrom", "Feature", "DpnII_start", "DpnII_end", "DpnII_index", "TSS")

genes3 <- unique(genes[, list(Chrom, cap_ENSG_ID, cap_ENST_ID, gen_DpnII_start, gen_DpnII_end, gen_DpnII_index)])

eqtls_genes <- unique(rbind(eqtls2, genes2))
length(unique(eqtls_genes$DpnII_index))

with_probes_eqtls <- eqtls2[success, on = c("Chrom", "DpnII_start", "DpnII_end"), nomatch = NULL]
with_probes_eqtls
length(unique(with_probes_eqtls$DpnII_index))
length(unique(with_probes_eqtls$Sequence))
length(unique(with_probes_eqtls$Feature))

# 1,458 eQTLs can be captured; 2,343 DpnII fragments

Proxies_captured <- eqtls3[success, on = c("Chrom", "DpnII_start", "DpnII_end"), nomatch = NULL]
length(unique(Proxies_captured$hg38Proxy_ID))
length(unique(Proxies_captured$hg38SNP_ID))

# 2,571 proxies

with_probes_genes <- genes2[success, on = c("Chrom", "DpnII_start", "DpnII_end"), nomatch = NULL]
with_probes_genes
length(unique(with_probes_genes$DpnII_index))
length(unique(with_probes_genes$Sequence))
length(unique(with_probes_genes$Feature))
length(unique(with_probes_genes$TSS))

# 8,533 genes can be captured; 15,921 DpnII fragments; 32,483 TSS

## Get distance of captured eQTLs to target genes, for Report

temp_with_probes_eqtls <- with_probes_eqtls
names(temp_with_probes_eqtls)[2] = "hg38SNP_ID"
cap_snps <- eqtls[temp_with_probes_eqtls, on = "hg38SNP_ID", nomatch = NULL, allow.cartesian = TRUE]
dist_cap_snps <- unique(cap_snps[, list(hg38SNP_ID, Gene, TSS, SNP_TSSdist)])
range(dist_cap_snps$SNP_TSSdist)
hist(dist_cap_snps$SNP_TSSdist, breaks = 50)

# Write a file for future information..

write.table(cap_snps, file = "~/eCHiC/design/final_design/detailed_final_captured_eQTLs.txt", 
              col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

# EDIT 19/04/2020 Writing a file showing exactly which TSS were captured directly

write.table(with_probes_genes, file = "~/eCHiC/design/final_design/detailed_final_captured_TSS.txt", 
            col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)


# Percentage of SNPs > 10kb 
ten <- dist_cap_snps[SNP_TSSdist > 10000]
length(unique(ten$hg38SNP_ID)) /length(unique(cap_snps$hg38SNP_ID)) *100

# Percentage SNPs > 20 kb
twenty <- dist_cap_snps[SNP_TSSdist > 20000]
length(unique(twenty$hg38SNP_ID)) /length(unique(cap_snps$hg38SNP_ID)) *100

# Percentage SNPs > 30kb
thirty <- dist_cap_snps[SNP_TSSdist > 30000]
length(unique(thirty$hg38SNP_ID)) /length(unique(cap_snps$hg38SNP_ID)) *100

# Percentage SNPs > 40 kb
forty <- dist_cap_snps[SNP_TSSdist > 40000]
length(unique(forty$hg38SNP_ID)) /length(unique(cap_snps$hg38SNP_ID)) *100

# Percentage SNPs > 50 kb
fifty <- dist_cap_snps[SNP_TSSdist > 50000]
length(unique(fifty$hg38SNP_ID)) /length(unique(cap_snps$hg38SNP_ID)) *100

# Percenrage SNPs < 250kb of gene promoters
close <- dist_cap_snps[SNP_TSSdist < 250000]
length(unique(close$hg38SNP_ID)) /length(unique(cap_snps$hg38SNP_ID)) *100

# How many control gene promoters are there?
# Do an antijoin between the final gene/probe file and eGenes

cap_genes <- with_probes_genes
names(cap_genes)[2] = "ENSG_ID"
eqtls

cap_non_target <- cap_genes[!eqtls, on = .(ENSG_ID)]
length(unique(cap_non_target$DpnII_index))
length(unique(cap_non_target$ENSG_ID))
length(unique(cap_non_target$TSS))

# 11,063 "non-target" gene promoter fragments (20,443 TSS corresponding to 6,760 genes) 

cap_target <- cap_genes[eqtls, on = .(ENSG_ID), allow.cartesian = TRUE, nomatch = NULL]
length(unique(cap_target$DpnII_index))
length(unique(cap_target$ENSG_ID))
length(unique(cap_target$i.TSS))
unique(cap_target)

# 5,435 "target" gene promoter fragments (6,213 TSS corresponding to 1,773 genes)


### Make design file for Agilent

# 1. Get Target IDs

eqtl_feat <- unique(eqtls[, list(Chrom, hg38SNP_ID, DpnII_start, DpnII_end, DpnII_index)])
names(eqtl_feat) = c("Chrom", "Target_ID_temp", "DpnII_start", "DpnII_end", "DpnII_index")
genes_feat <- unique(genes[, list(Chrom, cap_Gene, gen_DpnII_start, gen_DpnII_end, gen_DpnII_index)])
names(genes_feat) = c("Chrom", "Target_ID_temp", "DpnII_start", "DpnII_end", "DpnII_index")

test <- rbind(eqtl_feat, genes_feat)

ids <- aggregate(test$Target_ID_temp, list(test$DpnII_index), paste, collapse=";")
ids <- as.data.table(ids)
names(ids) = c("DpnII_index", "Target_ID")

with_Target_IDs <- test[ids, on = "DpnII_index"]
with_Target_IDs$test <- nchar(with_Target_IDs$Target_ID)
range(with_Target_IDs$test)

final <- unique(with_Target_IDs[, list(Chrom, DpnII_start, DpnII_end, DpnII_index, Target_ID)])
final

# 2. Get probe IDs 

final$temp <- paste(final$Chrom, final$DpnII_start, sep = ":")
final$Probe_ID <- paste(final$temp, final$DpnII_end, sep = "-")
final$test <- nchar(final$Probe_ID)
range(final$test)
frags_ids <- unique(final[, list(Chrom, DpnII_start, DpnII_end, Probe_ID, Target_ID)])


write.table(frags_ids, file = "~/eCHiC/design/final_design/ProbeIDs_TargetIDs_to_capturehg38.txt", 
            col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

# 3. Combine with probe sequences

agilent_design <- frags_ids[success, on = c("Chrom", "DpnII_start", "DpnII_end"), nomatch = NULL]


## Select one probe at random per restriction fragment

agilent_design_temp <- agilent_design[,.SD[sample(.N, min(1,.N))],by = Probe_ID]

agilent_design_temp

length(unique(agilent_design_temp$Probe_ID))
# 18,201 ... Sequence = 18,195. This is because some of the probes can potentially bind two regions in the genome.

dupls <- agilent_design %>% group_by(Sequence) %>% filter(n() > 1)

agilent_design$dups <- agilent_design$Probe_ID %in% dupls$Probe_ID 

dups2check <- agilent_design[dups == "TRUE"]


## Need to remove the duplicated seqeunces, if there are no alternatives.

write.table(dups2check, file = "~/eCHiC/design/probes/dup_probes2check.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

dup_probes2remove <- fread(file = "dup_probes2remove.txt", header = TRUE, sep = "\t")
probes2remove <- dup_probes2remove[Unique_in_UCSC == "FALSE"]
probes2remove

agilent_design_unique <- agilent_design[!agilent_design$Sequence %in% probes2remove$Sequence ,]

length(unique(agilent_design$Sequence))
length(unique(probes2remove$Sequence))
length(unique(agilent_design_unique$Sequence))

## Random sample the probes now

agilent_design_final <- agilent_design_unique[,.SD[sample(.N, min(1,.N))],by = Probe_ID]

length(agilent_design_final$Probe_ID)
length(agilent_design_final$Sequence)

#### 18,178 probes/fragments

## Make file for WashU

for_washu <- agilent_design_final[, list(Chrom, DpnII_start, DpnII_end, Target_ID)]
write.table(for_washu, file = "~/eCHiC/design/final_design/Frags_can_be_captured_hg38.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)


################## 4-column format

agilent_design_final_4col <- agilent_design_final[, list(Target_ID, Probe_ID, Sequence, "1")]
names(agilent_design_final_4col) = c("TargetID", "ProbeID", "Sequence", "Replication")

write.table(agilent_design_final_4col, file = "~/eCHiC/design/probes/eCHiC_4col.txt", col.names = TRUE, row.names = FALSE, quote= FALSE, sep = "\t")


## to try
#agilent_design_test_4col <- agilent_design_final_4col
#agilent_design_test_4col$ProbeID <- seq.int(nrow(agilent_design_test_4col))

#write.table(test, file = "~/eCHiC/design/probes/test_4col.txt", col.names = TRUE, row.names = FALSE, quote= FALSE, sep = "\t")


################## 6-column format


agilent_design_final_6col <- agilent_design_final[, list(Target_ID, Probe_ID, Sequence, "1", Strand, Chrom, Probe_start, Probe_end)]

agilent_design_final_6col$temp <- paste(agilent_design_final_6col$Chrom, agilent_design_final_6col$Probe_start , sep = ":")

agilent_design_final_6col$Coordinates = paste(agilent_design_final_6col$temp, agilent_design_final_6col$Probe_end, sep = "-")

names(agilent_design_final_6col) = c("TargetID", "ProbeID", "Sequence", "Replication", "Strand", "Chrom", 
                                     "Probe_start", "Probe_end", "temp", "Coordinates")

col6 <- agilent_design_final_6col[, list(TargetID, ProbeID, Sequence, Replication, Strand, Coordinates)]


write.table(col6, file = "~/eCHiC/design/probes/eCHiC_6col.txt", col.names = TRUE, row.names = FALSE, quote= FALSE, sep = "\t")


########################## 8-column format

agilent_design_final_8col <- agilent_design_final[, list(Target_ID, Probe_ID, Sequence, "1", Strand, Chrom, Probe_start, Probe_end, 
                                                         (Probe_start-1), (Probe_end-1))]

names(agilent_design_final_8col) = c("TargetID", "ProbeID", "Sequence", "Replication", "Strand", "Chromosome", "Probe_start", "Probe_end", 
                                     "Start", "Stop")

col8 <- agilent_design_final_8col[, list(TargetID, ProbeID, Sequence, Replication, Strand, Chromosome, Start, Stop)]

write.table(col8, file = "~/eCHiC/design/probes/eCHiC_8col.txt", col.names = TRUE, row.names = FALSE, quote= FALSE, sep = "\t")




################################Just a check ################################################

### do a comparison - use the full probeset

hg38_probes <- agilent_design2$Sequence

setwd("~/eCHiC/design/probes")

hg19_to_capture <- fread("~/eCHiC/design/final_design/hg19/ProbeIDs_TargetIDs_to_capture.txt")
hg19_success <- fread(file = "~/eCHiC/design/probes/hg19/probe_positions.txt", sep = "\t", header = FALSE)
names(hg19_success) = c("Chrom", "Probe_start", "Probe_end", "Strand", "Sequence", "DpnII_start", "DpnII_end")
hg19_success$Chrom = paste("chr", hg19_success$Chrom, sep="")

hg19_agilent_design <- hg19_to_capture[hg19_success, on = c("Chrom", "DpnII_start", "DpnII_end"), nomatch = NULL]


length(unique(hg19_agilent_design$Sequence))

hg19_probes <- hg19_agilent_design$Sequence

table(hg38_probes %in% hg19_probes)

# They share ~82% of sequences


## Try with lifting over the dpnII from hg19

hg19_agilent_design

library(liftOver)

chain <- import.chain("~/bin/liftover/hg19ToHg38.over.chain")
to_lift3 <- makeGRangesFromDataFrame(hg19_agilent_design, keep.extra.columns = TRUE, 
                                     ignore.strand = TRUE, seqnames.field = "Chrom", start.field = "DpnII_start", 
                                     end.field = "DpnII_end")


design_lifted <- liftOver(to_lift3, chain)
design_lifted <- unique(as.data.table(design_lifted))
design_lifted

hg19coords <- unique(design_lifted[, list(start, end)])
hg38coords <- unique(agilent_design[, list(DpnII_start, DpnII_end)])

table(hg38coords$DpnII_start %in% hg19coords$start)
table(hg38coords$DpnII_end %in% hg19coords$end)

# They share ~86% of DpnII fragment coordinates


######### Checking my hg38 probes against Lera's


myprobes <- success[, list(Chrom, Probe_start, Probe_end, Strand, Sequence, DpnII_start, DpnII_end)]

leraprobes <- fread("~/eCHiC/design/probes/from_Lera/probe_positions.txt", header = FALSE, sep = "\t")
leraprobes
names(leraprobes) = c("Chrom", "Probe_start", "Probe_end", "Strand", "Sequence", "DpnII_start", "DpnII_end")
leraprobes$Chrom = paste("chr", leraprobes$Chrom, sep = "")
leraprobes

test <- leraprobes[myprobes, on = "Sequence", mult = "first", nomatch = NULL]
test

leraprobes$inmine <- leraprobes$Sequence %in% myprobes$Sequence
temp <- leraprobes[inmine == "FALSE"]
null <- leraprobes[inmine == "FALSE"]
