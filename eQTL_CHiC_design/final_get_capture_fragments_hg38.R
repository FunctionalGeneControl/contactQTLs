####### This script was written to design a capture system for
####### eQTLs in monocytes
####### With an aim to sequence the eQTLs, or their proxies, directly within
####### the capture library ditags (DpnII fragments)


####### eQTL data source = eQTL metaanalysis - described in Roman Kreuzhuber's PhD thesis, available at https://www.repository.cam.ac.uk/items/d96c7c98-3068-4211-9725-dc0d65a0c22a

####### Helen Ray-Jones 01/07/2019 

library(data.table)
library(liftOver)
library(dplyr)
library(tidyr)
library(illuminaHumanv4.db)
library(IRanges)

setwd("~/eCHiC/design")
eqtls <- fread(file = "./source/cd14_eqtls_no_freq.txt", sep = "\t", header = TRUE)

eqtls$chrom = paste("chr", eqtls$chrom, sep = "")
names(eqtls) = c("Array_Address_Id", "Gene", "Chrom", "SNP_pos", "SNP_ID", "Beta", "Beta_wrt_minor", "P_value",
                 "FDR", "Indep_signal_i", "ENSG_ID", "N_samples", "MAF", "HWE_P_value", "Imputation_score")



proxies <- fread(file = "./LD/all_cd14_eqtls.maf1.ld", header = TRUE)
proxies <- proxies[, list(CHR_A, BP_A, SNP_A, MAF_A, BP_B, SNP_B, MAF_B, R2)]
names(proxies) = c("Chrom", "SNP_pos", "SNP_ID", "SNP_MAF", "Proxy_pos", "Proxy_ID", "Proxy_MAF", "R2")
proxies$Chrom <- sub("^", "chr", proxies$Chrom)


ld <- proxies[eqtls, on = c("Chrom", "SNP_ID", "SNP_pos"), allow.cartesian = TRUE, nomatch = NULL]

range(ld$SNP_MAF)
range(ld$Proxy_MAF)
range(ld$R2)

## Restrict MAF to 10% and R2 to 0.9

ld2 <- ld[(SNP_MAF >= 0.1) & (Proxy_MAF >= 0.1) & (R2 >= 0.9)]

range(ld2$SNP_MAF)
range(ld2$Proxy_MAF)
range(ld2$R2)


### Lift them over to hg38


### Trying the remap tool on NCBI
snp_to_remap <- unique(ld2[, list(Chrom, SNP_pos, SNP_ID)])
snp_to_remap$SNP_start <- snp_to_remap$SNP_pos - 1
snp_to_remap <- snp_to_remap[, list(Chrom, SNP_start, SNP_pos, SNP_ID)]
write.table(snp_to_remap, "~/eCHiC/design/snp_to_remap.txt", col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
proxy_to_remap <- unique(ld2[, list(Chrom, Proxy_pos, Proxy_ID)])
proxy_to_remap$Proxy_pos <- as.numeric(proxy_to_remap$Proxy_pos)
proxy_to_remap$Proxy_start <- proxy_to_remap$Proxy_pos - 1
proxy_to_remap <- proxy_to_remap[, list(Chrom, Proxy_start, Proxy_pos, Proxy_ID)]
write.table(proxy_to_remap, "~/eCHiC/design/proxy_to_remap.txt", col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)

# Upload these files and liftover with Remap


remap_snps <- fread("~/eCHiC/design/source/hg38_stuff/report_snp_to_remap.txt.xls.html", header = TRUE, sep = "\t", fill=TRUE)

no_remap_snps <- remap_snps[mapped_int=="NULL"]

remap_proxies <- fread("~/eCHiC/design/source/hg38_stuff/report_proxy_to_remap.txt.xls.html", header = TRUE, sep = "\t", fill = TRUE)

no_remap_proxies <- remap_proxies[mapped_int=="NULL"]

###### Using Remap I cannot liftover 19 proxies and one SNP but the others are ok.
###### 109,667 proxies and 5370 SNPs.


remap_snps_success <- remap_snps[mapped_int!="NULL"]
names(remap_snps_success)[1] = "snp_name"

remap_proxies_success <- remap_proxies[mapped_int!="NULL"]
names(remap_proxies_success)[1] = "proxy_name"

ld2$snp_name <- paste(ld2$Chrom, ld2$SNP_pos, sep = "_")
ld2$proxy_name <- paste(ld2$Chrom, ld2$Proxy_pos, sep = "_") 

with_hg38_snps <- remap_snps_success[ld2, on = "snp_name", nomatch = NULL]
with_hg38_snps <- with_hg38_snps[, list(Chrom, SNP_ID, SNP_pos, snp_name, mapped_start, SNP_MAF, Proxy_ID, Proxy_pos, proxy_name, Proxy_MAF,
                                        R2, Array_Address_Id, Gene, Beta, Beta_wrt_minor, P_value, FDR, Indep_signal_i, ENSG_ID, N_samples, MAF, HWE_P_value, 
                                        Imputation_score)]
names(with_hg38_snps)[2] = "hg19SNP_ID"
names(with_hg38_snps)[3] = "hg19SNP_pos"
names(with_hg38_snps)[5] = "hg38SNP_pos"
names(with_hg38_snps)[7] = "hg19Proxy_ID"
names(with_hg38_snps)[8] = "hg19Proxy_pos"

with_hg38_snps_and_proxies <- remap_proxies_success[with_hg38_snps, on = "proxy_name", nomatch = NULL]
with_hg38_snps_and_proxies <- with_hg38_snps_and_proxies[, list(Chrom, hg19SNP_ID, hg19SNP_pos, snp_name, hg38SNP_pos, SNP_MAF, hg19Proxy_ID, hg19Proxy_pos, proxy_name, mapped_start, 
                                                                Proxy_MAF, R2, Array_Address_Id, Gene, Beta, Beta_wrt_minor, P_value, FDR, Indep_signal_i, ENSG_ID, N_samples, 
                                                                MAF, HWE_P_value, 
                                                                Imputation_score)]


names(with_hg38_snps_and_proxies)[10] = "hg38Proxy_pos"



with_hg38_snps_and_proxies_format1 <- separate(with_hg38_snps_and_proxies, col = "hg19Proxy_ID", into = c("Proxy_Chrom", "remove_me", "Proxy_ref", "Proxy_alt"), 
                           sep = ":", remove = FALSE)

with_hg38_snps_and_proxies_format2 <- with_hg38_snps_and_proxies_format1
with_hg38_snps_and_proxies_format2$hg38Proxy_ID <- paste(with_hg38_snps_and_proxies_format2$Proxy_Chrom, 
                                                         with_hg38_snps_and_proxies_format2$hg38Proxy_pos, 
                                                         with_hg38_snps_and_proxies_format2$Proxy_ref, 
                                                         with_hg38_snps_and_proxies_format2$Proxy_alt, sep = ":")

with_hg38_snps_and_proxies_format3 <- separate(with_hg38_snps_and_proxies_format2, col = "hg19SNP_ID", into = c("SNP_Chrom", "remove_me2", "SNP_ref", "SNP_alt"), 
                                               sep = ":", remove = FALSE)
with_hg38_snps_and_proxies_format3$hg38SNP_ID <- paste(with_hg38_snps_and_proxies_format3$SNP_Chrom, 
                                                         with_hg38_snps_and_proxies_format3$hg38SNP_pos, 
                                                         with_hg38_snps_and_proxies_format3$SNP_ref, 
                                                         with_hg38_snps_and_proxies_format3$SNP_alt, sep = ":")

remapped <- with_hg38_snps_and_proxies_format3[, list(Chrom, hg19SNP_ID, hg38SNP_ID, hg38SNP_pos, SNP_MAF, hg19Proxy_ID, hg38Proxy_ID, hg38Proxy_pos, 
                                                      Proxy_MAF, R2, Array_Address_Id, Gene, Beta, Beta_wrt_minor, P_value, FDR, Indep_signal_i, 
                                                      ENSG_ID, N_samples, MAF, HWE_P_value, Imputation_score)]


# 5,370 SNPs and 109,659 Proxies remain (some proxies were removed becuase of that one SNP).


######### Get TSS information



# convert "Array Address IDs" to "Illumina probe IDs"

adrToIllumina = toTable(illuminaHumanv4ARRAYADDRESS)
adrToIllumina = adrToIllumina[, c("ArrayAddress", "IlluminaID")]
colnames(adrToIllumina) = c("Array_Address_Id", "Probe_Id")
illuminaToSymbol = toTable(illuminaHumanv4SYMBOLREANNOTATED)
adrToSymbol = merge(adrToIllumina, illuminaToSymbol, by.x="Probe_Id", by.y="IlluminaID")
adrToSymbol = adrToSymbol[,c("Array_Address_Id", "SymbolReannotated")]
colnames(adrToSymbol) = c("Array_Address_Id", "Symbol")
negIl = mappedLkeys(revmap(illuminaHumanv4REPORTERGROUPNAME)["negative"])
negAdr = mappedRkeys(illuminaHumanv4ARRAYADDRESS[negIl])



adrToIllumina2 <- as.data.table(adrToIllumina)
adrToSymbol2 <- as.data.table(adrToSymbol)


remapped$Array_Address_Id <- as.character(remapped$Array_Address_Id)
adrToSymbol2$Array_Address_Id <- as.character(adrToSymbol2$Array_Address_Id)
adrToIllumina2$Array_Address_Id <- as.character(adrToIllumina2$Array_Address_Id)

with_symbols <- adrToSymbol2[remapped, on = "Array_Address_Id"]

with_ill_IDs <- adrToIllumina2[remapped, on = "Array_Address_Id"] 
length(unique(with_ill_IDs$hg38SNP_ID))
# 5370


########## Get TSS by joining with Biomart file (ensembl genes 96, which I downloaded with TSS and 
########## corresponding v4Illumina IDs)


TSS <- fread(file = "~/eCHiC/design/source/hg38_stuff/Ensembl96_GRCh38_TSS_and_Illumina_probes.txt", 
                  header = TRUE, sep = "\t", na.strings="")

names(TSS) = c("ENSG_ID", "ENST_ID", "Chrom", "TSS", "Gene", "v4_Probe_ID")

TSS2 <- TSS[, list(ENSG_ID, ENST_ID, Gene, TSS, Probe_Id=v4_Probe_ID)]
TSS3 <- unique(TSS2)

with_probes_and_TSS <- TSS3[with_ill_IDs, on = c("Probe_Id", "ENSG_ID"), allow.cartesian = TRUE]

with_probes_and_TSS[is.na(Gene)]
temp <- with_probes_and_TSS[!is.na(Gene)]
length(unique(temp$hg38SNP_ID))

length(unique(with_probes_and_TSS$hg38SNP_ID)) # Check that the SNPs are still all there
# 5370

# Keep track of where gene info has come from

with_probes_and_TSS <- with_probes_and_TSS[, list(Chrom, hg19SNP_ID, hg38SNP_ID, hg38SNP_pos, SNP_MAF, 
                                                  hg19Proxy_ID, hg38Proxy_ID, hg38Proxy_pos, Proxy_MAF, 
                                                  R2, ENSG_ID, ENST_ID, TSS, i.Gene, Gene,
                                                  Array_Address_Id, Probe_Id,  
                                                  Beta, Beta_wrt_minor, P_value, FDR, Indep_signal_i, 
                                                  N_samples, MAF, HWE_P_value, Imputation_score)]

names(with_probes_and_TSS) = c("Chrom", "hg19SNP_ID", "hg38SNP_ID", "hg38SNP_pos", "SNP_MAF", 
                               "hg19Proxy_ID", "hg38Proxy_ID", "hg38Proxy_pos", "Proxy_MAF", 
                               "R2", "ENSG_ID", "ENST_ID", "TSS", "rep_Gene", "Gene",
                               "Array_Address_Id", "Probe_Id",  
                               "Beta", "Beta_wrt_minor", "P_value", "FDR", "Indep_signal_i", 
                               "N_samples", "MAF", "HWE_P_value", "Imputation_score")

with_probes_and_TSS[is.na(TSS)]
with_probes_and_TSS[is.na(Gene)]

## Some Illumina_V4 probe Ids do not have associated TSS's in this file. In these cases, will need to
## assign TSS based on the genes reported in the original eQTL dataset. This occurs for 172 SNPs.

no_TSS_match <- with_probes_and_TSS[is.na(TSS)]
no_TSS_match$Probe_maps2TSS <- "FALSE"
length(unique(no_TSS_match$hg38SNP_ID))

## Get the TSS for these SNPs/ENSG

all_TSS2 <- unique(TSS[, list(ENSG_ID, ENST_ID, TSS, Gene)])


no_TSS_match_plus_TSS <- no_TSS_match[all_TSS2, on="ENSG_ID", nomatch = NULL]
no_TSS_match_plus_TSS <- no_TSS_match_plus_TSS[, list(Chrom, hg19SNP_ID, hg38SNP_ID, hg38SNP_pos, SNP_MAF, 
                                                      hg19Proxy_ID, hg38Proxy_ID, hg38Proxy_pos, Proxy_MAF, 
                                                      R2, ENSG_ID, i.ENST_ID, i.TSS, rep_Gene, 
                                                      i.Gene, Array_Address_Id, Probe_Id, Probe_maps2TSS, 
                                                      Beta, Beta_wrt_minor, P_value, FDR, Indep_signal_i, 
                                                      N_samples, MAF, HWE_P_value, Imputation_score)]

names(no_TSS_match_plus_TSS) = c("Chrom", "hg19SNP_ID", "hg38SNP_ID", "hg38SNP_pos", "SNP_MAF", 
                                 "hg19Proxy_ID", "hg38Proxy_ID", "hg38Proxy_pos", "Proxy_MAF", 
                                 "R2", "ENSG_ID", "ENST_ID", "TSS", "rep_Gene", 
                                 "Gene", "Array_Address_Id", "Probe_Id", "Probe_maps2TSS", 
                                 "Beta", "Beta_wrt_minor", "P_value", "FDR", "Indep_signal_i", 
                                 "N_samples", "MAF", "HWE_P_value", "Imputation_score")


no_TSS_match_plus_TSS[is.na(TSS)] # Should be empty
no_TSS_match_plus_TSS[is.na(Gene)]# Should be empty

length(unique(no_TSS_match_plus_TSS$hg38SNP_ID))

# Check those that didn't make it

no_TSS_match$check <- !(no_TSS_match$hg38SNP_ID %in% no_TSS_match_plus_TSS$hg38SNP_ID)

## match these based on Gene.

to_match_on_gene <- no_TSS_match[check == "TRUE"]

to_match_on_gene <- to_match_on_gene[, list(Chrom, hg19SNP_ID, hg38SNP_ID, hg38SNP_pos, SNP_MAF, 
                                            hg19Proxy_ID, hg38Proxy_ID, hg38Proxy_pos, Proxy_MAF, 
                                            R2, ENSG_ID, rep_Gene, Array_Address_Id, Probe_Id, Probe_maps2TSS, 
                                            Beta, Beta_wrt_minor, P_value, FDR, Indep_signal_i, 
                                            N_samples, MAF, HWE_P_value, Imputation_score)]

names(to_match_on_gene)[12] = "Gene"

no_TSS_match_plus_TSS2 <- to_match_on_gene[all_TSS2, on="Gene", nomatch = NULL]
length(unique(no_TSS_match_plus_TSS2$hg38SNP_ID))
# rescue 24 SNPs. Combine

no_TSS_match_plus_TSS2$rep_Gene = no_TSS_match_plus_TSS2$Gene

no_TSS_match_plus_TSS2 <- no_TSS_match_plus_TSS2[, list(Chrom, hg19SNP_ID, hg38SNP_ID, hg38SNP_pos, SNP_MAF, 
                                                        hg19Proxy_ID, hg38Proxy_ID, hg38Proxy_pos, Proxy_MAF, 
                                                        R2, i.ENSG_ID, ENST_ID, TSS, rep_Gene, Gene, Array_Address_Id, 
                                                        Probe_Id, Probe_maps2TSS, 
                                                        Beta, Beta_wrt_minor, P_value, FDR, Indep_signal_i, 
                                                        N_samples, MAF, HWE_P_value, Imputation_score)]
names(no_TSS_match_plus_TSS2)[11] = "ENSG_ID"
rescued <- rbind(no_TSS_match_plus_TSS, no_TSS_match_plus_TSS2)



#no_TSS_match$check2 <- !(no_TSS_match$hg38SNP_ID %in% rescued$hg38SNP_ID)
#whatrthese <- no_TSS_match[check2 == "TRUE"] # 5 snps where the genes do not exist in GRCh38 ensembl gene set. 
# For example, ENSG00000138041 is deprecated and no longer in the database.

rescued$target_TSS_source <- "Ensembl"



rescued2 <- rescued[, list(Chrom, hg19SNP_ID, hg38SNP_ID, hg38SNP_pos, SNP_MAF, 
                                                      hg19Proxy_ID, hg38Proxy_ID, hg38Proxy_pos, Proxy_MAF, 
                                                      R2, ENSG_ID, ENST_ID, TSS, rep_Gene, 
                                                      Gene, Array_Address_Id, Probe_Id, Probe_maps2TSS, target_TSS_source,  
                                                      Beta, Beta_wrt_minor, P_value, FDR, Indep_signal_i, 
                                                      N_samples, MAF, HWE_P_value, Imputation_score)]


### Now for those where the probe DOES match TSS

TSS_match <- with_probes_and_TSS[!is.na(TSS)]

TSS_match$Probe_maps2TSS <- "TRUE"
TSS_match$target_TSS_source <- "Illumina"

TSS_match <- TSS_match[, list(Chrom, hg19SNP_ID, hg38SNP_ID, hg38SNP_pos, SNP_MAF, 
                              hg19Proxy_ID, hg38Proxy_ID, hg38Proxy_pos, Proxy_MAF, 
                              R2, ENSG_ID, ENST_ID, TSS, rep_Gene, 
                              Gene, Array_Address_Id, Probe_Id, Probe_maps2TSS, target_TSS_source,  
                              Beta, Beta_wrt_minor, P_value, FDR, Indep_signal_i, 
                              N_samples, MAF, HWE_P_value, Imputation_score)]


### Combine them back together

all_mapped <- rbind(TSS_match, rescued2)
all_mapped
length(unique(all_mapped$hg38SNP_ID))
length(unique(with_ill_IDs$hg38SNP_ID))
# Still only loss of 5.

all_mapped$SNP_ENSG <- paste(all_mapped$hg38SNP_ID, all_mapped$ENSG_ID, sep = "_")
all_mapped$SNP_TSSdist <- abs(all_mapped$hg38SNP_pos - all_mapped$TSS)

### Edit April 2020: Write this table to explore contacting/noncontacting SNPs
write.table(all_mapped, file = "./descriptive_files/eQTLs_proxies_and_TSSs_MAF0.1_Rs0.9_including_proximal_hg38.txt", sep = "\t", quote = FALSE, row.names = FALSE)
###

proximal <- all_mapped[SNP_TSSdist<10000]
to_remove <- unique(proximal$SNP_ENSG)

distal_snps <- unique(all_mapped[!to_remove, on = "SNP_ENSG"])
length(unique(distal_snps$hg38SNP_ID))
# 3,224 remaining so far

# now remove proxies near TSS...
###### Check proxies distance to TSS

distal_snps$Proxy_TSSdist <- abs(distal_snps$TSS - distal_snps$hg38Proxy_pos)

distal_snps[is.na(Proxy_TSSdist)] # Should be empty
distal_snps$Proxy_ENSG <- paste(distal_snps$hg38Proxy_ID, distal_snps$ENSG_ID, sep = "_")

proxies_to_remove <- distal_snps[Proxy_TSSdist < 10000]
proxies_to_remove <- unique(as.character(proxies_to_remove$Proxy_ENSG))

distal_SNPs_and_proxies <- distal_snps[!proxies_to_remove, on="Proxy_ENSG"]
distal_SNPs_and_proxies
distal_SNPs_and_proxies[is.na(Proxy_TSSdist)] # Should still be empty

range(distal_SNPs_and_proxies$Proxy_TSSdist) 

length(unique(as.character(distal_SNPs_and_proxies$hg38SNP_ID)))
# Still 3,224 


distal_SNPs_and_proxies <- 
  distal_SNPs_and_proxies[, list(Chrom, hg19SNP_ID, hg38SNP_ID, hg38SNP_pos, SNP_MAF, hg19Proxy_ID, 
                                 hg38Proxy_ID, hg38Proxy_pos, Proxy_MAF, R2, ENSG_ID, ENST_ID, TSS, 
                                 rep_Gene, Gene, Array_Address_Id, Probe_Id,
                                 Probe_maps2TSS, target_TSS_source, SNP_TSSdist, Proxy_TSSdist,
                                 Beta, Beta_wrt_minor, P_value, FDR, Indep_signal_i, 
                                 N_samples, MAF, HWE_P_value, Imputation_score, SNP_ENSG, Proxy_ENSG)]


## Make file of distal SNPs, could be useful

write.table(distal_SNPs_and_proxies, file = "~/eCHiC/design/TSS_overlap/Probe_distal_eqtls_hg38.txt", 
            col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)


############### Assign to DpnII fragments and remove those not within 100bp of DpnII cut site


dpn <- fread("~/eCHiC/design/source/hg38_stuff/Digest_GRCh38_DpnII_None_11-24-01_12-06-2019.txt", header = TRUE, sep = "\t")

dpnII <- dpn[, list(Chromosome, Fragment_Start_Position, Fragment_End_Position, Fragment_Number)]
names(dpnII) = c("Chrom", "DpnII_start", "DpnII_end", "DpnII_index")

# Side note ---- check DpnII fragment length. Need to remove those with no 5' or 3' restriction site

dpn_temp <- dpn
names(dpn_temp)[6] = "five_prime_Restriction_Site"
names(dpn_temp)[7] = "three_prime_Restriction_Site"
dpn_temp2 <- dpn_temp[(five_prime_Restriction_Site == "Re1") & (three_prime_Restriction_Site == "Re1")]
dpn_temp2$length <- dpn_temp2$Fragment_End_Position - dpn_temp2$Fragment_Start_Position
range(dpn_temp2$length)
dpn_lengths <- dpn_temp2[length < 4000]
range(dpn_lengths$length)
hist(dpn_lengths$length)
median(dpn_lengths$length)
median(dpn_temp2$length)
test <- dpn_temp2[length > 1000000]

# 266 bp (including the really long ones; otherwise 265 bp)

#

dpnII$DpnII_index <- seq.int(nrow(dpnII)) ## Assigning new index as old was not unique
dpnII$Chrom = paste("chr", dpnII$Chrom, sep = "")

distal_SNPs_and_proxies$hg38Proxy_end <- distal_SNPs_and_proxies$hg38Proxy_pos

setkey(distal_SNPs_and_proxies, Chrom, hg38Proxy_pos, hg38Proxy_end)
proxies_assigned <- foverlaps(dpnII, distal_SNPs_and_proxies, by.x = c("Chrom", "DpnII_start", "DpnII_end"), 
                              type = "any", nomatch = NULL)

proxies_assigned$DpnII_proxy_dist1 <- abs(proxies_assigned$hg38Proxy_end - proxies_assigned$DpnII_start)
proxies_assigned$DpnII_proxy_dist2 <- abs(proxies_assigned$hg38Proxy_end - proxies_assigned$DpnII_end)

close_DpnII <- proxies_assigned[DpnII_proxy_dist1<=100 | DpnII_proxy_dist2<=100]

length(unique(as.character(close_DpnII$hg38SNP_ID)))

# 2,723 eQTLs remain; and 29,306 proxies


############# Regbuild overlap


################# NOW USING NEW ENSEMBL BUILD #######################

# Source file was downloaded from Ensembl96 GRCh38 on 01.07/2019. The new regbuild includes updates described on the 
# Ensembl blog, published in Jan 2019. The Regulatory Build now covers 21% of the genome (previous build from 2015 was 10%)

# This analysis will filter the proxies, and NOT eQTLs, on overlap with the regulatory build


regions <- fread("~/eCHiC/design/source/hg38_stuff/Ensembl96_GRCh38_RegBuild.txt", sep = "\t", 
                    header = TRUE)
names(regions) = c("Chrom", "RegBuild_start", "RegBuild_end", "Feature", "Feature_descrip")
regions$Chrom <- sub("^", "chr", regions$Chrom )

setwd("~/eCHiC/design/RegBuild96_overlap")

############# Perform overlaps for proxies and eQTLs

##### Add a new column to eqtls, that says if the "proxy" is also the "eQTL"


close_DpnII$hg38Proxy_ID <- as.character(close_DpnII$hg38Proxy_ID)
close_DpnII$hg38SNP_ID <- as.character(close_DpnII$hg38SNP_ID)
close_DpnII$Match <- close_DpnII$hg38Proxy_ID == close_DpnII$hg38SNP_ID
table(close_DpnII$Match)

##### Perform foverlaps (allowing for RegBuild = NA) and then filter on Match = TRUE or RegBuild!=NA

setkey(regions, Chrom, RegBuild_start, RegBuild_end)

prox_over <- foverlaps(close_DpnII, regions, by.x=c("Chrom", "hg38Proxy_pos", "hg38Proxy_end"), 
                       type = "any")

length(unique(as.character(prox_over$hg38SNP_ID)))
# 2,723


regbuild_overlap1 <- prox_over[Match == "FALSE" & !is.na(Feature)]
regbuild_overlap2 <- prox_over[Match == "TRUE"]

reg_filter <- rbind(regbuild_overlap1, regbuild_overlap2)

##############################################################

length(unique(as.character(reg_filter$hg38SNP_ID)))
length(unique(as.character(reg_filter$hg38Proxy_ID)))
length(unique(as.character(reg_filter$ENSG_ID)))

# Keep 2,259 SNPs and 7,147 proxies (which includes eQTLs)


################ Look to see what the regions are... ##################

for_plot <- unique(reg_filter[, .(hg38Proxy_ID, Feature)])
barplot(table(for_plot$Feature),
        col = c("blue", "orange", "yellow", "red", "tomato", "gold"),
        main = "Proxy overlap with RegBuild (Ensembl96)", 
        ylab = "Proxy frequency",
        ylim = c(0,3000))

## Make files

write.table(reg_filter, file = "~/eCHiC/design/RegBuild96_overlap/overlaps/reg96_filter_finalhg38.txt",
            row.names=FALSE,sep="\t", quote = FALSE)

jpeg(filename = "./plots/reg96_overlap.jpeg",
     width=800,height=500,res=70)

barplot(table(for_plot$Feature),
        col = c("blue", "orange", "yellow", "red", "tomato", "gold"),
        main = "Proxy overlap with RegBuild (Ensembl96)", 
        ylab = "Proxy frequency", 
        ylim = c(0, 2000))

dev.off()

length(unique(as.character(reg_filter$DpnII_index)))

# 6,451.

# How many fragments per eQTL?

t <- unique(reg_filter[, .(hg38SNP_ID, DpnII_index)])
tab <- table(as.character(t$hg38SNP_ID))
range(tab)
median(tab)
hist(tab, breaks = 61,
     xlab = "No. of DpnII fragments",
     ylab = "Frequency of eQTLs",
     main = "Bait fragments per eQTL for RegBuild96 filtered dataset")
abline(v=10, col = "blue")

# Remove regions that have > 10 DpnII fragments

tab2 <- as.data.table(tab)

tab3 <- tab2[N<=10]
length(tab3$V1)
tab4 <- tab2[N>10]
length(tab4$V1)
# 114 regions to remove; 2,145 remain

ta <- unique(as.character(tab3$V1))

reg_filter_final <- reg_filter[ta, on = "hg38SNP_ID"]

length(unique(as.character(reg_filter_final$hg38SNP_ID)))
# 2,145


# Try above graph of DpnII frags again

t <- unique(reg_filter_final[, .(hg38SNP_ID, DpnII_index)])
tab <- table(as.character(t$hg38SNP_ID))
range(tab)
median(tab)
hist(tab, breaks = 10,
     xlab = "No. of DpnII fragments",
     ylab = "Frequency of eQTLs",
     main = "Bait fragments per eQTL for RegBuild96 filtered dataset; N>10 removed")


## Find distribution of DpnII fragments; how many should be merged?

# RegBuild filter

gap0 <- reduce(IRanges(start=t$DpnII_index, end = t$DpnII_index), min.gapwidth =1)
gap0t <- as.data.table(gap0)
table(gap0t$width)
range0 <- gap0t$width

gap1 <- reduce(IRanges(start=t$DpnII_index, end = t$DpnII_index), min.gapwidth =2)
gap1t <- as.data.table(gap1)
table(gap1t$width)
range1 <- gap1t$width

gap2 <- reduce(IRanges(start=t$DpnII_index, end = t$DpnII_index), min.gapwidth =3)
gap2t <- as.data.table(gap2)
table(gap2t$width)
range2 <- gap2t$width

gap3 <- reduce(IRanges(start=t$DpnII_index, end = t$DpnII_index), min.gapwidth =4)
gap3t <- as.data.table(gap3)
table(gap3t$width)
range3 <- gap3t$width

gap4 <- reduce(IRanges(start=t$DpnII_index, end = t$DpnII_index), min.gapwidth =5)
gap4t <- as.data.table(gap4)
table(gap4t$width)
range4 <- gap4t$width


boxplot (range0, range1, range2, range3, range4, 
         names = c(0:4), 
         xlab = "Max no. of intervening DpnII fragments allowed in merge",
         ylab = "No. of consecutive DpnII fragments",
         main = "Distribution of merged fragments for RegBuild96-filtered dataset")

boxplot (range0, range1, range2, range3, range4, 
         names = c(0:4), 
         xlab = "Max no. of intervening DpnII fragments allowed in merge",
         ylab = "No. of consecutive DpnII fragments",
         main = "Distribution of merged fragments for RegBuild96-filtered dataset",
         outline = FALSE)


################# DO NOT SAVE THE FINAL EQTL FILE YET; NEED TO FILTER ON PROXY SET LENGTH



############## Get eTSS and mTSS targets ############################################################


### Define boundaries of proxy sets


# Get max and min of proxy locations, per eQTL/gene

filter_eqs <- unique(reg_filter_final[, list(Chrom, hg38SNP_ID, hg38SNP_pos, hg38Proxy_ID, hg38Proxy_pos, SNP_ENSG)])

min <- filter_eqs[filter_eqs[ , .I[which.min(hg38Proxy_pos)], by = c("Chrom", "hg38SNP_ID", "hg38SNP_pos", "SNP_ENSG")]$V1]
names(min) = c("Chrom", "hg38SNP_ID", "hg38SNP_pos", "Min_proxy_ID", "Min_proxy_pos", "SNP_ENSG")
min[is.na(Min_proxy_pos)] # Should be empty

max <- filter_eqs[filter_eqs[ , .I[which.max(hg38Proxy_pos)], by = c("Chrom", "hg38SNP_ID", "hg38SNP_pos", "SNP_ENSG")]$V1]
names(max) = c("Chrom", "hg38SNP_ID", "hg38SNP_pos", "Max_proxy_ID", "Max_proxy_pos", "SNP_ENSG")
max[is.na(Max_proxy_pos)] # Should be empty

comb <- min[max, on=c("Chrom", "hg38SNP_ID", "hg38SNP_pos", "SNP_ENSG")]
comb[is.na(Max_proxy_pos)]
comb[is.na(Max_proxy_ID)]
comb[is.na(Min_proxy_pos)]
comb[is.na(Min_proxy_ID)]

comb$Proxy_set_length <- comb$Max_proxy_pos - comb$Min_proxy_pos

range(comb$Proxy_set_length)
hist(comb$Proxy_set_length)

# Impose a cutoff

hist(comb$Proxy_set_length, 
     ylab = "eQTL frequency",
     xlab = "Proxy set length, bp",
     main = "Proxy set lengths")
abline(v=200000, col = "blue")

to_remove <- comb[Proxy_set_length>200000]
length(to_remove$hg38SNP_ID)
# 36 to remove

eq <- comb[Proxy_set_length<=200000]
length(unique(as.character(eq$hg38SNP_ID)))

# 2,109 eQTLs remain

## Now write final eQTL file.

### Remove those eQTLs that had proxy set length > 200 Kb: 31 to remove (object = to_remove)

final_filtered_eqtls <- reg_filter_final[!(reg_filter_final$hg38SNP_ID %in% to_remove$hg38SNP_ID),]
length(unique(final_filtered_eqtls$hg38SNP_ID))
# 2109

# No. of Proxies:
length(unique(as.character(final_filtered_eqtls$hg38Proxy_ID)))
# 5075

# No. of eGenes:
length(unique(as.character(final_filtered_eqtls$ENSG_ID)))
# 1834

# No. of DpnII fragments:
length(unique(final_filtered_eqtls$DpnII_index))
# 4651

write.table(final_filtered_eqtls, file = "~/eCHiC/design/final_design/final_filtered_eqtls_illumina_july19_hg38.txt", 
            quote = FALSE, sep = "\t", row.names = FALSE)


#### Using defined windows of 50 Kb surrounding eGenes and Mirror genes.

# Get back all illumina implicated TSS's to join to min_max sets - in order to define eTSS and mTSS windows

## Join final_filtered_eqtls to eq

with_sets <- final_filtered_eqtls[eq, on=c("Chrom", "hg38SNP_ID", "hg38SNP_pos", "SNP_ENSG")]

length(unique(with_sets$hg38SNP_ID))
with_sets[is.na(hg38Proxy_ID)]
with_sets[is.na(hg38SNP_ID)]
with_sets[is.na(TSS)]

with_sets$TSS <- as.numeric(with_sets$TSS)
with_sets$Min_proxy_pos <- as.numeric(with_sets$Min_proxy_pos)
with_sets$Max_proxy_pos <- as.numeric(with_sets$Max_proxy_pos)


with_sets2 <- with_sets %>% 
  mutate(is_TSS_in_proxy_window = TSS >= Min_proxy_pos & TSS <= Max_proxy_pos)

table(with_sets2$is_TSS_in_proxy_window)

# FALSE TRUE
# 18955 3342

with_sets2 <- as.data.table(with_sets2)
TSS_in_proxy_window <- with_sets2[is_TSS_in_proxy_window=="TRUE"]
eQTL_TSS_to_remove <- unique(TSS_in_proxy_window[, list(hg38SNP_ID, TSS)])  
length(unique(as.character(eQTL_TSS_to_remove$hg38SNP_ID)))

# 185

# Get the eqtl_TSS combos and remove from dataset (PICK THESE UP LATER)

TSS_outside_proxy_window <- with_sets2[is_TSS_in_proxy_window=="FALSE"]

median(TSS_outside_proxy_window$Proxy_set_length)
range(TSS_outside_proxy_window$Proxy_set_length)
hist(TSS_outside_proxy_window$Proxy_set_length)

# Max length = 199,091 bp 

### Get eTSS and mirror TSS windows using the Fixed Kb approach

filt <- TSS_outside_proxy_window

# Find if proxy set is upstream (TSS - proxyset = +ve) or downstream (TSS - proxyset = -ve); separate

filt$up_down <- filt$TSS - filt$Max_proxy_pos

upstream <- filt[up_down >0]
length(unique(as.character(upstream$hg38SNP_ID)))

downstream <- filt[up_down < 0]
length(unique(as.character(downstream$hg38SNP_ID)))

# Define X, Y and Z for upstream - use the max_proxy_loc as reference

up <- upstream
down <- downstream

X <- up$TSS - up$Max_proxy_pos
range(X)
up$X <- X
test <- up[X<10000]
test

up_distal <- up[X>=10000]
range(up_distal$X)

Y <- 5000
Z <- 45000

up_distal$Y <- Y
up_distal$Z <- Z

up_distal$opt_mTSS_pos <- up_distal$Min_proxy_pos - up_distal$X

up_distal$egene_window_start <- up_distal$TSS - Y
up_distal$egene_window_end <- up_distal$TSS + Z
up_distal$mirror_window_start <- up_distal$Min_proxy_pos - up_distal$X - Z
up_distal$mirror_window_end <- up_distal$Min_proxy_pos - up_distal$X + Y

range(up_distal$mirror_window_end) # Some are negative because they are too near the start of the chr. leave in for now

# Define X, Y and Z for downstream - use the min_proxy_pos as reference

X <- down$TSS - down$Min_proxy_pos

range(X)

down$X <- X
test <- down[X>-10000]
test



down_distal <- down[X<=-10000]
range(down_distal$X)

Y <- 5000
Z <- 45000

down_distal$Y <- Y
down_distal$Z <- Z

down_distal$opt_mTSS_pos <- down_distal$Max_proxy_pos - down_distal$X

down_distal$egene_window_start <- down_distal$TSS - Z
down_distal$egene_window_end <- down_distal$TSS + Y
down_distal$mirror_window_start <- down_distal$Max_proxy_pos - down_distal$X - Y
down_distal$mirror_window_end <- down_distal$Max_proxy_pos - down_distal$X + Z



# Combine upstream and downstream again

result <- rbind(up_distal, down_distal)

# Check window length

result$egene_window_length <- result$egene_window_end - result$egene_window_start
range(result$egene_window_length) # Should be 50 Kb

result$mirror_window_length <- result$mirror_window_end - result$mirror_window_start
range(result$mirror_window_length) # Should be 50 Kb

##############################################################################

######## See how many TSS are included in this design

all_TSS3 <- as.data.table(TSS)
all_TSS3$TSS_end <- all_TSS3$TSS 
all_TSS3$Chrom <- paste("chr", all_TSS3$Chrom, sep = "")
all_TSS4 <- all_TSS3[, list(Chrom, TSS, TSS_end, ENSG_ID, Gene, ENST_ID)]
names(all_TSS4) <- c("Chrom", "cap_TSS_start", "cap_TSS_end", "cap_ENSG_ID", "cap_Gene", "cap_ENST_ID")
all_TSS4

### Overlap for the eGene window

setkey(result, Chrom, egene_window_start, egene_window_end)
eGenes <- foverlaps(all_TSS4, result, by.x=c("Chrom", "cap_TSS_start", "cap_TSS_end"), 
                    type = "any", nomatch = NULL)

### Overlap for the mirror window

setkey(result, Chrom, mirror_window_start, mirror_window_end)
mirrors <- foverlaps(all_TSS4, result, by.x=c("Chrom", "cap_TSS_start", "cap_TSS_end"), 
                     type = "any", nomatch = NULL)

### Find out how many genes/TSS in total to target

eGenes$captured_gene_type <- "in_eGene_window"
length(unique(eGenes$cap_ENSG_ID))
length(unique(eGenes$cap_TSS_start))
# 5680 eGenes, 24,772 TSS


mirrors$captured_gene_type <- "in_mGene_window"
length(unique(mirrors$cap_ENSG_ID))
length(unique(mirrors$cap_TSS_start))
# 4616 mirror genes, 15,014 TSS

all_TSS_in_windows <- unique(rbind(eGenes, mirrors))
length(unique(all_TSS_in_windows$cap_ENSG_ID))
length(unique(all_TSS_in_windows$cap_TSS_start))
# Total = 9182 genes, 35,926 TSS

# How many DpnII fragments?

dpnII_for_genes <- dpnII
names(dpnII_for_genes) = c("Chrom", "gen_DpnII_start", "gen_DpnII_end", "gen_DpnII_index")

setkey(all_TSS_in_windows, Chrom, cap_TSS_start, cap_TSS_end)
genes_dpnII <- foverlaps(dpnII_for_genes, all_TSS_in_windows, by.x=c("Chrom", "gen_DpnII_start", "gen_DpnII_end"),
                         type="any", nomatch = NULL)

length(unique(genes_dpnII$gen_DpnII_index))

# 18,392 DpnII fragments for the genes outside of the proxy sets

# Now need to add in those regions where the eTSS falls within the proxy set

TSS_in_proxy_window


## 1. Get all TSS in these windows (foverlaps on min proxy - max proxy with TSS file)
## 2. Get all proxies in these windows (just based on the SNP_ID)
## 3. Get distance of every single TSS in window to every single proxy in window
## 4. If a TSS is <10Kb of any proxy, remove it (currency = "chr_TSS"; even if it is distal from some SNPs you want to remove it)

# 1. Get TSS
setkey(TSS_in_proxy_window, Chrom, Min_proxy_pos, Max_proxy_pos)
targets <- unique(foverlaps(all_TSS4, TSS_in_proxy_window, by.x=c("Chrom", "cap_TSS_start", "cap_TSS_end"),
                            type="any", nomatch = NULL))

# 2. Get proxies


final_filtered_eqtls2 <- unique(final_filtered_eqtls[, list(Chrom, hg38SNP_ID, hg38SNP_pos, SNP_MAF, hg38Proxy_ID, 
                                                            hg38Proxy_pos, Proxy_MAF, R2)])


targets2 <- unique(targets[, list(Chrom, hg38SNP_ID, hg38SNP_pos, SNP_MAF, Min_proxy_ID, Min_proxy_pos, Max_proxy_ID, Max_proxy_pos, 
                                  Proxy_set_length, cap_TSS_start, cap_TSS_end, cap_ENSG_ID, cap_Gene, cap_ENST_ID)])

targets3 <- targets2[final_filtered_eqtls2, on =c("Chrom", "hg38SNP_ID", "hg38SNP_pos", "SNP_MAF"), allow.cartesian = TRUE, nomatch = NULL]


# 3. Get proxy_cTSS distances

targets3$cap_TSS_proxy_dist <- abs(targets3$cap_TSS_start - targets3$hg38Proxy_pos)
targets3$cap_Chr_TSS <- paste(targets3$Chrom, targets3$cap_TSS_start, sep = "_")

# 4.

close <- targets3[cap_TSS_proxy_dist < 10000]

length(unique(as.character(close$cap_Chr_TSS)))
# 1132 cap TSS to remove


targets4 <- targets3[!(targets3$cap_Chr_TSS %in% close$cap_Chr_TSS),]


range(targets4$cap_TSS_proxy_dist)
# This is fine: >10Kb

length(unique(as.character(targets4$cap_Chr_TSS)))
# 2,125 TSS


### How many DpnII fragments?


setkey(targets4, Chrom, cap_TSS_start, cap_TSS_end)
TSS_in_windows_assigned <- foverlaps(dpnII_for_genes, targets4, by.x=c("Chrom", "gen_DpnII_start", "gen_DpnII_end"),
                                     type="any", nomatch = NULL)


length(unique(TSS_in_windows_assigned$gen_DpnII_index))

# 934 DpnII fragments

############### COMBINE ALL GENE TSS DPNII FRAGS ############

without <- genes_dpnII[, list(Chrom, cap_TSS_start, cap_TSS_end, cap_ENSG_ID, cap_Gene, 
                              cap_ENST_ID, gen_DpnII_start, gen_DpnII_end, gen_DpnII_index)]

within <- TSS_in_windows_assigned[, list(Chrom, cap_TSS_start, cap_TSS_end, cap_ENSG_ID, cap_Gene, 
                                         cap_ENST_ID, gen_DpnII_start, gen_DpnII_end, gen_DpnII_index)]


all_gene_dpnII <- unique(rbind(without, within))
length(unique(all_gene_dpnII$gen_DpnII_index))

# 19,089 fragments


############### COMBINE WITH EQTL AND PROXY DPNII

# use final_filtered_eqtls

all_proxies_dpnII_frags <- unique(final_filtered_eqtls[, list(Chrom, DpnII_start, DpnII_end, DpnII_index)])
all_gene_dpnII_frags <- unique(all_gene_dpnII[, list(Chrom, gen_DpnII_start, gen_DpnII_end, gen_DpnII_index)])
names(all_gene_dpnII_frags) = c("Chrom", "DpnII_start", "DpnII_end", "DpnII_index")

capture_design_dpnII_frags <- unique(rbind(all_proxies_dpnII_frags, all_gene_dpnII_frags))

capture_design_dpnII_frags


######################## 23,667 DpnII fragments in total ... 

## Make nice file for gene targets

in_proxy_windows <- TSS_in_windows_assigned[, list(Chrom, hg38SNP_ID, hg38SNP_pos, SNP_MAF, hg38Proxy_ID, hg38Proxy_pos, Proxy_MAF, 
                                                   R2, Min_proxy_ID, Min_proxy_pos, Max_proxy_ID, Max_proxy_pos, 
                                                   Proxy_set_length, cap_TSS_start, cap_TSS_end, cap_ENSG_ID, 
                                                   cap_Gene, cap_ENST_ID, gen_DpnII_start, gen_DpnII_end, gen_DpnII_index)]
in_proxy_windows$captured_gene_type <- "In_proxy_window"

genes_dpnII2 <- genes_dpnII[, list(Chrom, hg38SNP_ID, hg38SNP_pos, SNP_MAF, hg38Proxy_ID, hg38Proxy_pos, Proxy_MAF, 
                                   R2, Min_proxy_ID, Min_proxy_pos, Max_proxy_ID, Max_proxy_pos, 
                                   Proxy_set_length, cap_TSS_start, cap_TSS_end, cap_ENSG_ID, 
                                   cap_Gene, cap_ENST_ID, gen_DpnII_start, gen_DpnII_end, gen_DpnII_index, captured_gene_type)]

all_captured_genes <- unique(rbind(in_proxy_windows, genes_dpnII2))

length(unique(as.character(all_captured_genes$gen_DpnII_index)))

write.table(all_captured_genes, file="~/eCHiC/design/final_design/all_captured_genes_finalhg38.txt", 
            col.names = TRUE, sep = "\t", row.names = FALSE, quote = FALSE)


########################

# Write files to check on WashU

# DpnII

write.table(capture_design_dpnII_frags, file = "~/eCHiC/design/final_design/dpnIIhg38.bed", 
            col.names = FALSE, sep = "\t", row.names = FALSE, quote = FALSE)

# eQTLs

temp <- final_filtered_eqtls
temp$SNP_Gene <- paste(temp$hg38SNP_ID, temp$Gene, sep = "_")
eqtls.bed <- unique(temp[, list(Chrom, hg38SNP_pos, SNP_end=(hg38SNP_pos+1), SNP_Gene)])

write.table(eqtls.bed, file = "~/eCHiC/design/final_design/eqtlshg38.bed", 
            col.names = FALSE, sep = "\t", row.names = FALSE, quote = FALSE)

# Proxies

proxies.bed <- unique(temp[, list(Chrom, hg38Proxy_pos, Proxy_end=(hg38Proxy_pos+1), hg38SNP_ID)])
proxies.bed # 5,388

write.table(proxies.bed, file = "~/eCHiC/design/final_design/proxieshg38.bed", 
            col.names = FALSE, sep = "\t", row.names = FALSE, quote = FALSE)


# Targeted TSS


TSS.bed <- unique(all_captured_genes[, list(Chrom, cap_TSS_start, cap_TSS_end=(cap_TSS_end+1), cap_Gene)])

write.table(TSS.bed, file = "~/eCHiC/design/final_design/TSShg38.bed", 
            col.names = FALSE, sep = "\t", row.names = FALSE, quote = FALSE)

# Proxy sets - note- will add 1 to the "end" point to make it compatible for WashU

proxy_sets.bed <- unique(eq[, list(Chrom, Min_proxy_pos, Max_proxy_pos, hg38SNP_ID)])
proxy_sets.bed$Max_proxy_pos = proxy_sets.bed$Max_proxy_pos +1

write.table(proxy_sets.bed, file = "~/eCHiC/design/final_design/proxy_setshg38.bed", 
            col.names = FALSE, sep = "\t", row.names = FALSE, quote = FALSE)

# eGene windows

ewindows.bed <- unique(result[, list(Chrom, egene_window_start, egene_window_end, hg38SNP_ID)])

write.table(ewindows.bed, file = "~/eCHiC/design/final_design/ewindowshg38.bed", 
            col.names = FALSE, sep = "\t", row.names = FALSE, quote = FALSE)

# mGene windows

mwindows.bed <- unique(result[, list(Chrom, mirror_window_start, mirror_window_end, hg38SNP_ID)])
mwindows.bed <- mwindows.bed[mirror_window_start >0 & mirror_window_end >0]

write.table(mwindows.bed, file = "~/eCHiC/design/final_design/mwindowshg38.bed", 
            col.names = FALSE, sep = "\t", row.names = FALSE, quote = FALSE)


# RegBuild
reg_regions <- unique(temp[, list(Chrom, RegBuild_start, RegBuild_end, Feature)])
reg_regions.bed <- reg_regions[!is.na(RegBuild_start)]
reg_regions.bed$Feature <- gsub(" ", "_", reg_regions.bed$Feature)

write.table(reg_regions.bed, file = "~/eCHiC/design/final_design/reg_regionshg38.bed", 
            col.names = FALSE, sep = "\t", row.names = FALSE, quote = FALSE)







######################## END ###############################

