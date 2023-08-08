### Load in the allelic counts from phASER, which were generated per gene interaction
### Note that the allelic counts file seem to ignore the blacklisting option and blacklisting is only in the haplo counts files
### We can blacklist manually if we want to
### Per sample, load the AS counts and the per-eGene-interaction counts.

library(data.table)
library(dplyr)
library(tidyr)

#### Parse the output from phASER
setwd("~/HRJ_monocytes/AS_CHiC/BaseQTL/phASER_out")

# This file details which Genes are called in which file.
# Will need to match the Gene no. (from the name of the file) with the SNP ID (hg38) and then match it to this file.
Genes <- fread("~/HRJ_monocytes/AS_CHiC/BaseQTL/scripts/bedpe_for_chic_egenes_with_Gene_number.txt")
names(Genes) = c("SNPChr", "SNP_Dpnstart", "SNP_Dpnend", "TSSChr", "hg38TSS_5Kstart", "hg38TSS_5Kend", "ID", "Gene_number")
# Will need the hg38 proxy IDs, which are here:
all_info <- fread("~/HRJ_monocytes/eqtls/snps_went_into_loopingQTL_analysis_with_gene_targets.txt")

### Make a file that has hg19Proxy_ID, hg38Proxy_ID, ENSG, Gene and Gene_number
Genes_short <- unique(Genes[, .(ID, Gene_number)])
Genes_short[, c("DpnID", "TSSID", "hg19Proxy_ID", "ENSG_ID", "Gene") := tstrsplit(ID, split = "_")]
Genes_short[, c("DpnID", "TSSID", "ID") := NULL]
Genes_short <- unique(Genes_short)
all_info_short <- unique(all_info[, .(hg19Proxy_ID, hg38Proxy_ID, hg38Proxy_pos, Chr, ENSG_ID, Gene)])

Gene_info <- Genes_short[all_info_short, on = c("hg19Proxy_ID", "ENSG_ID", "Gene"), nomatch = NULL]

##############################################################################

# Make individual files and one large matrix per analysis

###### Testing ######
#mycounts <- fread("./S025NM-01_Gene1_noBlacklist.allelic_counts.txt")
#Gene_info <- Gene_info
#Gene_number <- 1
#####################

get_counts <- function(counts_list, Gene_info, Gene_number) { # Gene_number should be in format "GeneN"
  all_counts <- data.table()
  for(file in counts_list) {
    mycounts <- fread(file, sep = "\t")
    names(mycounts)[1] = "Chr"
    names(mycounts)[2] = "hg38Proxy_pos"
    id <- paste("_", Gene_number, "_noBlacklist.allelic_counts.txt", sep = "")
    myname <- sub(id, "", file)
    
    ## Limit the Gene_info file to just the gene that we are interested in
    
    our_gene_number <- sub("Gene", "Gene_", Gene_number)
    our_gene <- Gene_info[Gene_number == our_gene_number]

    ######## Get the hg19 SNP ID and the Gene info, intersecting based on chr:hg38pos
    with_info <- mycounts[our_gene, on = c("Chr", "hg38Proxy_pos"), nomatch = NULL]
    #with_info[which(duplicated(with_info[, hg19Proxy_ID])), ] Now empty

    with_info[, Sample := myname]
    
    all_counts <- rbind(all_counts, with_info, fill= TRUE)
  }
  return(all_counts)
}

all_genes <- data.table()
for(gene in 1:33) {
        Gene_counts <- list.files("./", pattern = paste0("\\_Gene", gene, "_noBlacklist.allelic_counts.txt$"))
	Genename <- paste0("Gene", gene)
	Gene_results <- get_counts(Gene_counts, Gene_info, Genename)
	all_genes <- rbind(all_genes, Gene_results)
}


# Should be in the format: feature, SNP, hg38Pos, allele, counts

###### To test with Dplyr

matrix1 <- as.data.table(all_genes %>% pivot_longer (
  cols = c("refAllele", "altAllele"), 
  names_to = "Allele_type",
  values_to = "Allele",
))

matrix1[Allele_type == "refAllele", Counts := refCount]
matrix1[Allele_type == "altAllele", Counts := altCount]
matrix1[, feature := paste(hg19Proxy_ID, ENSG_ID, Gene, sep = "_")]
matrix2 <- unique(matrix1[, .(feature, hg19Proxy_ID, hg38Proxy_pos, Allele, Counts, Sample)]) # check if this matters that we didn't filter to these cols?

matrix2 <- as.data.table(matrix1 %>% pivot_wider (
  id_cols = c("feature", "hg19Proxy_ID", "hg38Proxy_pos", "Allele"), 
  names_from = Sample,
  values_from = Counts,
  values_fill = 0
))

setnames(matrix2, c("hg19Proxy_ID", "hg38Proxy_pos", "Allele"), c("SNP", "hg38SNP_pos", "allele")) ## matches the last AS_ATAC format for Elena


####### How many can we use??
egenes_matrix2 <- copy(matrix2)
egenes_matrix2[, rs := rowSums(.SD), .SDcols = 5:38]
test <- egenes_matrix2[rs >= 50]
quantile(egenes_matrix2$rs)
test$hg19Proxy_ID <- gsub("\\_[A,C,G,T]", "", x = test$SNP)
length(unique(test$hg19Proxy_ID)) # 261 (previously 794 proxy snps) BUT we do not include homozygotes!
test2 <- egenes_matrix2[rs == 0]
test2$hg19Proxy_ID <- gsub("\\_[A,C,G,T]", "", x = test2$SNP)
length(unique(test2$hg19Proxy_ID)) # 3567 (previosuly 3307 proxy snps)
########

######## Write files.

setorder(matrix2, feature)

setwd("~/HRJ_monocytes/AS_CHiC/BaseQTL/allele_counts")
write.table(matrix2, file = "mono_34reps_CHiC_allele_counts_phASER.txt", 
            quote = FALSE, row.names = FALSE, sep = "\t")

############# phASER does not handle homozygotes so need to think about how to represent them in the table.
############# 


################################################################################################
############### Now to get the TOTAL counts per DpnII fragment for those interactions of interest.
################################################################################################

# To get total counts, can use the bedpe files generated when searching for eQTL-eGene interactions.

setwd("~/HRJ_monocytes/AS_CHiC/BaseQTL/bedpe")


###### Testing ######
#mycounts <- fread("./S025NM-01_eQTLs_to_eGenes.Gene_1.bedpe")
#Gene_info <- Gene_info
#Gene_number <- "Gene1"
#####################

# The bedpe has the info that we need in it (SNP, Gene, etc)

collate_total_counts <- function(filelist, Gene_info, Gene_number) { # Gene_number should be in format "GeneN"
  all_counts <- data.table()
  for(file in filelist) {
    mycounts <- fread(file)
    setnames(mycounts, c("V1", "V7", "V17"), c("Chr", "readID", "ID"))
    
    small <- unique(mycounts[, .(Chr, readID, ID)])
    
    small[, c("DpnID", "TSSID", "hg19Proxy_ID", "ENSG_ID", "Gene") := tstrsplit(ID, split = "_")]
    with_hg38 <- small[Gene_info, on = c("hg19Proxy_ID", "ENSG_ID", "Gene", "Chr"), nomatch = NULL] # We do not lose any lines
    with_hg38[, c("DpnID", "TSSID", "hg38Proxy_ID", "Gene_number") := NULL]
    
    with_hg38[, feature := paste(hg19Proxy_ID, ENSG_ID, Gene, sep = "_")]
    to_tally <- unique(with_hg38[, .(feature, hg19Proxy_ID, hg38Proxy_pos, readID)])

    to_tally[, `:=` (count = .N), by = feature]
    
    total <- unique(to_tally[, .(feature, hg19Proxy_ID, hg38Proxy_pos, count)])
    
    mygene <- sub("Gene", "Gene_", Gene_number)
    id <- paste("_eQTLs_to_eGenes.", mygene, ".bedpe", sep = "")
    myname <- sub(id, "", file)
    total$Sample <- myname

    all_counts <- rbind(all_counts, total, fill= TRUE)
  }
  return(all_counts)
}

all_totals <- data.table()
for(gene in 1:33) {
  Gene_counts <- list.files("./", pattern = paste0("\\_eQTLs_to_eGenes.Gene_", gene, ".bedpe$"))
  Genename <- paste0("Gene", gene)
  Gene_results <- collate_total_counts(Gene_counts, Gene_info, Genename)
  all_totals <- rbind(all_totals, Gene_results)
}
fwrite(all_totals, file = "all_totals_temp.txt", 
       sep = "\t", quote = F, row.names = F, col.names = T)

####### Make the large matrix

totals_mat <- as.data.table(all_totals %>% pivot_wider (
  id_cols = c("feature", "hg19Proxy_ID", "hg38Proxy_pos"), 
  names_from = Sample,
  values_from = count,
  values_fill = 0
))

setnames(totals_mat, c("hg19Proxy_ID", "hg38Proxy_pos"), c("SNP", "hg38SNP_pos")) ### 
setorder(totals_mat, feature)


####### How many can we use??
totals_mat2 <- copy(totals_mat)
totals_mat2[, rs := rowSums(.SD), .SDcols = 4:37]
test <- totals_mat2[rs >= 50]
quantile(totals_mat2$rs)
length(unique(test$SNP)) # 5118 SNPs, previously 5196 
########


######## Write file

setwd("~/HRJ_monocytes/AS_CHiC/BaseQTL/allele_counts")
write.table(totals_mat, file = "mono_34reps_CHiC_total_counts_from_bedpe.txt", 
            quote = FALSE, row.names = FALSE, sep = "\t")

