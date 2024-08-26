########## Modifying the homer scripts to find cases where snps create new DpnII sites.
########## For every SNP in the genome, see if it (alone or in combination with a SNP within 4bp) creates a DpnII cut site.
########## Have to make two fasta files: one for "REF", one for "ALT".
########## Make for all genotyped SNPs and then probe using the sig results from Leo - finding SNPs in LD which are in same original frag.

#### Useful: the Biostrings Vignette, which is here https://bioconductor.org/packages/devel/bioc/vignettes/Biostrings/inst/doc/BiostringsQuickOverview.pdf

library(Biostrings)
library(data.table)
setwd("~/HRJ_monocytes/findmotifs")

##### I previously already made 40bp around all genotyped snps. See below for what was done.

##### 1. Read in the Genome as fasta files (DNA) 
#chr22 <- readDNAStringSet("~/HRJ_monocytes/external_data/GRCh38/chr22.fa")

##### 2. Get 40bp regions around each SNP.
options(scipen=999)
#genos <- fread("~/HRJ_monocytes/AS_ATAC/BaseQTL/snps/all.rsq3.maf05.hg38.dose.finalRenamed.newIDs.sites.vcf")
#genos22 <- genos[`#CHROM` == 22 & POS < 16851356]
#genos22_list <- c(genos22$POS-20)
#genos22_for_range <- lapply(genos22_list, IRanges, width = 41) # this combined with prev (pos-20) gets us the coords we need
#at22 <- as.list(genos22_for_range)
#unlist(extractAt(chr22, at22)) # nope

##### OK, make bed files per chrom that can be used with bedtools for getfasta.
##### Bed file should be 0-based (POS here is 1-based)
#setwd("~/HRJ_monocytes/findmotifs")
#dir.create("snp_bed")
#setwd("snp_bed")
#genos2 <- copy(genos)
#genos2[, pos0 := POS-21]
#genos2[, pos1 := POS+20]
#genos_bed <- unique(genos2[, .(`#CHROM`, pos0, pos1, ID)])
#names(genos_bed) = c("Chr", "start", "stop", "ID")
#chroms <- seq(1:22)
#for(chr in chroms) {
#  mychrom <- genos_bed[Chr == chr]
#  myname=paste0("all.rsq3.maf05.hg38.", chr, ".40bp.bed")
#  fwrite(mychrom, file = myname, sep = "\t", quote = F, row.names = F, col.names = F)
#}

##### 3. Now running bedtools

##### 4. Read in the results and change the snp, per significant dataset.

### btw, check that the snp matches the reference genome. Yes, seems to.
setwd("~/HRJ_monocytes/findmotifs/snp_fasta")

chroms <- seq(1:22)

myseq <- data.table()
for(i in chroms) {
  myname <- paste0("all.rsq3.maf05.hg38.", i, ".40bp.txt")
  mychrom <- fread(myname, header = FALSE)
  myseq <- rbind(myseq, mychrom)
}

#myseq <- fread("./all.rsq3.maf05.hg38.1.40bp.txt")
names(myseq) = c("SNP_loc", "Sequence")

### Split the "SNP_loc" column, and make the allele. We should do REF as well as ALT just to make sure.
myseq[, c("SNP", "hg38loc") := tstrsplit(SNP_loc, split = "::")]
myseq[, SNP_loc := NULL]
myseq[, c("REF", "ALT") := tstrsplit(SNP, split = ":", keep = c(3,4))]

### Get the location of the SNP itself.
myseq[, c("Chr", "coords") := tstrsplit(hg38loc, split = ":")]
myseq[, c("seq_start", "seq_stop") := tstrsplit(coords, split = "-")]
myseq[, seq_start := as.numeric(seq_start)]
myseq[, seq_stop := as.numeric(seq_stop)]
#### Get the SNP location in hg38 (NOTE here I am not talking about proxies. Talking about all genotyped SNPs.)
myseq[, hg38SNP_loc := seq_start + 21]
#### Save this table as the REF.
fwrite(myseq, file = "./ALL_genotyped_SNPs_REF_40bp_fasta_with_locations.txt", sep = "\t", 
       quote = F, row.names = F, col.names = T)
myseq <- fread("./ALL_genotyped_SNPs_REF_40bp_fasta_with_locations.txt")

######## Now, mutate the SNPs within the same 40bp sequence and see if a new DpnII site is generated.
######## First off, find the ones that ALREADY have a DpnII cut site that intersect the SNP.
######## The cut site needs to be: GATC and the SNP is at position 21
myseq[, seq1 := substr(Sequence, 18, 21)]
myseq[, seq2 := substr(Sequence, 19, 22)]
myseq[, seq3 := substr(Sequence, 20, 23)]
myseq[, seq4 := substr(Sequence, 21, 24)]
### Any of these four seqs, if containing GATC, have the ref SNP in the cut site
dpnII_REF <- myseq[seq1 == "GATC" | seq2 == "GATC" | seq3 == "GATC" | seq4 == "GATC"]
fwrite(dpnII_REF, file = "~/HRJ_monocytes/findmotifs/find_new_dpnII/all_geno_SNP_DpnII_REF.txt", 
       sep = "\t", quote = F, row.names = F, col.names = T)

######################### Function to overlap multiple snps into a seq
#### Now re-intersect to get the SNP overlap when they share a sequence
#### First intersect the myseq table with your chosen snps, before running. e.g. seqtable[SNP %in% eqs$snp_id]
#### Creates all possible shared seqs for snps in the list...
#### Doesnt work unless you group by hg19SNP_ID first! Now incorporating into the function below.

#get_seq_overlaps <- function(myseq) {
#  snps[, hg38Proxy_pos2 := hg38Proxy_pos]
#  snps <- unique(myseq[, .(SNP, REF, ALT, Chr, hg38Proxy_pos)])
#  locations <- unique(myseq[, .(Sequence, hg38loc, Chr, seq_start, seq_stop)])
#  snp_overlap <- foverlaps(snps, locations, by.x = c("Chr", "hg38Proxy_pos", "hg38Proxy_pos2"), nomatch = NULL)
#  snp_overlap[, seqPos := hg38Proxy_pos - seq_start + 1]
#  setkey(locations, Chr, seq_start, seq_stop)
#  seqtable <- unique(snp_overlap[, .(Sequence, SNP, REF, ALT, seqPos)])
#  return(seqtable)
#}

#####################################


##### 1. Running the function for SNP sets
leo_sets <- readRDS("~/HRJ_monocytes/leo_triplets/input/alex_inputs/Helen_geno_matrix_list.Rds")
mynames1 <- seq(1:8972)
mynames2 <- rep("set_", 8972)
mynames <- as.list(paste0(mynames2, mynames1))
names(leo_sets) = mynames
#


my_table <- data.table()

for( i in seq_along(leo_sets)){
  if(nrow(leo_sets[[i]]) > 0) { ### because we have an empty table in the list, no. 2273
    my_set <- as.data.table(leo_sets[[i]], keep.rownames = T)
    my_set[, SNP := paste(rn, i, sep = "_")]
    my_set2 <- my_set[, .(SNP)]
  }
  my_table <- rbind(my_table, my_set2)
}



#columnselect<-function(dF){
#  df[,c("SNP")]
#}
#mylist <- lapply(X=leo_sets, FUN=columnselect)



### Get the hg38 locations and remove the sample IDs
#sets <- do.call(rbind.data.frame, c(leo_sets)) # unfortunately this doesnt preserve the snp name when there is only one snp in the set

my_table[, c("hg19_ID", "set") := tstrsplit(SNP, split = "_", fixed = TRUE, fill = NA)] 
fwrite(my_table, file = "~/HRJ_monocytes/leo_triplets/input/alex_inputs/Helen_SNP_sets.txt", 
       sep = "\t", quote = F, row.names = F, col.names = T)

# how many snps per set?
mytable <- unique(my_table)
mytabletally <- as.data.table(mytable %>% group_by(set) %>% tally())
range(mytabletally$n) # 1 to 197
median(mytabletally$n) # 27
hist(mytabletally$n)

### Get the hg38 locations
locs <- as.data.table(read.csv("~/HRJ_monocytes/leo_triplets/input/alex_inputs/snp_chromlocations.csv"))
setkey(locs, annot)
sets_locs <- my_table[locs, on = c(hg19_ID = "annot"), nomatch = NULL]
sets_locs[, Chr := as.numeric(substr(chrom, 4,5))]

# will group by set in the mutate seqs function
# get the REF sequences from above.
myseq <- fread("~/HRJ_monocytes/findmotifs/snp_fasta/ALL_genotyped_SNPs_REF_40bp_fasta_with_locations.txt")
setkey(myseq, SNP, Chr, hg38SNP_loc)
sets_seqs1 <- sets_locs[myseq, on = c(hg19_ID = "SNP", "Chr", position = "hg38SNP_loc"), nomatch = NULL]


## Change the genotypes if the SNPs are in the same SNP set. Because we don't know about LD. But surely they will be in LD if that close.
## Now requires that the alleles are listed as REF and ALT
mutate_seqs_nt <- function(myseq) {
  final_refalt <- data.table()
  set_list <- unique(as.list(myseq$set))
  for(SET in set_list) { # Split by the SNP set
    my_set <- myseq[set == SET]
    snps <- unique(my_set[, .(hg19_ID, REF, ALT, Chr, position, set)])
    snps[, position2 := position]
    locations <- unique(my_set[, .(Sequence, position, Chr, seq_start, seq_stop)])
    locations[, mod_seq_start := seq_start + 1]
    setkey(locations, Chr, mod_seq_start, seq_stop)
    
    ### The following overlap will get all the snps in LD within the 40bp fasta region
    snp_overlap <- foverlaps(snps, locations, by.x = c("Chr", "position", "position2"), nomatch = NULL, type = "within")
    snp_overlap[, seqPos := i.position - seq_start] # used to be seq_start + 1 but I had the positions one off.
    # Used i.position bc this is the position of the SNP in question (whereas position is just the middle SNP.)
    seqtable <- unique(snp_overlap[, .(Sequence, hg19_ID, REF, ALT, seqPos, set)])
    
    
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
      mysequence_id <- paste(c(mytable$hg19_ID), collapse = "_") # per sequence, get all the snp IDs pasted
      
      mysequence_REF_char <- as.character(mysequence_REF)
      mysequence_ALT_char <- as.character(mysequence_ALT)
      
      current_set <- as.character(unique(mytable$set))
      
      # get a column with all the snps in that particular sequence, and keep the set ID
      seqs <- data.table(REFSeq = mysequence_REF_char, ALTSeq = mysequence_ALT_char, id = mysequence_id, set = current_set)
      refalt <- rbind(refalt, seqs, fill = TRUE)
    }
    #return(refalt)
    final_refalt <- unique(rbind(final_refalt, refalt, fill = TRUE))
  }
  return(final_refalt)
}


#### Run the mutate function - testing on two sets.
sets_seqs <- unique(sets_seqs1[, .(set, hg19_ID, Chr, position, REF, ALT, Sequence, seq_start, seq_stop)])
#sets_seqs_tester <- sets_seqs[set == "set_1" | set == "set_2"]
#eqs_mutated <- mutate_seqs_nt(sets_seqs_tester) 
#sets_seqs_tester <- sets_seqs[set == "155"]
#eqs_mutated <- mutate_seqs_nt(sets_seqs_tester)
#eqs_mutated[id %like% "10:17004814:G:A"]

eqs_mutated <- mutate_seqs_nt(sets_seqs)


#### The above works
setwd("~/HRJ_monocytes/findmotifs/find_new_dpnII")
fwrite(eqs_mutated, file = "LeoSNPsets.REFALTseqs.txt", 
       row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

#### Then need to look for DNPII sites.
new_seqs <- fread("~/HRJ_monocytes/findmotifs/find_new_dpnII/LeoSNPsets.REFALTseqs.txt")
## In fact here, look in both the REF and the ALT.
new_seqs[, REFseq1 := substr(REFSeq, 18, 21)]
new_seqs[, REFseq2 := substr(REFSeq, 19, 22)]
new_seqs[, REFseq3 := substr(REFSeq, 20, 23)]
new_seqs[, REFseq4 := substr(REFSeq, 21, 24)]
### Any of these four seqs, if containing GATC, have the ref SNP in the cut site for the reference allele
new_seqs[REFseq1 == "GATC" | REFseq2 == "GATC" | REFseq3 == "GATC" | REFseq4 == "GATC", Dpn_site_REF := TRUE]
new_seqs[, ALTseq1 := substr(ALTSeq, 18, 21)]
new_seqs[, ALTseq2 := substr(ALTSeq, 19, 22)]
new_seqs[, ALTseq3 := substr(ALTSeq, 20, 23)]
new_seqs[, ALTseq4 := substr(ALTSeq, 21, 24)]
### Any of these four seqs, if containing GATC, have the ref SNP in the cut site for the alternative allele
new_seqs[ALTseq1 == "GATC" | ALTseq2 == "GATC" | ALTseq3 == "GATC" | ALTseq4 == "GATC", Dpn_site_ALT := TRUE]

fwrite(new_seqs, file = "~/HRJ_monocytes/findmotifs/find_new_dpnII/LeoSNPsets_DpnII_REFALT.txt", 
       sep = "\t", quote = F, row.names = F, col.names = T)





