# Get LD proxy SNPs for cQTLs
# $ plink2 --pfile  ~/analysis/plink2/all_hg38 'vzs' --r2-unphased 'allow-ambiguous-allele' --ld-window-kb 50 --ld-window-r2 0.8 --ld-snp-list GUESS_BaseQTL_contactQTLs_hg19_rsids.txt  --out consensus_set_LD_50k_revision_plink2 --threads 8 

library(data.table)
setDTthreads(2)
ld = fread("consensus_set_LD_50k_revision_plink2.vcor")
# Note plink2 doesn't include input SNPs without SNPs in LD anymore
# Also unlike in plink v1.x ID_A is never duplicated in ID_B, so need to add them explicitly
length(unique(ld$ID_A)) #[1] 591
length(unique(ld$ID_B)) #[1] 5148
samesnpdt = data.table(CHROM_A=rep(NA_real_,length(unique(ld$ID_A))), POS_A=rep(NA_real_,length(unique(ld$ID_A))), ID_A=unique(ld$ID_A), CHROM_B=rep(NA_real_,length(unique(ld$ID_A))), POS_B=rep(NA_real_,length(unique(ld$ID_A))), ID_B=unique(ld$ID_A), UNPHASED_R2=1)
setnames(samesnpdt, 1, "#CHROM_A")
input = fread("GUESS_BaseQTL_contactQTLs_hg19_rsids.txt")
lostrsid = input[! rsid%in%ld$ID_A]$rsid
lostsnpdt = data.table(CHROM_A=rep(NA_real_,length(lostrsid)), POS_A=rep(NA_real_,length(lostrsid)), ID_A=lostrsid, CHROM_B=rep(NA_real_,length(lostrsid)), POS_B=rep(NA_real_,length(lostrsid)), ID_B=lostrsid, UNPHASED_R2=1)
setnames(lostsnpdt, 1, "#CHROM_A")
ld1 = rbindlist(list(ld, samesnpdt, lostsnpdt))
nrow(ld1) # [1] 14877
fwrite(ld1, "GUESS_BaseQTL_contactQTLs_hg19_rsids_LD08.txt", sep="\t")

# Merge with gwas catalog
gwas = fread("gwas_catalog_v1.0-associations_e109_r2023-06-03.tsv") # downloaded from the GWAS catalog
ldgwas = merge(ld1, gwas, by.x="ID_B", by.y="SNPS")
length(unique(ldgwas$ID_B)) # [1] 299 unique SNPs
# ...in LD with
length(unique(ldgwas$ID_A))
# ...[1] 234 original triplet + contact eQTLs

fwrite(ldgwas, "gwascat_consensus_LD_r2unphased_08_revision_phase3.txt", sep="\t")