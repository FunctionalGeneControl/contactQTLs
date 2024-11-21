library(data.table)
deepsea = fread("deepSea_q01_base_triplet.txt")
enf = fread("q01_base_triplet_revised.txt")
tfbw = fread("cQTL_intersect_with_Remap_Tob_Max_ChIP_wide.txt")
dim(tfbw) # 205 297 - so evidence of 296 TFBSs in monocytes (+rsid) at 205/641 cQTLs

enf = fread("q01_base_triplet_revised.txt")
dim(enf) # 641 663 - so 662 TFs + the rsid column
ds = fread("deepSea_q01_base_triplet.txt") # note the TF names are already changed to be consistent with Enformer
dim(ds) # 641 158 - so 157 TFs + the rsid column

# make sure the column order is the same for all of three of them and convert them all to matrices
# also just retain the overlapping set of TFs
TFBSnames = colnames(tfbw)[-1]
enfnames = colnames(enf)[-1]
dsnames = colnames(ds)[-1]

pertnames = enfnames[enfnames %in% dsnames]
length(pertnames) # 113 - so for 113 TFs we have evidence of perturbation from both approaches
tfenfnames = enfnames[enfnames %in% TFBSnames]
length(tfenfnames) # 179 - so for 179 TFs we have evidence for binding and Enformer
tfdsnames = dsnames[dsnames %in% TFBSnames]
length(tfdsnames) # 77 - so for 77 TFs we have evidence for binding and DeepSea

consnames = TFBSnames[TFBSnames%in% enfnames & TFBSnames%in% dsnames]
length(consnames) # 72 - so for 72 TFs we have evidence for both binding and perturbation from both approaches

unionnames = TFBSnames[TFBSnames %in% c(enfnames, dsnames)]
length(unionnames) # 184 - so for 184/296 we have perturbation evidence from at least one approach

tfbm = as.matrix(tfbw[,-1])
rownames(tfbm) = tfbw$rsid
enfm = as.matrix(enf[,-1])
rownames(enfm) = enf$id
dsm = as.matrix(ds[,-1])
rownames(dsm) = ds$rsid

boundEnf = matrix(nrow=nrow(tfbm), ncol=length(tfenfnames))
rownames(boundEnf) = rownames(tfbm) # important!
colnames(boundEnf) = sort(tfenfnames)
enfm_b = enfm[rownames(boundEnf),] # important!
for(tf in colnames(boundEnf)){
  boundEnf[,tf] = (enfm_b[,tf]==TRUE) & (tfbm[,tf] != "")
}
sort(colSums(boundEnf), decreasing = T)
# SPI1   CEBPB   STAT3   IKZF1    NFIC    JUND   NR2F1   NR2F2     FOS   TCF12 BHLHE40  STAT5A    CREM    JUNB     MAX     TBP 
# 33      17      15      14      13      11      11      11      10      10       9       9       8       7       7       7 
nrow(boundEnf[rowSums(boundEnf)>0,]) # 62 - i.e., 9.6% total; 30% bound
write.table(boundEnf, "bound_enformer_cQTLs.txt", quote=F, sep="\t")

boundDs = matrix(nrow=nrow(tfbm), ncol=length(tfdsnames))
rownames(boundDs) = rownames(tfbm) # important!
colnames(boundDs) = sort(tfdsnames)
dsm_b = dsm[rownames(boundDs),] # important!
for(tf in colnames(boundDs)){
  boundDs[,tf] = (dsm_b[,tf]==TRUE) & (tfbm[,tf] != "")
}
sort(colSums(boundDs), decreasing = T)
# SPI1   CEBPB   IKZF1   STAT3   MEF2A    NFIC     MAX     SRF   TCF12    ELF1 BHLHE40   CREB1     FOS    MXI1   FOXM1    JUND 
# 30      23      23      23      20      17      16      16      16      14      13      13      13      13      12      12 
nrow(boundDs[rowSums(boundDs)>0,]) # 79 - i.e., 12.3% total; 38.5% bound
write.table(boundDs, "bound_deepSea_cQTLs.txt", quote=F, sep="\t")

# Used for downstream analysis
boundCons = matrix(nrow=nrow(tfbm), ncol=length(consnames))
rownames(boundCons) = rownames(tfbm) # important!
colnames(boundCons) = sort(consnames)
dsm_cb = dsm[rownames(boundCons),] # important!
enfm_cb = enfm[rownames(boundCons),] # important!
for(tf in colnames(boundCons)){
  boundCons[,tf] = (dsm_cb[,tf]==TRUE) & (enfm_cb[,tf]==TRUE) & (tfbm[,tf] != "")
}
sort(colSums(boundCons), decreasing = T)
# SPI1   CEBPB   STAT3   IKZF1     FOS    NFIC   TCF12     TBP BHLHE40    JUND     MAX    MXI1    PBX3  STAT5A   MEF2A     MYC 
# 19      13      12       9       8       7       7       6       5       5       5       5       5       5       4       4 
nrow(boundCons[rowSums(boundCons)>0,]) # 49 - i.e., 7.6% total; 23.9% bound
write.table(boundCons, "bound_deepSea_Enformer_cons_cQTLs.txt", quote=F, sep="\t")

# Used for the supplementary figure 
boundUnion = matrix(nrow=nrow(tfbm), ncol=length(unionnames))
rownames(boundUnion) = rownames(tfbm) # important!
colnames(boundUnion) = sort(unionnames)
dsm_ub = dsm[rownames(boundUnion),] # important!
enfm_ub = enfm[rownames(boundUnion),] # important!
for(tf in colnames(boundUnion)){
  if(tf %in% colnames(dsm_ub) & tf %in% colnames(enfm_ub)){
    boundUnion[,tf] = ((dsm_ub[,tf]==TRUE) | (enfm_ub[,tf]==TRUE)) & (tfbm[,tf] != "")
  }else{
    if(tf %in% colnames(dsm_ub)){
      boundUnion[,tf] = (dsm_ub[,tf]==TRUE) & (tfbm[,tf] != "")
    }else{ # then i should be in colnames(enfm_ub) by design
      boundUnion[,tf] = (enfm_ub[,tf]==TRUE) & (tfbm[,tf] != "")
    }
  }
}
sort(colSums(boundUnion), decreasing = T)
# SPI1   IKZF1   CEBPB   STAT3    NFIC   MEF2A   TCF12    JUND     MAX BHLHE40     SRF    ELF1     FOS   CREB1   FOXM1     JUN 
# 44      28      27      26      23      21      19      18      18      17      16      15      15      14      14      14 
nrow(boundUnion[rowSums(boundUnion)>0,]) # 87 - i.e., 13.5% total; 42.4% bound
write.table(boundUnion, "bound_deepSea_Enformer_union_cQTLs.txt", quote=F, sep="\t")

### Supplementary Figure 7B

jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}

DeepSea_Enformer_Jac = function(tf){
  posEnf = rownames(boundEnf[boundEnf[,tf]==TRUE,])
  posDeepSea = rownames(boundDs[boundDs[,tf]==TRUE,])
  Jac = jaccard(posEnf, posDeepSea)
  
  c(Jac, sum(boundUnion[, tf]))
}

DeepSea_Enformer_JacVec = Vectorize(DeepSea_Enformer_Jac, vectorize.args = "tf")
DeepSea_Enformer_JacVec(c("SPI1", "CEBPB", "STAT3"))

res = DeepSea_Enformer_JacVec(colnames(boundCons[,colSums(boundCons)>0]))
plot(res[2,], res[1,], pch=19, cex=1, col="maroon", ylab="Jaccard index (Enformer/DeepSea)", xlab="Union of predicted TF perturbations (Enformer+DeepSea)", 
     main="Agreement between DeepSea and Enformer")
text(res[2,], res[1,], labels=colnames(res), pos=4, cex=0.5)

### Figure 6B, left
### The final pie chart is assembled manually from the components plotted here
cQTLsummary = matrix(nrow=nrow(enf), ncol=5) # note that nrow(enf)=nrow(deepsea)=n(cQTLs)=641
rownames(cQTLsummary) = rownames(enfm)
colnames(cQTLsummary) = c("TFbound", "Enformer", "DeepSea", "Union", "Consensus")
cQTLsummary[rownames(tfbm), "TFbound"] = "yes"
cQTLsummary[rownames(boundEnf[rowSums(boundEnf)>0,]), "Enformer"] = "yes"
cQTLsummary[rownames(boundDs[rowSums(boundDs)>0,]), "DeepSea"] = "yes"
cQTLsummary[rownames(boundUnion[rowSums(boundUnion)>0,]), "Union"] = "yes"
cQTLsummary[rownames(boundUnion[rowSums(boundCons)>0,]), "Consensus"] = "yes"
cQTLsummary[is.na(cQTLsummary[,1]),1] = "no"
cQTLsummary[is.na(cQTLsummary[,2]),2] = "no"
cQTLsummary[is.na(cQTLsummary[,3]),3] = "no"
cQTLsummary[is.na(cQTLsummary[,4]),4] = "no"
cQTLsummary[is.na(cQTLsummary[,5]),5] = "no"
cQTLsummary = as.data.frame(cQTLsummary)
cQTLsummary$TFbound = factor(cQTLsummary$TFbound, levels=c("no", "yes"))
cQTLsummary$DeepSea = factor(cQTLsummary$DeepSea, levels=c("no", "yes"))
cQTLsummary$Enformer = factor(cQTLsummary$Enformer, levels=c("no", "yes"))
cQTLsummary$Either = factor(cQTLsummary$Union, levels=c("no", "yes"))
cQTLsummary$Consensus = factor(cQTLsummary$Consensus, levels=c("no", "yes"))

library(webr)
library(ggplot2)
PieDonut(cQTLsummary, aes(TFbound, Union), explode = T, start=1.5*pi/2)
PieDonut(cQTLsummary, aes(TFbound, Consensus), explode = T, start=1.5*pi/2)

#### Figure 6B, right

nc =  table(rowSums(boundCons))
ncbin = cut(as.numeric(names(nc)), breaks=c(-0.1,0,1,2,3,5,6,8,10,15,20,50))
nc_per_bin = tapply(nc,ncbin,sum)
barplot(nc_per_bin[-1], col="maroon", cex.names=0.8, las=2, 
        names.arg = c("1","2","3","4-5","6", "7-8","9-10","11-15","16-20", ">20"), 
        xlab="# perturbed TFs", ylab="# cQTLs", beside=T)

#### Figure 6C

mt1 = data.frame(tf=colnames(boundCons)[colSums(boundCons)>=5], 
                       n=colSums(boundCons)[colSums(boundCons)>=5])
mt1 = mt1[order(mt1$n),]
barplot(mt1$n, names.arg = mt1$tf, 
        cex.names=0.7, col = "#5F9C5F", las =2, ylab="# Perturbed cQTLs", horiz = T)

#### Overlap with published SPI1/PU.1 and CTCF tfQTLs

# Supplementary data from Watt et al., Nat Comms 2021
pu.1 = readxl::read_xlsx("41467_2021_22548_MOESM5_ESM.xlsx", sheet = 1)
pu1qtl = pu.1$rsid[pu.1$lFDR<0.05] # 5465

# Accounting for LD (R2>0.99). 
# Using such a narrow LD window because both tfQTLs and SPI1-perturbing cQTLs 
# were called in the vicinity of PU.1 ChIP peaks
# So increasing the LD window will just reduce the specificity of this analysis
ld1 = fread("GUESS_BaseQTL_contactQTLs_hg19_rsids_LD08.txt")
ld1 = ld1[UNPHASED_R2>0.99]
spi1ld = ld1[ID_A %in% rownames(boundCons[boundCons[,"SPI1"]==TRUE,])] 
predSPI1 = length(unique(spi1ld$ID_A)) # 19
predSPI1_tfQTL = length(unique(spi1ld[ID_B %in% pu1qtl]$ID_A)) # 5

non_pu1_cQTLs = ld1[!ID_A %in% rownames(boundCons[boundCons[,"SPI1"]==TRUE,])]
nonpredSPI1 = length(unique(non_pu1_cQTLs$ID_A)) # 622
nonpredSPI1_tfQTL = length(unique(non_pu1_cQTLs[ID_B %in% pu1qtl]$ID_A)) # 29

fisher.test(matrix(c(predSPI1_tfQTL, predSPI1-predSPI1_tfQTL,
                     nonpredSPI1_tfQTL, nonpredSPI1-nonpredSPI1_tfQTL),
                   nrow=2, byrow=T))$p.value # 0.002

# Supplementary data from Ding et al., PLoS Genet 2014
ctcfqtl = fread("pca1_01_cls.tsv")$VARIANT_ID
ctcf_cQTLs = rownames(boundCons[boundCons[,"CTCF"]==TRUE,])
ctcf_cQTLs[ctcf_cQTLs%in%ctcfqtl] # [1] "rs2353678" "rs7146599"

#### Number of cQTLs predicted to perturb the binding of either one of the three top TFs
nrow(boundCons[boundCons[,"SPI1"]==T | boundCons[,"CEBPB"]==T | boundCons[,"STAT3"]==T,]) # 31
