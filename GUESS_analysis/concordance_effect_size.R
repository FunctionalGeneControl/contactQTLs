library(data.table)
setwd("~/OneDrive - Imperial College London/eQTL_capture_HiC/")
# conc = fread("GUESS_conc_disc.txt")
# GUESS_chic_beta_std = scale(conc$GUESS_chic_beta)
# GUESS_atac_beta_std = scale(conc$GUESS_ATAC_beta)
# GUESS_RNA_beta_std = scale(conc$GUESS_RNA_beta)
# means = rowMeans(cbind(abs(GUESS_chic_beta_std), abs(GUESS_atac_beta_std), abs(GUESS_RNA_beta_std)))
# boxplot(list(Concordant=means[conc$overall_concordance==TRUE], Discordant=means[conc$overall_concordance==FALSE]))


guess = fread("Table_S9_GUESS_Rev.txt") # 919
guess1 = guess[-grep(",", guess[["FDR_rsids"]])] # 869
effects = strsplit(guess1[["Post-proc. best model post. E(Beta|Y)"]], ",")
# now, where we have >3 elements, it's because we have data for SNPs that didn't pass FDR
# find which one is the FDR snp in them and retain only those
for(i in 1:length(effects)){
  if(length(effects[[i]])>3){
    allSnps = strsplit(guess1[["Post-proc. best model SNP"]][i], ",")[[1]]
    k=which(allSnps==guess1[["Post-proc. FDR mPPI SNP"]][i])
    start = 3*(k-1)+1
    effects[[i]] = c(effects[[i]][start],effects[[i]][start+1], effects[[i]][start+2])
  }
}
effectsmat = t(sapply(effects, as.numeric))
conc = apply(effectsmat, 1, function(x)sign(x[1])==sign(x[2])&sign(x[2])==sign(x[3]))
table(conc)
# FALSE  TRUE 
# 469   400 
# plotting standardised effect sizes - for this we also need SE

se = strsplit(guess1[["Post-proc. best model post. SE(Beta|Y)"]], ",")
# now, where we have >3 elements, it's because we have data for SNPs that didn't pass FDR
# find which one is the FDR snp in them and retain only those
for(i in 1:length(se)){
  if(length(se[[i]])>3){
    allSnps = strsplit(guess1[["Post-proc. best model SNP"]][i], ",")[[1]]
    k=which(allSnps==guess1[["Post-proc. FDR mPPI SNP"]][i])
    start = 3*(k-1)+1
    se[[i]] = c(se[[i]][start],se[[i]][start+1], se[[i]][start+2])
  }
}
semat = t(sapply(se, as.numeric))

# Supplementary Figure 5A
boxplot(list(Concordant=rowMeans(abs(effectsmat[conc==T,]/semat[conc==T,])), 
             Discordant=rowMeans(abs(effectsmat[conc==F,]/semat[conc==F,]))), 
        ylab="mean abs standardised effect size across 3 modalities", xlab="Effect concordance across 3 modalities")

cat(rowMeans(abs(effectsmat[conc==T,]/semat[conc==T,])), sep="\n", file="supl_fig5a_concordant.txt")
cat(rowMeans(abs(effectsmat[conc==F,]/semat[conc==F,])), sep="\n", file="supl_fig5a_discordant.txt")
               