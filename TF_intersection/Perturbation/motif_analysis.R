library(data.table); setDTthreads(4)
library(motifbreakR)
library(SNPlocs.Hsapiens.dbSNP155.GRCh37) # dbSNP155 in hg19
library(BSgenome.Hsapiens.UCSC.hg19)     # hg19 genome

### motifbreakR run and pruning

cQTLs = fread("GUESS_BaseQTL_contactQTLs_hg19_rsids.txt")

snps.all = snps.from.rsid(rsid = cQTLs$rsid,
                          dbSNP = SNPlocs.Hsapiens.dbSNP155.GRCh37,
                          search.genome = BSgenome.Hsapiens.UCSC.hg19)
results.all<- motifbreakR(snpList = snps.all, filterp = TRUE,
                          pwmList = subset(MotifDb, 
                                  dataSource %in% c("HOCOMOCOv11-core-A", "HOCOMOCOv11-core-B", "HOCOMOCOv11-core-C")),
                          threshold = 1e-4,
                          method = "ic",
                          bkg = c(A=0.25, C=0.25, G=0.25, T=0.25),
                          BPPARAM = BiocParallel::SerialParam())

# Make sure we only focus on the allelic change tested 
results.all.dt = as.data.frame(results.all, row.names = NULL)
setDT(results.all.dt)
res = merge(results.all.dt, cQTLs, by.x="SNP_id", by.y="rsid", allow.cartesian = TRUE)
res1 = res[ALT==gsub("\\S+\\:\\S+\\:\\S+\\:(\\S+)", "\\1", get("hg19QTL_ID"))]

nrow(results.all.dt) # [1] 2097
nrow(res1) # [1] 1171

saveRDS(res1[, 1:28], "results_allTriplets_dt_pruned_alleles_revised.Rds") 

### Downstream analysis for Figure 6

tfbw = fread("cQTL_intersect_with_Remap_Tob_Max_ChIP_wide.txt")
tfbm = as.matrix(tfbw[,-1])
rownames(tfbm) = tfbw$rsid

resMotif = readRDS("results_allTriplets_dt_pruned_alleles_revised.Rds")

boundMotifs = resMotif[apply(resMotif[, c("SNP_id", "geneSymbol")],1, function(x){
  if(x[2]%in%colnames(tfbm) & x[1]%in%rownames(tfbm)) { tfbm[x[1],x[2]]!="" }
  else{ FALSE} 
} )]
nrow(resMotif)
nrow(boundMotifs)

# How many cQTLs are predicted to perturb SPI1, STAT3 or CEBPB?
spi1.total = nrow(boundCons[boundCons[,"SPI1"]==TRUE,])
stat3.total = nrow(boundCons[boundCons[,"STAT3"]==TRUE,])
cebpb.total = nrow(boundCons[boundCons[,"CEBPB"]==TRUE,])

# How many cQTLs perturb the cognate motif of these factors?
spi1.cognate = nrow(boundMotifs[SNP_id %in% rownames(boundCons[boundCons[,"SPI1"]==TRUE,]) & geneSymbol=="SPI1"])
stat3.cognate = nrow(boundMotifs[SNP_id %in% rownames(boundCons[boundCons[,"STAT3"]==TRUE,]) & geneSymbol=="STAT3"])
cebpb.cognate = nrow(boundMotifs[SNP_id %in% rownames(boundCons[boundCons[,"CEBPB"]==TRUE,]) & geneSymbol=="CEBPB"])

# How many cQTLs that perturb these factors only perturb non-cognate motifs?
# Note we can't do nrow like above because there may be multiple motifs per cQTL
spi1.onlync = length(unique(boundMotifs[SNP_id %in% rownames(boundCons[boundCons[,"SPI1"]==TRUE,]) & !SNP_id %in% boundMotifs[geneSymbol=="SPI1"]$SNP_id ]$SNP_id))
stat3.onlync = length(unique(boundMotifs[SNP_id %in% rownames(boundCons[boundCons[,"STAT3"]==TRUE,]) & !SNP_id %in% boundMotifs[geneSymbol=="STAT3"]$SNP_id ]$SNP_id))
cebpb.onlync = length(unique(boundMotifs[SNP_id %in% rownames(boundCons[boundCons[,"CEBPB"]==TRUE,]) & !SNP_id %in% boundMotifs[geneSymbol=="CEBPB"]$SNP_id ]$SNP_id))

### Figure 6D
motifmat = matrix(c(stat3.cognate, stat3.total-stat3.cognate-stat3.onlync, stat3.onlync,
                    cebpb.cognate, cebpb.total-cebpb.cognate-cebpb.onlync, cebpb.onlync,
                    spi1.cognate, spi1.total-spi1.cognate-spi1.onlync, spi1.onlync), ncol=3, nrow=3)

barplot(motifmat, horiz = T, beside = F, col=c("#E3521A", "white", "darkgrey"), names.arg=c("STAT3", "CEBPB", "SPI1"), cex.names=0.9)

### Figure 6E
stat3_nc = sort(table(boundMotifs[SNP_id %in% rownames(boundCons[boundCons[,"STAT3"]==TRUE,]) & !SNP_id %in% boundMotifs[geneSymbol=="STAT3"]$SNP_id ]$geneSymbol))
barplot(stat3_nc, col="darkgrey", cex.names=0.9, horiz=T)
