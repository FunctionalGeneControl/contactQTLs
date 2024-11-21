library(data.table)
setDTthreads(4)

processDeepSea = function(sample_prefix){
  
  ref = fread(paste0(sample_prefix,".ref_predictions.tsv"))
  print(head(ref[,1:9]))
  alt = fread(paste0(sample_prefix,".alt_predictions.tsv"))
  print(head(alt[,1:9]))
  
  refmat = as.matrix(ref[, -(1:9)])
  rownames(refmat) = ref$name
  print(head(refmat[,1:10]))
  
  altmat = as.matrix(alt[, -(1:9)])
  rownames(altmat) = alt$name
  print(head(altmat[,1:10]))
  
  # Remove sample name to just retain TF name
  tf = toupper(gsub("\\S+\\|(\\S+)\\|\\S+", "\\1", colnames(altmat)))
  tf[tf=="PU.1"] <- "SPI1"
  tf[tf=="C-MYC"] <- "MYC"
  tf[tf=="C-JUN"] <- "JUN"
  tf[tf=="C-FOS"] <- "FOS"
  tf[tf=="GATA-1"] <- "GATA1"
  tf[tf=="GATA-2"] <- "GATA2"   
  tf[tf=="EGFP-GATA2"] <- "GATA2"   
  tf[tf=="EGFP-FOS"] <- "FOS"
  tf[tf=="EGFP-JUNB"] <- "JUNB"
  tf[tf=="EGFP-JUND"] <- "JUND"
  tf[tf=="PAX5-C20"] <- "PAX5"   
  tf[tf=="PAX5-N19"] <- "PAX5"
  tf[tf=="NF-E2"] <- "NFE2"
  tf[tf=="NF-YA"] <- "NFYA"
  tf[tf=="NF-YB"] <- "NFYB"
  
  print(all(colnames(altmat)==colnames(refmat))) # TRUE
  colnames(refmat) = colnames(altmat) = tf
  
  print(dim(refmat)) #  1464  2002
  print(dim(altmat)) #  1464  2002
  
  # Filter out non-TFs
  
  tf.filt = !tf %in% c("H3K27ME3",            "H3K36ME3"           ,
                       "H3K4ME1",             "H3K4ME3",             "H3K9AC",              "H3K9ME3"        ,     "DNASE.ALL.PEAKS"    ,
                       "DNASE.FDR0.01.HOT",   "DNASE.FDR0.01.PEAKS" ,"DNASE.HOT"  ,         "H2AK5AC"    ,         "H2A.Z"   ,           
                       "H2BK120AC"    ,       "H2BK12AC"     ,       "H2BK15AC"      ,      "H2BK20AC"      ,      "H2BK5AC" ,           
                       "H3K14AC"      ,       "H3K18AC"       ,      "H3K23AC"      ,       "H3K23ME2"     ,       "H3K27AC"  ,          
                       "H3K4AC"           ,   "H3K4ME2"     ,        "H3K56AC"     ,        "H3K79ME1"        ,    "H3K79ME2"     ,      
                       "H4K20ME1"      ,      "H4K5AC"          ,    "H4K8AC"      ,        "H4K91AC"   ,          "H4K12AC"   ,         
                       "H3T11PH" ,            "H2AK9AC"     ,        "H3K9ME1"   , "DNASE",  "POL2"      ,          "POL2-4H8"    ,        
                       "POL2(B)"     ,        "POL2(PHOSPHOS2)"  ,   "POL3"     , "HDAC1"        ,       "HDAC2"   ,            "HDAC6", "EGFP-HDAC8",
                       "P300")        
  
  refmat.f = refmat[, tf.filt]
  altmat.f = altmat[, tf.filt]
  
  print(dim(refmat.f)) #  1464  595
  print(dim(altmat.f)) #  1464  595
  
  sadsar = (altmat.f-refmat.f)*(log(altmat.f+1,2)-log(refmat.f+1,2))
  print(dim(sadsar)) # 1464  595
  
  # combine per TF by max sadsar
  sadsarMax = matrix(nrow=nrow(sadsar), ncol=length(unique(colnames(sadsar))))
  rownames(sadsarMax) = rownames(sadsar)
  colnames(sadsarMax) = unique(colnames(sadsar))
  
  for(tf in colnames(sadsarMax)){
    print(tf)
    if(length(colnames(sadsar)[colnames(sadsar)==tf])>1){
      sadsarMax[, tf] <- apply(sadsar[,colnames(sadsar)==tf],1,max)
    }else{
      sadsarMax[, tf] <- sadsar[, tf]
    }
  }
  
  print(dim(sadsarMax)) # 1464 157
  
  sadsarMaxDT = as.data.table(sadsarMax)
  sadsarMaxDT[, rsid:=rownames(sadsarMax)]
  setcolorder(sadsarMaxDT, c("rsid",colnames(sadsarMax)))
  print(head(sadsarMaxDT[, 1:10]))
  
  fwrite(sadsarMaxDT, paste0(sample_prefix, ".sadsarMax.txt"), sep="\t")
}

# Files obtained from the running the online DeepSea platform on a random 1% 1000 Genomes SNPs
# (Chunked by <20000 SNPs per file)
processDeepSea("f6afee88-da11-4e07-85a3-024c796f949b_random_1KG_deepSea_input_hg38_chr1")
processDeepSea("b407d4ed-54a4-4fea-b832-4d12f0752b61_random_1KG_deepSea_input_hg38_chr2_chunk_2")
processDeepSea("dcb24f42-e08e-4e06-876c-ed35aa8def72_random_1KG_deepSea_input_hg38_chr10")
processDeepSea("d76931bd-2d96-42e6-990c-0282b4c3bd42_random_1KG_deepSea_input_hg38_chr6")
processDeepSea("b51f1a1a-dea4-45d1-ac54-e277e65ddefa_random_1KG_deepSea_input_hg38_chr9")
processDeepSea("a6321439-babf-41b3-a183-00d49f69ca06_random_1KG_deepSea_input_hg38_chr21")
processDeepSea("99f52cff-faff-4f97-be5e-959296e10593_random_1KG_deepSea_input_hg38_chr2_chunk_1")
processDeepSea("95eab4bc-e323-4fad-a5c8-1767cfb785ad_random_1KG_deepSea_input_hg38_chr7")
processDeepSea("8c44ac89-9e6a-4778-891a-752443692051_random_1KG_deepSea_input_hg38_chr5")
processDeepSea("8200d6e1-606e-44f6-ac66-719fa8b706ff_random_1KG_deepSea_input_hg38_chr12")
processDeepSea("73e63a28-a079-434a-a4fc-cc6915b5a3e2_random_1KG_deepSea_input_hg38_chr17")
processDeepSea("7246c4b5-14b8-4cb0-9a26-1f143798e120_random_1KG_deepSea_input_hg38_chr19")
processDeepSea("6f640d91-ebc8-44d9-a0dd-cd7128a7e654_random_1KG_deepSea_input_hg38_chr8")
processDeepSea("698edb8d-ac01-4bf8-b979-157417f2056b_random_1KG_deepSea_input_hg38_chr14")
processDeepSea("5b68f740-b40f-41dd-8197-5795a5aaeca0_random_1KG_deepSea_input_hg38_chr16")
processDeepSea("57f138d2-eccb-46ab-a99c-43d9a1e527ae_random_1KG_deepSea_input_hg38_chr4")
processDeepSea("54fb1f84-6f0b-43bc-b8fc-8118f442712e_random_1KG_deepSea_input_hg38_chr11")
processDeepSea("531fa29b-b3c5-4c5e-b8a4-517d90e94f3e_random_1KG_deepSea_input_hg38_chr13")
processDeepSea("44e1ede7-d5e8-4093-9175-65f315eb70d0_random_1KG_deepSea_input_hg38_chr15")
processDeepSea("3c445db8-5c1c-47d2-b17c-0a9599913d3c_random_1KG_deepSea_input_hg38_chr20")
processDeepSea("31beca07-faa5-4b22-ac5d-658ff68fc9dd_random_1KG_deepSea_input_hg38_chr22")
processDeepSea("0f615a80-b42d-41c3-9864-c44d6f6f60ae_random_1KG_deepSea_input_hg38_chr18")
processDeepSea("0fc215cd-b90b-4461-a824-5138cbac0723_random_1KG_deepSea_input_hg38_chr3")

ssMaxFiles = list.files(pattern="sadsarMax")
ssMaxList = vector("list")
for (f in ssMaxFiles){
  print(f)
  ssMaxList[[f]] = fread(f)
}
ssMax= rbindlist(ssMaxList)
dim(ssMax) # 251656    157
rm(ssMaxList)

thresholds = apply(ssMax[,-1], 2, function(x)quantile(x,0.99))
fwrite(data.table(tf=names(top1p), top1p_sadsar=top1p), "deepSea_sadsar_top1p_thresholds_2.txt", sep="\t")

### Now process the cQTLs result obtained from the online platform

processDeepSea("fd47b487-0571-468f-9d94-9db21c3cd030_GUESS_BaseQTL_contactQTLs_hg19_rsids_DeepSEA_input")
sadsar = fread("DeepSEA_16Apr/fd47b487-0571-468f-9d94-9db21c3cd030_GUESS_BaseQTL_contactQTLs_hg19_rsids_DeepSEA_input.sadsarMax.txt")
thresholds = fread("deepSea_sadsar_top1p_thresholds_2.txt")

sadsar.binary = copy(sadsar)
for(thistf in names(sadsar)[-1]){
  print(thistf)
  sadsar.binary[[thistf]] = sadsar[[thistf]]>thresholds[tf==thistf]$top1p_sadsar
  print(sum(sadsar.binary[[thistf]]))
}
fwrite(sadsar.binary, "deepSea_q01_base_triplet.txt", sep="\t")