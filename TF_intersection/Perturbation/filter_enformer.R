library(data.table)
setDTthreads(4)

filterEnformer = function(a_corr, idcol="V1"){
  
  a_corr_chip = copy(a_corr[,names(a_corr)%like%"CHIP:", with=F])
  a_corr_chip = as.matrix(a_corr_chip)
  
  print(ncol(a_corr)) # 4926
  print(ncol(a_corr_chip)) # 3678
  
  rownames(a_corr_chip) = a_corr[[idcol]]
  colnames(a_corr_chip)<-gsub("CHIP:","", colnames(a_corr_chip)) 
  colnames(a_corr_chip)<-gsub("\\:.+","", colnames(a_corr_chip)) 
  
  a_corr_chip = a_corr_chip[, colnames(a_corr_chip)!="."]
  print(ncol(a_corr_chip)) # 3558
  
  colnames(a_corr_chip) <- gsub("eGFP-","", colnames(a_corr_chip)) 
  colnames(a_corr_chip) <- gsub("3xFLAG-","", colnames(a_corr_chip)) 
  colnames(a_corr_chip) <- gsub("h","", colnames(a_corr_chip)) # to deal with hBMAL1 and hHIF1A
  colnames(a_corr_chip)[colnames(a_corr_chip)=="CEBPb"] = "CEBPB"
  
  # remove histone mods based on viewing unique(colnames(a_corr_chip))
  
  a_corr_chip = a_corr_chip[, -which(colnames(a_corr_chip)%in%c("H2AFZ", "H2AK5ac",
                                                                "H2AK9ac",    "H2BK120ac"   ,    "H2BK12ac"     , 
                                                                "H2BK15ac"   ,     "H2BK20ac"    ,   
                                                                "H2BK5ac"  ,       "H3F3A",          
                                                                "H3K14ac"   ,      "H3K18ac"    ,     "H3K23ac"     ,  
                                                                "H3K23me2"  ,      "H3K27Ac"   ,      "H3K27ac",        
                                                                "H3K27me3"  ,      "H3K36me3"    ,    "H3K4ac"   ,    
                                                                "H3K4me1"     ,    "H3K4me2"    ,     "H3K4me3"  ,      
                                                                "H3K56ac" ,        "H3K79me1"   ,     "H3K79me2"    ,  
                                                                "H3K9ac"    ,      "H3K9me1"    ,     "H3K9me2"  ,   "H3T11p",   
                                                                "H3K9me3"       ,  "H3T11ph"      ,   "H4K12ac"    ,   
                                                                "H4K20me1"    ,    "H4K5ac"     ,     "H4K8ac"   ,      
                                                                "H4K91ac"))]
  
  # remove co-factors: EZH3, EED and YY1, also HDACs, POL2/RNAP II, EP300 and BRDs
  
  a_corr_chip = a_corr_chip[, -which(colnames(a_corr_chip)%in%  c("EZH2", "EZH2phosphoT487","EZH2pospoT487", "EED", "POLR2A","POLR2ApospoS5","POLR2ApospoS2",     "POLR2AphosphoS2","POLR2AphosphoS5", "POLR2B","POLR2G", "POLR2H",    "HDAC1", "HDAC2", "HDAC3", "HDAC6", "HDAC8", "BRD4", "BRD9", "YY1", "YY2", "EP300", "RNAPII")) ]
  
  # remove weird stuff
  
  a_corr_chip = a_corr_chip[, -which(colnames(a_corr_chip)%in%c("abcam", "active"))]
  
  print(ncol(a_corr_chip)) # 1433
  sort(unique(colnames(a_corr_chip)))
  
  print(length(unique(colnames(a_corr_chip))))  # 662
  
  # take max across identically named columns 
  
  acmb = matrix(nrow=nrow(a_corr_chip), ncol=length(unique(colnames(a_corr_chip))))
  rownames(acmb)<-rownames(a_corr_chip) 
  colnames(acmb)<-unique(colnames(a_corr_chip))
  
  for(tf in colnames(acmb)){
    print(tf)
    if(length(colnames(a_corr_chip)[colnames(a_corr_chip)==tf])>1){
      # note in this case we can't just do a_corr_chip[, tf] - it'll pick just one (the first one, I think)
      acmb[, tf] <- apply(a_corr_chip[,colnames(a_corr_chip)==tf],1,max)
    }else{
      acmb[, tf] <- a_corr_chip[, tf]
    }
  }
  
  return(acmb)
  
}

message("Processing base_triplet")
a = fread("sadsar_base_triplet_revised.txt")
acmb = filterEnformer(a)
acmb = as.data.frame(acmb)
acmb[,"id"] = rownames(acmb)
setDT(acmb)
setcolorder(acmb, c(ncol(acmb),1:(ncol(acmb)-1)))
fwrite(acmb, "sadsar_base_triplet_filtmax_revised.txt", sep="\t")

message("Processing random [hold on tight]")
r = fread("sadsar_random.txt.gz")
rcmb = filterEnformer(r)
rcmb = as.data.frame(rcmb)
rcmb[,"id"] = rownames(rcmb)
setDT(rcmb)
setcolorder(rcmb, c(ncol(rcmb),1:(ncol(rcmb)-1)))
fwrite(rcmb, "sadsar_random_filtmax.txt", sep="\t")
