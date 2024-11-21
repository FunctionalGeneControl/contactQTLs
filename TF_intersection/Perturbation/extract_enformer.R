library(data.table)
library(rhdf5)

setDTthreads(4)

baseQTLs=unlist(fread("df_sig_scores_baseqtl.csv")[,1])
triplets = unlist(fread("p_scores_GUESS_28Feb24.csv")[,1])
base_triplets = unique(c(baseQTLs, triplets))
length(unique(c(baseQTLs, triplets))) # 641
base_triplets = unlist(fread("GUESS_BaseQTL_contactQTLs_hg19_rsids.txt")$rsid)

## found that 12 SNPs on chr12 go under different IDs
# pos = h5read("1000G.MAF_threshold=0.005.12.h5", name="pos")
# snps = h5read("1000G.MAF_threshold=0.005.12.h5", name="snp")
# look up the missing SNPs positions in dbSNP, and then look up SNP id in the enformer dataset by position, such as:
# obtain the recode ID dictionary
recodeIDs = c(rs7304465=	"ss1388058849",
              rs6488819= "ss1388058902",
              rs7958511=	"ss1388058951",
              rs10732571= "ss1388058995",
              rs7309151= "ss1388059064",
              rs7303600=	"ss1388059094",
              rs7299161=	"ss1388059156",
              rs7299583=	"ss1388059165",
              rs7303774=	"ss1388059207",
              rs7487016=	"ss1388059298",
              rs4883466= "ss1388059317",
              rs4883468=	"ss1388059324")

base_triplets[base_triplets %in% names(recodeIDs)] = recodeIDs

length(base_triplets) # 641

getSarSad = function(infile, snp_index, outfile, snps_list, datasets_list, recode=NULL){
  
  message("Extracting sar")
  sar = h5read(infile, name="SAR", index=list(NULL, snp_index))
  sar = t(sar)
  colnames(sar) <- datasets_list
  rownames(sar) <- snps_list[snp_index]
  message("Extracting sad")
  sad = h5read(infile, name="SAD", index=list(NULL, snp_index))
  sad = t(sad)
  colnames(sad) = datasets_list
  rownames(sad) = snps_list[snp_index]
  sadsar = sad*sar
  
  if(!is.null(recode)){
    if(any(rownames(sadsar)%in%recode)){
      message("reverting recoded IDs...")
      rownames(sadsar)[rownames(sadsar) %in% recode] <- names(recode[recode %in% rownames(sadsar)])
    }
  }
  write.table(sadsar, outfile, sep="\t", row.names=T, col.names=T, quote=F)
  
}

for(chr in 1:22){

  h5name = paste0("1000G.MAF_threshold=0.005.",chr,".h5")
  message("Extracting from", h5name)
  
  snps = h5read(h5name, name="snp")
  datasets = h5read(h5name, name="target_labels")
  
  bt_index = which(snps %in% base_triplets)
  message("Extracting for base+triplets")
  getSarSad(h5name, bt_index, paste0("sadsar_chr",chr,"_base_triplet_revised.txt"), snps, datasets, recodeIDs)
  
  set.seed(123)
  rand_index = sample(1:length(snps), length(snps)*0.01)
  message("Extracting random 1%")
  getSarSad(h5name, rand_index, paste0("sadsar_chr",chr,"_random.txt"), snps, datasets)
  
}

a = vector("list"); 
for (chr in 1:22){ a[[chr]] = fread(paste0("sadsar_chr", chr,"_base_triplet_revised.txt")) }; 
a = rbindlist(a)
a = a[!duplicated(V1)]
nrow(a) #  641 - as it should be
fwrite(a, "sadsar_base_triplet_revised.txt", sep="\t")

r = vector("list");
for (chr in 1:22){ r[[chr]] = fread(paste0("sadsar_chr", chr,"_random.txt")) };
r = rbindlist(r)
fwrite(r, "sasdsar_random.txt", sep="\t")

