library(R.4Cker)
library(ggplot2)
library(DESeq2)
library(data.table)
setwd("/data/cmn_vamal/4C/THBS1/4Cker")

#enz_file=read.table("./reducedGenome/hg38_dpnii_flanking_sites_47_unique_2.bed", stringsAsFactors = FALSE)
enz_file=read.table("./design/dpnII_AseI.rmap", stringsAsFactors = FALSE)

my_obj = createR4CkerObjectFromFiles(files = c("./mapped_with_pipe4C/counts/M4_T_rm_self_und.bedGraph", 
					       "./mapped_with_pipe4C/counts/M9_T_rm_self_und.bedGraph",
					      "./mapped_with_pipe4C/counts/M13_T_rm_self_und.bedGraph", 
					       "./mapped_with_pipe4C/counts/M4_C_rm_self_und.bedGraph", 
					       "./mapped_with_pipe4C/counts/M9_C_rm_self_und.bedGraph", 
					       "./mapped_with_pipe4C/counts/M13_C_rm_self_und.bedGraph"), 
				     bait_chr="chr15", 
				     bait_coord= 39315242, 
				     bait_name = "THBS1_QTLs", 
				     primary_enz = "GATC", 
				     samples = c("M4", "M9", "M13", "M4", "M9", "M13"), 
				     conditions = c("rs2033937_T", "rs2033937_C"), 
				     replicates = c(3,3),
				     species = "hg", 
				     output_dir = "./mapped_with_pipe4C/4Cker_out_rs2033937", 
				     enz_file=enz_file)

nb_results=nearBaitAnalysis(my_obj,k=30) # they reccommend 5 but increasing it to 30 lookes better. I think it's to do with fragments dropping out.

pdf(file="./mapped_with_pipe4C/4Cker_out_rs2033937/rs2033937_nearBaitAnalysis.pdf", width = 7, height = 4)
ggplot(nb_results$norm_counts_avg, aes(x=Coord, y=Count, colour=Condition))+
	xlim(38800000,39800000)+
	theme_bw()+
	geom_vline(xintercept = 39581000, color = "gray40", linetype = "dashed")+
	geom_line(alpha = 0.5)+xlab(paste("Chromosome coordinates (", my_obj@bait_chr, ")", sep =""))+
	ylab("Normalized counts")+
	ggtitle(paste("Near bait analysis (", my_obj@bait_name, " bait)", sep = ""))
dev.off()

### Running differential analysis using paired DESeq.
source("~/helen_Rfunctions.R")

pdf(file="./mapped_with_pipe4C/4Cker_out_rs2033937/rs2033937_dif_highP_paired.pdf", width = 8, height = 4)
res_df = differentialAnalysis_rlog(obj=my_obj,
			norm_counts_avg=nb_results$norm_counts_avg,
			windows=nb_results$window_counts,
			conditions=c("rs2033937_T", "rs2033937_C"),
			region="nearbait",
			coordinates=NULL,
			pval=0.05, baitloc = 39315242, 
			oemin = 39581079, oemax = 39599466, # highlighting the THBS1 gene
                        plotStart = 39000000, plotEnd = 39700000, 
                        yMin = 0, yMax = 17) 
dev.off()


res_dt <- as.data.table(res_df)
fwrite(res_dt, file = "./mapped_with_pipe4C/4Cker_out_rs2033937/rs2033937_dif_highP_results_pvals.txt", sep = "\t", quote = F, row.names = F, col.names = T)
# padj is sig from 39563594 to 39605081 (encompasses THBS1)


#########################################


### plot one sample versus another.
library(data.table)
M4T <- as.data.table(my_obj@data_nearbait[1])
M9T <- as.data.table(my_obj@data_nearbait[2])
M13T <- as.data.table(my_obj@data_nearbait[3])
M4C <- as.data.table(my_obj@data_nearbait[4])
M9C <- as.data.table(my_obj@data_nearbait[5])
M13C <- as.data.table(my_obj@data_nearbait[6])

# Region THBS1: 
#chr15:39581079-39599466
startgene <- 39581079
end <- 39599466
# get promoter as well. 300bp upstream of start
start <- startgene - 300

get_total_reads <- function(sample, start, end) {
	mygene <- sample[V3 < end & V3 > start | V2 > start & V2 < end]
	mygene[, V5 := sum(V4)] # summing across the whole of the gene
	mycount <- unique(mygene$V5)
	return(mycount)
}

M4T_THBS1 <- get_total_reads(M4T, start, end)
M9T_THBS1 <- get_total_reads(M9T, start, end)
M13T_THBS1 <- get_total_reads(M13T, start, end)
M4C_THBS1 <- get_total_reads(M4C, start, end)
M9C_THBS1 <- get_total_reads(M9C, start, end)
M13C_THBS1 <- get_total_reads(M13C, start, end)

### make a graph
toplot <- data.table(sample = c("M4", "M4", "M9", "M9", "M13", "M13"), 
		     Normalised_reads_across_THBS1 = c(M4T_THBS1, M4C_THBS1, M9T_THBS1, M9C_THBS1, M13T_THBS1, M13C_THBS1), 
		     Allele = c("T", "C", "T", "C", "T", "C"))

pdf(file = "./mapped_with_pipe4C/4Cker_out_rs2033937/rs2033937_alleles_sample_comparison.pdf")
p <- ggplot(toplot, aes(x = sample, y = Normalised_reads_across_THBS1, fill = Allele)) +
	geom_bar(position="dodge", stat="identity") +
       theme(text = element_text(size = 15))
print(p)
dev.off()

T <- c(toplot[Allele == "T", Normalised_reads_across_THBS1])
C <- c(toplot[Allele == "C", Normalised_reads_across_THBS1])
t.test(T, C, paired = TRUE, alternative = "two.sided")

### Left to right ratio for near bait
# my bait is chr15:39315234-39315580

baitstart <- 39315234
baitend <- 39315580
library(stringr)
get_left_right <- function(counts, baitleft, baitright) {
	left <- counts[V3 < baitstart]
	right <- counts[V2 > baitend]
	left[, V5 := sum(V4)]
	right[, V5 := sum(V4)]
	leftsum <- unique(left$V5)
        rightsum <- unique(right$V5)
	myname <- deparse(substitute(counts))
	myallele <- str_sub(myname,start=-1) 
	mysample <- str_replace_all(myname, myallele, "")
	reads <- data.table("left" = leftsum, "right" = rightsum, 
			    "Allele" = myallele, "Sample" = mysample)
	return(reads)
}

M4T_LR <- get_left_right(M4T)
M4C_LR <- get_left_right(M4C)
M9T_LR <- get_left_right(M9T)
M9C_LR <- get_left_right(M9C)
M13T_LR <- get_left_right(M13T)
M13C_LR <- get_left_right(M13C)

LR <- rbind(M4T_LR, M4C_LR, M9T_LR, M9C_LR, M13T_LR, M13C_LR)
LR[, RL_ratio := right/left]
LR[, sample := c("M4", "M4", "M9", "M9", "M13", "M13")]
LR[, allele := c("T", "C", "T", "C", "T", "C")]

fwrite(LR, file = "./mapped_with_pipe4C/4Cker_out_rs2033937/rightLeftRatio_numbers.txt", sep = "\t", quote = F, row.names = F, col.names = TRUE)

### Make a graph
pdf(file = "./mapped_with_pipe4C/4Cker_out_rs2033937/RightLeftRatio.pdf")
p <- ggplot(LR, aes(x = Sample, y = RL_ratio, fill = Allele)) +
        geom_bar(position="dodge", stat="identity") +
       theme(text = element_text(size = 15))
print(p)
dev.off()


pdf(file = "./mapped_with_pipe4C/4Cker_out_rs2033937/RightLeftRatio_linegraph.pdf", width = 4, height = 6)
p <- ggplot(LR, aes(x = allele, y = RL_ratio, group = sample)) + 
geom_line(color = "black") + geom_point(size = 2, color = "red") + theme_minimal() +
        ylab("Right-left ratio, 4C-seq reads") +        
		theme(text = element_text(size = 15), axis.title.y = element_text(angle = 90, vjust = 0.5))
print(p)
dev.off()

