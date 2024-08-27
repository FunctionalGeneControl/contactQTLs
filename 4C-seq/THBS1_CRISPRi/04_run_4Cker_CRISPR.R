library(R.4Cker)
library(ggplot2)
library(DESeq2)
library(data.table)


### now trying using the files that were mapped using pipe4C.
### the bed file "flanking sites ..." can be replaced with the rmap, using the first three columns.
setwd("/data/cmn_vamal/4C/THBS1/4Cker")
#enz_file=read.table("./reducedGenome/hg38_dpnii_flanking_sites_47_unique_2.bed", stringsAsFactors = FALSE)
enz_file=read.table("./design/dpnII_AseI.rmap", stringsAsFactors = FALSE)

my_obj = createR4CkerObjectFromFiles(files = c("./mapped_with_pipe4C/counts/ScrR1_THBS1_rm_self_und.bedGraph",
                                               "./mapped_with_pipe4C/counts/ScrR2_THBS1_rm_self_und.bedGraph", 
					       "./mapped_with_pipe4C/counts/ScrR3_THBS1_rm_self_und.bedGraph", 
                                               "./mapped_with_pipe4C/counts/SNPsR1_THBS1_rm_self_und.bedGraph",
					       "./mapped_with_pipe4C/counts/SNPsR2_THBS1_rm_self_und.bedGraph", 
					       "./mapped_with_pipe4C/counts/SNPsR3_THBS1_rm_self_und.bedGraph"), 
				     bait_chr="chr15", 
				     bait_coord= 39315242, 
				     bait_name = "THBS1_QTLs", 
				     primary_enz = "GATC", 
				     samples = c("A", "B", "C", "A", "B", "C"), 
				     conditions = c("scrambled", "SNPs"), 
				     replicates = c(3,3),
				     species = "hg", 
				     output_dir = "./mapped_with_pipe4C/4Cker_out_SNPs", 
				     enz_file=enz_file)

## this was a test to see what happens when we do not remove proximal ligations. Not much changed.
#my_obj = createR4CkerObjectFromFiles(files = c("./mapped_with_pipe4C/counts/ScrR1_THBS1.bedGraph",
#                                              "./mapped_with_pipe4C/counts/ScrR2_THBS1.bedGraph", 
#					       "./mapped_with_pipe4C/counts/ScrR3_THBS1.bedGraph", 
#                                              "./mapped_with_pipe4C/counts/SNPsR1_THBS1.bedGraph",
#					       "./mapped_with_pipe4C/counts/SNPsR2_THBS1.bedGraph", 
#					       "./mapped_with_pipe4C/counts/SNPsR3_THBS1.bedGraph"), 
#				     bait_chr="chr15", 
#				     bait_coord= 39315242, 
#				     bait_name = "THBS1_QTLs", 
#				     primary_enz = "GATC", 
#				     samples = c("1_scrambled", "2_scrambled", "3_scrambled", "1_SNPs", "2_SNPs", "3_SNPs"), 
#				     conditions = c("scrambled", "SNPs"), 
#				     replicates = c(3,3),
#				     species = "hg", 
#				     output_dir = "./mapped_with_pipe4C/4Cker_out_SNPs", 
#				     enz_file=enz_file)                                    

#### below: do this using my paired analysis.

nb_results=nearBaitAnalysis(my_obj,k=30) # they recommend 5 but 30 works better here
pdf(file="./mapped_with_pipe4C/4Cker_out_SNPs/nearBaitAnalysis_K30.pdf", width = 7, height = 4)
ggplot(nb_results$norm_counts_avg, aes(x=Coord, y=Count, colour=Condition))+
	xlim(38800000,39800000)+
	theme_bw()+
	geom_vline(xintercept = 39581000, color = "gray40", linetype = "dashed")+
	geom_line(alpha = 0.5)+xlab(paste("Chromosome coordinates (", my_obj@bait_chr, ")", sep =""))+
	ylab("Normalized counts")+
	ggtitle(paste("Near bait analysis (", my_obj@bait_name, " bait)", sep = ""))
dev.off()


# now trying with rlog - This looks better
source("~/helen_Rfunctions.R")

pdf(file="./mapped_with_pipe4C/4Cker_out_SNPs/dif_highP_rlog_K30.pdf", width = 8, height = 4)
res_df = differentialAnalysis_rlog(obj=my_obj,
			norm_counts_avg=nb_results$norm_counts_avg,
			windows=nb_results$window_counts,
			conditions=c("scrambled", "SNPs"),
			region="nearbait",
			coordinates=NULL,
			pval=0.05, baitloc = 39315242, 
			oemin = 39581079, oemax = 39599466, # highlighting the THBS1 gene
                        plotStart = 39000000, plotEnd = 39700000, 
                        yMin = 0, yMax = 17) 
dev.off()
res_dt <- as.data.table(res_df)
fwrite(res_dt, file = "./mapped_with_pipe4C/4Cker_out_SNPs/rs2033937_CRISPR_dif_highP_results_pvals.txt", sep = "\t", quote = F, row.names = F, col.names = T)



#### plot the read counts at THBS1 and the R/L ratio.
scr1 <- as.data.table(my_obj@data_nearbait[1])
scr2 <- as.data.table(my_obj@data_nearbait[2])
scr3 <- as.data.table(my_obj@data_nearbait[3])
snps1 <- as.data.table(my_obj@data_nearbait[4])
snps2 <- as.data.table(my_obj@data_nearbait[5])
snps3 <- as.data.table(my_obj@data_nearbait[6])

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

scr1_THBS1 <- get_total_reads(scr1, start, end)
scr2_THBS1 <- get_total_reads(scr2, start, end)
scr3_THBS1 <- get_total_reads(scr3, start, end)
snps1_THBS1 <- get_total_reads(snps1, start, end)
snps2_THBS1 <- get_total_reads(snps2, start, end)
snps3_THBS1 <- get_total_reads(snps3, start, end)

### make a graph
toplot <- data.table(sample = c("Scr1", "Scr2", "Scr3", "SNPs1", "SNPs2", "SNPs3"),
                     Normalised_reads_across_THBS1 = c(scr1_THBS1, scr2_THBS1, scr3_THBS1, snps1_THBS1, snps2_THBS1, snps3_THBS1))

pdf(file = "./mapped_with_pipe4C/4Cker_out_SNPs/THBS1_CRISPRi_sample_comparison.pdf")
p <- ggplot(toplot, aes(x = sample, y = Normalised_reads_across_THBS1)) +
        geom_bar(position="dodge", stat="identity") +
       theme(text = element_text(size = 15))
print(p)
dev.off()

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
        mysample <- deparse(substitute(counts))
        reads <- data.table("left" = leftsum, "right" = rightsum,
                            "Sample" = mysample)
        return(reads)
}
scr1_LR <- get_left_right(scr1)
scr2_LR <- get_left_right(scr2)
scr3_LR <- get_left_right(scr3)
snps1_LR <- get_left_right(snps1)
snps2_LR <- get_left_right(snps2)
snps3_LR <- get_left_right(snps3)

LR <- rbind(scr1_LR, scr2_LR, scr3_LR, snps1_LR, snps2_LR, snps3_LR)
LR[, RL_ratio := right/left]
fwrite(LR, file = "./mapped_with_pipe4C/4Cker_out_SNPs/rightLeftRatio_numbers.txt", sep = "\t", quote = F, row.names = F, col.names = T)

### Make a graph
pdf(file = "./mapped_with_pipe4C/4Cker_out_SNPs/RightLeftRatio.pdf")
p <- ggplot(LR, aes(x = Sample, y = RL_ratio)) +
        geom_bar(position="dodge", stat="identity") +
       theme(text = element_text(size = 15))
print(p)
dev.off()

#### The direction consistent with the allele specific direction.
LR[, condition := c("Control", "Control", "Control", "CRISPRi", "CRISPRi", "CRISPRi")]
LR[, Replicate := c(1, 2, 3, 1, 2, 3)]
pdf(file = "./mapped_with_pipe4C/4Cker_out_SNPs/RightLeftRatio.pdf")
p <- ggplot(LR, aes(x = Replicate, y = RL_ratio, fill = condition)) +
        geom_bar(position="dodge", stat="identity") +
       theme(text = element_text(size = 15))
print(p)
dev.off()

LR[, condition := c("Control", "Control", "Control", "CRISPRi", "CRISPRi", "CRISPRi")]
LR[, Replicate := c(1, 2, 3, 1, 2, 3)]
pdf(file = "./mapped_with_pipe4C/4Cker_out_SNPs/RightLeftRatio_linegraph_CRISPR_SNPs.pdf", width = 4, height = 6)
p <- ggplot(LR, aes(x = condition, y = RL_ratio, group = Replicate)) +
        geom_line(color = "black") + geom_point(size = 2, color = "red") + theme_minimal() +
       theme(text = element_text(size = 15), axis.title.y = element_text(angle = 90, vjust = 0.5)) + 
       ylab("Right-left ratio") + xlab("Condition")
print(p)
dev.off()

### Plotting the whole region in CRISPRi was done using Plotgardener.






