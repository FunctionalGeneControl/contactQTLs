library(R.4Cker)
library(ggplot2)
library(DESeq2)
library(data.table)
source("~/helen_Rfunctions.R")
setwd("/data/cmn_vamal/4C/PCK2/4Cker")
enz_file=read.table("./design/dpnII_hpych4v.rmap", stringsAsFactors = FALSE)

my_obj = createR4CkerObjectFromFiles(files = c("./mapped_with_pipe4C/counts/M4_PCK2_G_rm_self_und.bedGraph",
                                               "./mapped_with_pipe4C/counts/M9_PCK2_G_rm_self_und.bedGraph",
                                               "./mapped_with_pipe4C/counts/M13_PCK2_G_rm_self_und.bedGraph",
                                               "./mapped_with_pipe4C/counts/M4_PCK2_A_rm_self_und.bedGraph",
                                               "./mapped_with_pipe4C/counts/M9_PCK2_A_rm_self_und.bedGraph",
                                               "./mapped_with_pipe4C/counts/M13_PCK2_A_rm_self_und.bedGraph"),
			             bait_chr="chr14",
                                     bait_coord= 24058527,
                                     bait_name = "PCK2_QTLs",
                                     primary_enz = "GATC",
                                     samples = c("M4", "M9", "M13", "M4", "M9", "M13"),
                                     conditions = c("rs7146599_G", "rs7146599_A"),
                                     replicates = c(3,3),
                                     species = "hg",
                                     output_dir = "./mapped_with_pipe4C/4Cker_out_rs7146599",
                                     enz_file=enz_file)

# G is ref, A is alt. G should interact more strongly with PCK2.
# PCK2 gene (with promoter on L): chr14:24,094,311-24,104,125

nb_results=nearBaitAnalysis(my_obj,k=30) # they reccommend 5 but increasing it to 30 lookes better - I think it's to do with fragments dropping out.

pdf(file="./mapped_with_pipe4C/4Cker_out_rs7146599/nearBaitAnalysis.pdf", width = 7, height = 4)
ggplot(nb_results$norm_counts_avg, aes(x=Coord, y=Count, colour=Condition))+
        xlim(24000000,24120000)+
	theme_bw()+
	geom_vline(xintercept = 24094400, color = "gray40", linetype = "dashed")+
        geom_line(alpha = 0.5)+xlab(paste("Chromosome coordinates (", my_obj@bait_chr, ")", sep =""))+
        ylab("Normalized counts")+
        ggtitle(paste("Near bait analysis (", my_obj@bait_name, " bait)", sep = ""))
dev.off()


### Paired analysis, displaying rlog in the graph
pdf(file="./mapped_with_pipe4C/4Cker_out_rs7146599/dif_highP_Paired_rlog.pdf", width = 7, height = 4)
res_df = differentialAnalysis_rlog(obj=my_obj,
                        norm_counts_avg=nb_results$norm_counts_avg,
                        windows=nb_results$window_counts,
                        conditions=c("rs7146599_G", "rs7146599_A"),
                        region="nearbait",
                        coordinates=c("23980000", "24140000"),
			pval=0.05, 
                        baitloc = 24058463, 
                        oemin = 24094311, 
                        oemax = 24104125, 
                        plotStart = 24000000, 
                        plotEnd = 24130000, 
                        yMin = 0, 
                        yMax = 15)
dev.off() # 14 significnatly different regions
# I want to get the p-val though:
fwrite(res_df, file = "./mapped_with_pipe4C/4Cker_out_rs7146599/dif_highP_results_pvals.txt", sep = "\t", quote = F, row.names = F, col.names = T)

### in the paired analysis the p val drops at PCK2 (p = 0.04328212). But padj at start 24095353 = 0.08890382 
### we get sig result at the 3' end of PCK2

### Now plot per sample

M4G <- as.data.table(my_obj@data_nearbait[1])
M9G <- as.data.table(my_obj@data_nearbait[2])
M13G <- as.data.table(my_obj@data_nearbait[3]) 
M4A <- as.data.table(my_obj@data_nearbait[3])
M9A <- as.data.table(my_obj@data_nearbait[4])
M13A <- as.data.table(my_obj@data_nearbait[6])

# Region PCK2:
#chr14:24094311-24104125
mychr <- "chr14"
start <- 24094311
end <- 24104125


get_total_reads <- function(sample, start, end) {
        mygene <- sample[V1 == mychr & V3 < end & V3 > start | V1 == mychr & V2 > start & V2 < end]
        mygene[, V5 := sum(V4)] # summing across the whole of the gene
        mycount <- unique(mygene$V5)
        return(mycount)
}

M4G_PCK2 <- get_total_reads(M4G, start, end)
M9G_PCK2 <- get_total_reads(M9G, start, end)
M13G_PCK2 <- get_total_reads(M13G, start, end)
M4A_PCK2 <- get_total_reads(M4A, start, end)
M9A_PCK2 <- get_total_reads(M9A, start, end)
M13A_PCK2 <- get_total_reads(M13A, start, end)

### make a graph
toplot <- data.table(sample = c("M4", "M4", "M9", "M9", "M13", "M13"),
                     Normalised_reads_across_PCK2 = c(M4G_PCK2, M4A_PCK2, M9G_PCK2, M9A_PCK2, M13G_PCK2, M13A_PCK2),
                     Allele = c("G", "A", "G", "A", "G", "A"))

pdf(file = "./mapped_with_pipe4C/4Cker_out_rs7146599/PCK2_alleles_sample_comparison.pdf")
p <- ggplot(toplot, aes(x = sample, y = Normalised_reads_across_PCK2, fill = Allele)) +
        geom_bar(position="dodge", stat="identity") +
       theme(text = element_text(size = 15))
print(p)
dev.off()


