#!/usr/bin/env Rscript

#### Get the allele counts allowing for parallelised submission of samples. 

suppressMessages(library(argparser))
args = commandArgs(trailingOnly=T)

p <- arg_parser("Getting allele counts", name="Rscript get_allele_counts_for_parallel.R", hide.opts = TRUE)
p <- add_argument(p, arg="<bed>",
                  help="path to bed file")
p <- add_argument(p, arg="<peaks>",
                  help="path to peaks file in format Chr, hg38Proxy_zero, hg38Proxy_pos, hg19SNP_ID, PeakChr, Peak_start, Peak_end, Peak_name, Dist2Peak")
p = add_argument(p, arg="--outdir",
                 help="Full path to output directory", default=".")
p = add_argument(p, arg="--verbose",
                 help = "Flag specifying whether to print process steps.", flag=TRUE)

opts = parse_args(p, args)

bed = opts[["<bed>"]]
peaks=opts[["<peaks>"]]
outdir = opts[["outdir"]]
verb = opts[["verbose"]]

suppressMessages(library(data.table))
suppressMessages(library(dplyr))

########### 
filename <- basename(bed)
id <- ".bed"
myname <- sub(id, "", filename)

sink(file = paste0(outdir, "/", myname, ".report.txt"))

# Will need to intersect the snps with the peaks found in this bed file (see script 5_match_SNPs_with_peaks.sh; this bed file used "newIDs"):
alleles <- fread(peaks)
names(alleles) = c("Chr", "hg38Proxy_zero", "hg38Proxy_pos", "hg19SNP_ID", "PeakChr", "Peak_start", "Peak_end", "Peak_name", "Dist2Peak")

mycounts <- fread(bed)
names(mycounts) = c("Chr", "readStart", "readEnd", "readID", "qual", "strand")
#temp <- mycounts[1:1000000, ]

# Intersect the counts with the peaks
setkey(mycounts, Chr, readStart, readEnd)
with_peaks <- foverlaps(alleles, mycounts, by.x = c("PeakChr", "Peak_start", "Peak_end"), nomatch = NULL)

filename <- basename(bed)
id <- ".bed"
myname <- sub(id, "", filename)


with_peaks[, feature1 := paste(Peak_start, Peak_end, sep = "-")]
with_peaks[, feature := paste(PeakChr, feature1, sep = ":")]

#### Not yet tested

print("Counting SNPs")
to_tally <- unique(with_peaks[, .(hg19SNP_ID, feature, hg38Proxy_pos, readID)])

overall_counts <- as.data.table(to_tally %>% 
                                  group_by(feature, hg19SNP_ID, hg38Proxy_pos) %>% distinct() %>% tally())
overall_counts[, Sample := myname]
setnames(overall_counts, c("hg19SNP_ID", "hg38Proxy_pos"), c("SNP", "hg38SNP_pos"))


write.table(overall_counts, file = paste0(outdir, "/", myname, "_total_counts.txt"), 
              quote = FALSE, row.names = FALSE, sep = "\t")
  
  
