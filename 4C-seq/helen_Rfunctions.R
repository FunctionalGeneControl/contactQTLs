### Below code is modified from 4Cker tool:
https://github.com/rr1859/R.4Cker

### 4Cker function for finding differential interactions with a paired analysis, plotting the results with a green rectangle at the gene (or OE region)

differentialAnalysis_helen <- function(obj,norm_counts_avg, windows,conditions, region,coordinates=NULL, pval, baitloc, oemin, oemax, plotStart, plotEnd){
  if(.Platform$OS.type=="windows"){
    quartz<-function() windows()
  }
  pval_options=c(0.01, 0.05,0.1)
  if(length(conditions) != 2)
    stop("Only 2 conditions can be analyzed")
  
  colData_df = data.frame(unlist(lapply(1:length(obj@conditions), function(j) c(rep(obj@conditions[j], obj@replicates[j])))))
  
  rownames(colData_df) = colnames(windows)[-c(1:4)]
  colnames(colData_df) = c("condition")
  # add a subject column for paired analysis
  colData_df$subject <- unlist(lapply(1:length(obj@samples), function(j) sub("_[^_]+$", "", obj@samples[j])))
  
  condition = factor(unlist(lapply(1:length(obj@conditions), function(j) c(rep(obj@conditions[j], obj@replicates[j])))))
  # add subject here too
  subject <- factor(unlist(lapply(1:length(obj@samples), function(j) sub("_[^_]+$", "", obj@samples[j]))))

  #DESeq steps
  ### MAKING A PAIRED ANALYSIS
  dds = DESeqDataSetFromMatrix(countData=windows[, -c(1:4)],
                               colData = colData_df,
                               design = ~ subject + condition) ### Added subject here.
  
  dds <- estimateSizeFactors(dds)
  dds <- estimateDispersions(dds, fitType = "local")
  dds <- nbinomWaldTest(dds)
  res <- results(dds, c("condition", conditions[1], conditions[2]))
  if(length(which(res$padj < pval)) == 0){
    print("No significant changes..trying a higher p-value...")
    pval_row=which(pval_options == pval)
    while(!length(which(res$padj < pval)) == 0 & pval_row <=3){
      pval_row=pval_row+1
    }
    if(length(which(res$padj < pval_options[pval_row])) == 0)
      stop("No significant interactions")
    else{
      pval=pval_options[pval_row]
      print("Using ",pval , "...")
    }
  }
  norm_counts = counts(dds, normalized = TRUE)
  norm_counts_log = log(norm_counts+1,10)
  condition1_row = which(obj@conditions == conditions[1])
  condition2_row = which(obj@conditions == conditions[2])
  cols_conditions=NULL
  j=1
  for(i in obj@replicates){
    cols_conditions = rbind(cols_conditions, c(j, (j+i-1)))
    j=j+i
  }
  sig_rows = rep("not_sig", nrow(windows))
  sig_rows[which(res$padj < pval)] = "sig"
  if(region == "nearbait"){
    plot_df = data.frame(coord=c(rowMeans(windows[,2:3]),rowMeans(windows[,2:3])),
                         counts=c(rowMeans(norm_counts[,cols_conditions[condition1_row,1]:cols_conditions[condition1_row,2]]),
                                  rowMeans(norm_counts[,cols_conditions[condition2_row,1]:cols_conditions[condition2_row,2]])),
                         conditions=c(rep(conditions[1], nrow(windows)), rep(conditions[2], nrow(windows))),
                         sig=c(sig_rows, sig_rows))
    if(!is.null(coordinates)){
      plot_df=plot_df[which(plot_df[,1] >= coordinates[1] & plot_df[,1] <= coordinates[2]),]
    }
    
    # gene region
    rect <- data.frame(xmin=oemin, xmax=oemax, ymin=-Inf, ymax=Inf)

    print(ggplot(plot_df, aes(x=coord, y=counts, colour=conditions))+geom_line()+
      theme_bw()+
      xlab(paste("Chromosome coordinates (", obj@bait_chr, ")", sep =""))+
      ylab("Normalised counts")+geom_point(data=subset(plot_df,sig=="not_sig"), shape=1, size=0.5)+
      geom_point(data=subset(plot_df,sig=="sig")) + 
      # range 
      xlim(plotStart,plotEnd) +
      # add the bait (SNP)
      geom_vline(xintercept = baitloc, linetype = "dashed") +
      # add the gene region
      geom_rect(data=rect, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
              color="transparent", fill = "seagreen",
              alpha=0.3,
              inherit.aes = FALSE))
  }
  if(region == "cis"){
    plot_df = data.frame(coord=c(rowMeans(windows[,2:3]),rowMeans(windows[,2:3])),
                         counts=c(rowMeans(norm_counts_log[,cols_conditions[condition1_row,1]:cols_conditions[condition1_row,2]]),
                                  rowMeans(norm_counts_log[,cols_conditions[condition2_row,1]:cols_conditions[condition2_row,2]])),
                         conditions=c(rep(conditions[1], nrow(windows)), rep(conditions[2], nrow(windows))),
                         sig=c(sig_rows, sig_rows))
    print(ggplot(plot_df, aes(x=coord, y=counts, colour=conditions))+
      theme_bw()+
      xlab(paste("Chromosome coordinates (", obj@bait_chr, ")", sep =""))+
      ylab("Normalized counts")+geom_point(data=subset(plot_df,sig=="not_sig"), shape=1, size=0.5)+
      geom_vline(xintercept = baitloc, linetype = "dashed") +
      geom_vline(xintercept = oeloc, linetype = "dashed", color = "gray40") +
      geom_point(data=subset(plot_df,sig=="sig")))
  }
  sig_merge_windows=merge_windows(windows[which(res$padj < pval), 1:3])
  print(paste("BED file of significant domains saved in ", obj@output_dir, sep = ""))
  write.table(sig_merge_windows,
    paste(obj@output_dir, obj@bait_name, "_", conditions[1], "_", conditions[2], "_", region, "_pval", pval,"_diff.bed", sep = ""),
    quote=FALSE, col.names=FALSE, row.names=FALSE, sep = "\t")
  return(cbind(windows[,1:3],data.frame(res)))
}


#### Can we use different normalisation methods for visualisation?
## first trying rlog from deSeq2
## note that this doesnt affect the differential analysis, which is always run on the raw counts.
## The normalisation looks much better between samples on the plot.
differentialAnalysis_rlog <- function(obj,norm_counts_avg, windows,conditions, region,coordinates=NULL, pval, baitloc, oemin, oemax, plotStart, plotEnd, yMin = 0, yMax){
  if(.Platform$OS.type=="windows"){
    quartz<-function() windows()
  }
  pval_options=c(0.01, 0.05,0.1)
  if(length(conditions) != 2)
    stop("Only 2 conditions can be analyzed")
  
  colData_df = data.frame(unlist(lapply(1:length(obj@conditions), function(j) c(rep(obj@conditions[j], obj@replicates[j])))))
  
  rownames(colData_df) = colnames(windows)[-c(1:4)]
  colnames(colData_df) = c("condition")
  # add a subject column for paired analysis
  colData_df$subject <- unlist(lapply(1:length(obj@samples), function(j) sub("_[^_]+$", "", obj@samples[j])))
  
  condition = factor(unlist(lapply(1:length(obj@conditions), function(j) c(rep(obj@conditions[j], obj@replicates[j])))))
  # add subject here too
  subject <- factor(unlist(lapply(1:length(obj@samples), function(j) sub("_[^_]+$", "", obj@samples[j]))))

  #DESeq steps
  ### MAKING A PAIRED ANALYSIS
  dds = DESeqDataSetFromMatrix(countData=windows[, -c(1:4)],
                               colData = colData_df,
                               design = ~ subject + condition) ### Added subject here.
  
  dds <- estimateSizeFactors(dds)
  dds <- estimateDispersions(dds, fitType = "local")
  dds <- nbinomWaldTest(dds)
  res <- results(dds, c("condition", conditions[1], conditions[2]))
  if(length(which(res$padj < pval)) == 0){
    print("No significant changes..trying a higher p-value...")
    pval_row=which(pval_options == pval)
    while(!length(which(res$padj < pval)) == 0 & pval_row <=3){
      pval_row=pval_row+1
    }
    if(length(which(res$padj < pval_options[pval_row])) == 0)
      stop("No significant interactions")
    else{
      pval=pval_options[pval_row]
      print("Using ",pval , "...")
    }
  }
  norm_counts = counts(dds, normalized = TRUE)
  norm_counts_log = log(norm_counts+1,10)
  ## I have added the below; the rlog_counts then replaces norm_counts in the "nearbait" analysis below
  ## Note that the cis analysis always plotted the log normalised counts.
  non_norm_counts = counts(dds, normalized = FALSE)
  rlog_counts = rlog(non_norm_counts)
  ##

  condition1_row = which(obj@conditions == conditions[1])
  condition2_row = which(obj@conditions == conditions[2])
  cols_conditions=NULL
  j=1
  for(i in obj@replicates){
    cols_conditions = rbind(cols_conditions, c(j, (j+i-1)))
    j=j+i
  }
  sig_rows = rep("not_sig", nrow(windows))
  sig_rows[which(res$padj < pval)] = "sig"
  if(region == "nearbait"){
    plot_df = data.frame(coord=c(rowMeans(windows[,2:3]),rowMeans(windows[,2:3])),
                         counts=c(rowMeans(rlog_counts[,cols_conditions[condition1_row,1]:cols_conditions[condition1_row,2]]), # replaced norm_counts with rlog_counts
                                  rowMeans(rlog_counts[,cols_conditions[condition2_row,1]:cols_conditions[condition2_row,2]])), # replaced norm_counts with rlog_counts
                         conditions=c(rep(conditions[1], nrow(windows)), rep(conditions[2], nrow(windows))),
                         sig=c(sig_rows, sig_rows))
    if(!is.null(coordinates)){
      plot_df=plot_df[which(plot_df[,1] >= coordinates[1] & plot_df[,1] <= coordinates[2]),]
    }
    
    # gene region
    rect <- data.frame(xmin=oemin, xmax=oemax, ymin=-Inf, ymax=Inf)

    print(ggplot(plot_df, aes(x=coord, y=counts, colour=conditions))+geom_line(size = 0.5, alpha = 0.7)+
      theme_bw()+
      xlab(paste("Chromosome coordinates (", obj@bait_chr, ")", sep =""))+
      ylab("rlog normalised \ncounts")+
      #geom_point(aes(shape = conditions), size = 0.7) + 
      geom_point(data=subset(plot_df,sig=="not_sig"), size=0.3)+
      #scale_shape_manual(values=c(1,2)) +
      geom_point(data=subset(plot_df,sig=="sig")) +
      labs(colour = "Condition") + #, shape = "Condition") +
      # range 
      xlim(plotStart,plotEnd) +
      ylim(yMin, yMax) +
      # add the bait (SNP)
      geom_vline(xintercept = baitloc, linetype = "dashed") +
      # add the gene region
      geom_rect(data=rect, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
              color="transparent", fill = "seagreen",
              alpha=0.2,
              inherit.aes = FALSE) + 
      # theme to adjust labels etc
      theme(axis.title.y = element_text(angle = 0, vjust = 0.5)))
  }

  if(region == "cis"){
    plot_df = data.frame(coord=c(rowMeans(windows[,2:3]),rowMeans(windows[,2:3])),
                         counts=c(rowMeans(norm_counts_log[,cols_conditions[condition1_row,1]:cols_conditions[condition1_row,2]]),
                                  rowMeans(norm_counts_log[,cols_conditions[condition2_row,1]:cols_conditions[condition2_row,2]])),
                         conditions=c(rep(conditions[1], nrow(windows)), rep(conditions[2], nrow(windows))),
                         sig=c(sig_rows, sig_rows))
    print(ggplot(plot_df, aes(x=coord, y=counts, colour=conditions))+
      theme_bw()+
      xlab(paste("Chromosome coordinates (", obj@bait_chr, ")", sep =""))+
      ylab("Normalized counts")+geom_point(data=subset(plot_df,sig=="not_sig"), shape=1, size=0.5)+
      geom_vline(xintercept = baitloc, linetype = "dashed") +
      geom_vline(xintercept = oeloc, linetype = "dashed", color = "gray40") +
      geom_point(data=subset(plot_df,sig=="sig")))
  }
  sig_merge_windows=merge_windows(windows[which(res$padj < pval), 1:3])
  print(paste("BED file of significant domains saved in ", obj@output_dir, sep = ""))
  write.table(sig_merge_windows,
    paste(obj@output_dir, obj@bait_name, "_", conditions[1], "_", conditions[2], "_", region, "_pval", pval,"_diff.bed", sep = ""),
    quote=FALSE, col.names=FALSE, row.names=FALSE, sep = "\t")
  print(paste("File of normalised counts saved in ", obj@output_dir, sep = ""))
  write.table(plot_df,
	     paste(obj@output_dir, obj@bait_name, "_", conditions[1], "_", conditions[2], "_", region, "_rlog_norm_counts.txt", sep = ""), 
	     quote=FALSE, col.names=FALSE, row.names=FALSE, sep = "\t")


  return(cbind(windows[,1:3],data.frame(res)))
}

