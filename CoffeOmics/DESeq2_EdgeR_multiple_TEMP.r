#!/usr/bin/env Rscript

# Performs Differential Expression analysis with DESeq2 and EdgeR
# ARGS are the sample numbers to compare
# ex: 'script.r 1v2 1v3 10v11' will use DESeq2 and EdgeR to compare samples 1v2, 1v3 and 10v11

# Required packages
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("DESeq2", version = "3.8")
# BiocManager::install("edgeR", version = "3.8")
# BiocManager::install("apeglm", version = "3.8")  # shrinkage
# BiocManager::install("ReportingTools", version = "3.8")  # html report
# install.packages("pheatmap")

# get command line arguments
args <- commandArgs(trailingOnly=TRUE)

# check for arguments
if (length(args)==0) {
	stop("At least one argument must be supplied in the form <x>v<y>", call.=FALSE)
}

# folder containing all htseq_counts (should be current folder)
directory <- getwd()
cat("selected directory is:", directory, "\n")

# ask if folder is correct
cat("Is this the folder containing the htseq_counts files? (y/n): ")
read_stdin <- file("stdin")
answer <- readLines(read_stdin, 1)
close(read_stdin)
if (answer!="y") {
        stop("Exiting...", call.=FALSE)
}

# loading libraries
cat("Loading libraries...\n")
suppressMessages(library("edgeR"))  # edgeR must be loaded before DESeq2
suppressMessages(library("DESeq2"))
suppressMessages(library("pheatmap"))
suppressMessages(library("ReportingTools"))
suppressMessages(library("RColorBrewer"))
suppressMessages(library("gplots"))

compare_ab <- function(a, b){
	cat("\nComparing samples", a, "and", b, "\n")

	# filter sample files to the ones defined in a and b
	# expects 3 replicates for each sample
	a_grep <- paste("_", a, "[A-C]", sep="")
	b_grep <- paste("_", b, "[A-C]", sep="")
	ab_grep <- paste(a_grep, b_grep, sep="|")
	sampleFiles <- grep(ab_grep, list.files(directory), value=TRUE)
	# if sample a is at the end, reverse vector
	if (grepl(a_grep, sampleFiles[6])) {
	  sampleFiles <- rev(sampleFiles)
	}
	cat("sampleFiles:\n")
	print(sampleFiles)

	# condition is a vs b
	sampleCondition <- c(a, a, a, b, b, b)
	cat("sampleCondition:\n")
	print(sampleCondition)

	sampleTable <- data.frame(sampleName = sampleFiles,
		fileName = sampleFiles,
		condition = sampleCondition)

	cat("Reading htseq_count files...\n")
	deseq2_HTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
		directory = directory,
		design= ~ condition)
	edger_HTSeq <- readDGE(sampleFiles, group=sampleCondition)

	# filter low reads count genes (<15, as per the original MGX analysis)
	keep <- rowSums(counts(deseq2_HTSeq)) >= 15
	deseq2_HTSeq <- deseq2_HTSeq[keep,]
	keep <- filterByExpr(edger_HTSeq)
	edger_HTSeq <- edger_HTSeq[keep, , keep.lib.sizes=FALSE]

	######
	# DE #
	######
	cat("STARTING DE\n")
	a_vs_b_title <- paste(a, b, sep="vs")

	# de
	deseq2_de <- suppressMessages(DESeq(deseq2_HTSeq))
	deseq2_res <- results(deseq2_de)
	edger_HTSeq <- estimateTagwiseDisp(estimateCommonDisp(calcNormFactors(edger_HTSeq)))
	edger_de <- exactTest(edger_HTSeq)
	edger_de_fdr <- topTags(edger_de, n = "inf")
	edger_de_fdr$table <- edger_de_fdr$table[!rownames(edger_de_fdr$table) %in% 
	                                     c("__no_feature", "__ambiguous", "__too_low_aQual", "__not_aligned", "__alignment_not_unique"), ]
	
	# deseq2 and edger intersection
	png(paste(a_vs_b_title, "deseq2_vs_edger_venn.png", sep="_"), width=960, height=960)
	dev.control('enable')
	venn::venn(list(deseq2_res@rownames, rownames(edger_de_fdr)), snames=c('deseq2', 'edger'),
	     zcolor='style')
	dev.copy(postscript, paste(a_vs_b_title, "deseq2_vs_edger_venn.eps", sep="_"),
	         width=960, height=960, onefile=TRUE, horizontal=FALSE)
	dev.off()
	dev.off()
	
	# filter deseq2 genes to only be the intersection
	intersection_genes <- intersect(deseq2_res@rownames, rownames(edger_de_fdr))
	deseq2_de <- deseq2_de[intersection_genes]
	deseq2_res <- results(deseq2_de)

	# shrinkage of effect size
	resShrink <- suppressMessages(lfcShrink(deseq2_de, coef=resultsNames(deseq2_de)[2], type="apeglm"))

	################
	# Output Files #
	################
	cat("Outputing files:\n")

	# heatmap 1
	cat("Heatmaps...\n")
	ntd <- normTransform(deseq2_de)
	select <- order(rowMeans(counts(deseq2_de, normalized=TRUE)), decreasing=TRUE)[1:20]
	df <- as.data.frame(colData(deseq2_de)[, c("condition", "sizeFactor")])
	png(paste(a_vs_b_title, "heatmap.png", sep="_"), width=960, height=960)
	dev.control('enable')
	pheatmap(assay(ntd)[select, ], cluster_rows=FALSE, show_rownames=FALSE,
		cluster_cols=FALSE, annotation_col=df)
	dev.copy(postscript, paste(a_vs_b_title, "heatmap.eps", sep="_"), width=960, height=960,
		onefile=TRUE, horizontal=FALSE)
	dev.off()
	dev.off()

	# heatmap 2
	png(paste(a_vs_b_title, "heatmap2_20.png", sep="_"), width=960, height=960)
	dev.control('enable')
	heatmap.2(assay(ntd)[select, ], trace='none', margins=c(3, 8), dendrogram='column',
		labCol=substr(colnames(assay(ntd)[select, ]), start=9, stop=10),
		key.title='Color Key', key.xlab='Gene Counts', key.ylab='# Genes w/ count')
	dev.copy(postscript, paste(a_vs_b_title, "heatmap2_20.eps", sep="_"), width=960,
		height=960, onefile=TRUE, horizontal=FALSE)
	dev.off()
	dev.off()

	png(paste(a_vs_b_title, "heatmap2_all.png", sep="_"), width=960, height=960)
	dev.control('enable')
	heatmap.2(assay(ntd), trace='none', margins=c(3, 1), dendrogram='column',
		labCol=substr(colnames(assay(ntd)), start=9, stop=10), labRow='',
		key.title='Color Key', key.xlab='Gene Counts', key.ylab='# Genes w/ count')
	dev.copy(postscript, paste(a_vs_b_title, "heatmap2_all.eps", sep="_"), width=960,
		height=960, onefile=TRUE, horizontal=FALSE)
	dev.off()
	dev.off()

	# heatmap of the sample-to-sample distances
	cat("Heatmap of the sample-to-sample distances...\n")
	vsd <- vst(deseq2_de, blind=FALSE)
	sampleDists <- dist(t(assay(vsd)))
	sampleDistMatrix <- as.matrix(sampleDists)
	rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$sizeFactor, sep="-")
	colnames(sampleDistMatrix) <- NULL
	colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
	png(paste(a_vs_b_title, "heatmap_sample_to_sample.png", sep="_"), width=960, height=960)
	dev.control('enable')
	pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists,
		clustering_distance_cols=sampleDists, col=colors)
	dev.copy(postscript, paste(a_vs_b_title, "heatmap_sample_to_sample.eps", sep="_"),
		width=960, height=960, onefile=TRUE, horizontal=FALSE)
	dev.off()
	dev.off()

	# cook's distance boxplot
	cat("Cook's boxplot...\n")
	png(paste(a_vs_b_title, "cooks_boxplot.png", sep="_"), width=960, height=960)
	dev.control('enable')
	par(mar=c(8,5,2,2))
	boxplot(log10(assays(deseq2_de)[["cooks"]]), range=0, las=2)
	dev.copy(postscript, paste(a_vs_b_title, "cooks_boxplot.eps", sep="_"), width=960,
		height=960, onefile=TRUE, horizontal=FALSE)
	dev.off()
	dev.off()

	# log2 fold changes plot
	cat("Log2 fold changes plot...\n")
	png(paste(a_vs_b_title, "log2_fold_changes.png", sep="_"),
		width=960, height=960)
	dev.control('enable')
	plotMA(deseq2_res, main=a_vs_b_title, ylim=c(-10,10))
	dev.copy(postscript, paste(a_vs_b_title, "log2_fold_changes.eps", sep="_"), width=960,
		height=960, onefile=TRUE, horizontal=FALSE)
	dev.off()
	dev.off()

	# log2 fold changes box plot
	cat("Log2 fold changes box plot...\n")
	png(paste(a_vs_b_title, "log2_fold_changes_boxplot.png", sep="_"),
	    width=960, height=960)
	dev.control('enable')
	boxplot_fig <- deseq2_res[["log2FoldChange"]]
	boxplot(deseq2_res[["log2FoldChange"]], names=a_vs_b_title)
	dev.copy(postscript, paste(a_vs_b_title, "log2_fold_changes_boxplot.eps", sep="_"),
		width=960, height=960, onefile=TRUE, horizontal=FALSE)
	dev.off()
	dev.off()

	# shrunken log2 fold changes plot
	cat("Shrunken log2 fold changes plot...\n")
	png(paste(a_vs_b_title, "shrunken_log2_fold_changes.png", sep="_"),
		width=960, height=960)
	dev.control('enable')
	plotMA(resShrink, main=paste(a_vs_b_title, "shrunken"), ylim=c(-10,10))
	dev.copy(postscript, paste(a_vs_b_title, "shrunken_log2_fold_changes.eps", sep="_"),
		width=960, height=960, onefile=TRUE, horizontal=FALSE)
	dev.off()
	dev.off()

	# HTML Report
	# cat("Html report...\n")
	# htmlreport <- HTMLReport(shortName=paste(a_vs_b_title, 'html_report', sep="_"),
	# 	title=paste(a_vs_b_title, 'DESeq2 html report', sep=" "),
	# 	reportDirectory=paste("./", a_vs_b_title, "_html_report", sep=""))
	# publish(deseq2_de, htmlreport, pvalueCutoff=0.05, annotation.db="org.Mm.eg.db",
	# 	factor=sampleCondition, reportDir=paste("./", a_vs_b_title,
	# 	"_html_report", sep=""))
	# finish(htmlreport)

	# data
	cat("Data files...\n")
	write(mcols(deseq2_res)$description, file=paste(a_vs_b_title,
		"log2_fold_changes_variables_and_tests.txt", sep="_"))
	write(mcols(resShrink)$description, file=paste(a_vs_b_title,
		"shrunken_log2_fold_changes_variables_and_tests.txt", sep="_"))
	write.csv(as.data.frame(deseq2_res), file=paste(a_vs_b_title, "results.csv", sep="_"))

	# Finished
	cat("Finished DE of", a, "vs", b, "\n")
	cat("------------------------------------\n")
	
	return(c(boxplot_fig, a_vs_b_title))
}

foldchange_boxplots = list()
a_vs_b_titles = character(length = 0)

# for each comparison, run DESeq2
for (comparison in args) {
    # split args
    avb <- strsplit(comparison, "v")
    a <- avb[[1]][1]
    b <- avb[[1]][2]
    
    # a and b must be numbers as string
    if (suppressWarnings(is.na(as.numeric(a))) || suppressWarnings(is.na(as.numeric(b)))) {
        stop("a and b must be numbers", call.=FALSE)
    }
    
    # run DESeq2
    multi_figs <- compare_ab(a, b)
    
    # build multi figures vetors
    foldchange_boxplots = c(foldchange_boxplots, list(as.numeric(head(multi_figs, n=-1))))
    a_vs_b_titles = c(a_vs_b_titles, tail(multi_figs, n=1))
}

# build multi figures
cat("Building multi figures...", "\n")

# log2_fold_changes_boxplot
png("log2_fold_changes_boxplot.png", width=960, height=960)
dev.control('enable')
do.call(boxplot, list(foldchange_boxplots, names=a_vs_b_titles))
dev.copy(postscript, "log2_fold_changes_boxplot.eps", width=960, height=960, onefile=TRUE,
         horizontal=FALSE)
dev.off()
dev.off()
