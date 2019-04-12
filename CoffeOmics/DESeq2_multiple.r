#!/usr/bin/env Rscript

# Performs Differential Expression analysis with DESeq2
# ARGS are the sample numbers to compare
# ex: 'script.r 1v2 1v3 10v11' will use DESeq2 to compare samples 1v2, 1v3 and 10v11

# Required packages
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("DESeq2", version = "3.8")
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
cat("Is this the folder containing the htseq_counts files? (y/n) ")
read_stdin <- file("stdin")
answer <- readLines(read_stdin, 1)
close(read_stdin)
if (answer!="y") {
        stop("Exiting...", call.=FALSE)
}

# loading libraries
cat("\nLoading libraries...\n")
suppressMessages(library("DESeq2"))
suppressMessages(library("pheatmap"))
suppressMessages(library("ReportingTools"))
suppressMessages(library("RColorBrewer"))
suppressMessages(library("gplots"))

DESeq2_compare_ab <- function(a, b){
	cat("Comparing samples", a, "and", b, "\n")

	# filter sample files to the ones defined in a and b
	ab <- paste(paste("_", a, "[A-C]", sep=""),
		paste("_", b, "[A-C]", sep=""), sep="|")
	sampleFiles <- grep(ab, list.files(directory), value=TRUE)
	cat("sampleFiles:\n")
	print(sampleFiles)

	# condition is a vs b
	sampleCondition <- c(a, a, a, b, b, b)  # always 3 replicates for this experiment
	cat("sampleCondition:\n")
	print(sampleCondition)

	sampleTable <- data.frame(sampleName = sampleFiles,
		fileName = sampleFiles,
		condition = sampleCondition)

	cat("Reading htseq_count files...\n")
	ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
	  directory = directory,
		design= ~ condition)

	# filter low reads count genes (<15, as per the original MGX analysis)
	keep <- rowSums(counts(ddsHTSeq)) >= 15
	ddsHTSeq <- ddsHTSeq[keep,]

	######
	# DE #
	######
	cat("\nSTARTING DE\n")
	dds <- suppressMessages(DESeq(ddsHTSeq))
	res <- results(dds)

	# shrinkage of effect size
	resShrink <- suppressMessages(lfcShrink(dds, coef=resultsNames(dds)[2], type="apeglm"))

	# output files
	cat("Outputing files:\n")
	a_vs_b_title <- paste(a, b, sep="vs")

	# heatmap 1
	cat("Heatmaps...\n")
	ntd <- normTransform(dds)
	select <- order(rowMeans(counts(dds, normalized=TRUE)), decreasing=TRUE)[1:20]
	df <- as.data.frame(colData(dds)[, c("condition", "sizeFactor")])
	png(paste(a_vs_b_title, "heatmap.png", sep="_"), width=960, height=960)
	pheatmap(assay(ntd)[select, ], cluster_rows=FALSE, show_rownames=FALSE,
	  cluster_cols=FALSE, annotation_col=df)
	dev.off()

	# heatmap 2
	png(paste(a_vs_b_title, "heatmap2.png", sep="_"), width=960, height=960)
	heatmap.2(assay(ntd)[select, ], trace='none', margins=c(3, 8),
	  labCol=substr(colnames(assay(ntd)[select, ]), start=9, stop=10))
	dev.off()

	# heatmap of the sample-to-sample distances
	cat("Heatmap of the sample-to-sample distances...\n")
	vsd <- vst(dds, blind=FALSE)
	sampleDists <- dist(t(assay(vsd)))
	sampleDistMatrix <- as.matrix(sampleDists)
	rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$sizeFactor, sep="-")
	colnames(sampleDistMatrix) <- NULL
	colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
	png(paste(a_vs_b_title, "heatmap_sample_to_sample.png", sep="_"),
		width=960, height=960)
	pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists,
		clustering_distance_cols=sampleDists, col=colors)
	dev.off()

	# cook's distance boxplot
	cat("Cook's boxplot...\n")
	png(paste(a_vs_b_title, "cooks_boxplot.png", sep="_"), width=960, height=960)
	par(mar=c(8,5,2,2))
	boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2)
	dev.off()

	# log2 fold changes plot
	cat("Log2 fold changes plot...\n")
	png(paste(a_vs_b_title, "log2_fold_changes.png", sep="_"),
		width=960, height=960)
	plotMA(res, main=a_vs_b_title, ylim=c(-10,10))
	dev.off()

	# shrunken log2 fold changes plot
	cat("Shrunken log2 fold changes plot...\n")
	png(paste(a_vs_b_title, "shrunken_log2_fold_changes.png", sep="_"),
		width=960, height=960)
	plotMA(resShrink, main=paste(a_vs_b_title, "shrunken"), ylim=c(-10,10))
	dev.off()

	# HTML Report
	cat("Html report...\n")
	htmlreport <- HTMLReport(shortName=paste(a_vs_b_title, 'html_report', sep="_"),
		title=paste(a_vs_b_title, 'DESeq2 html report', sep=" "),
		reportDirectory=paste("./", a_vs_b_title, "_html_report", sep=""))
	publish(dds, htmlreport, pvalueCutoff=0.05, annotation.db="org.Mm.eg.db",
		factor=sampleCondition, reportDir=paste("./", a_vs_b_title,
		"_html_report", sep=""))
	finish(htmlreport)

	# data
	cat("Data files...\n")
	write(mcols(res)$description, file=paste(a_vs_b_title,
		"log2_fold_changes_variables_and_tests.txt", sep="_"))
	write(mcols(resShrink)$description, file=paste(a_vs_b_title,
		"shrunken_log2_fold_changes_variables_and_tests.txt", sep="_"))
	write.csv(as.data.frame(res), file=paste(a_vs_b_title, "results.csv", sep="_"))

	# Finished
	cat("Finished DE of", a, "vs", b, "\n")
}


# for each comparison, run DESeq2
for (comparison in args) {
	# split args
	avb <- strsplit(comparison, "v")
	a <- avb[[1]][1]
	b <- avb[[1]][2]

	# a and b must be numbers as string
	# string comparison must be a < b, because linux orders files this way
	if (a >= b) {
		stop("a must be < b (string comparison)", call.=FALSE)
	}

	# run DESeq2
	DESeq2_compare_ab(a, b)
}
