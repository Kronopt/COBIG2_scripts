#!/usr/bin/env Rscript

# http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/118-principal-component-analysis-in-r-prcomp-vs-princomp/

# install.packages('factoextra')
library(factoextra)


htseq_counts_folder <- getwd()


###
# read htseq-counts files
###

read_htseq <- function(main_folder, file_name) {
  read.csv(paste(main_folder, "/", file_name, sep=''),
           header=FALSE,
           row.names=1,
           col.names=c("gene", substring(file_name, 9, 10)),
           sep='\t')
}

htseq_count_list <- list.files(path=htseq_counts_folder, pattern="*.htseq_counts*")
for (i in 1:length(htseq_count_list)) {
  if (!exists("dataset")){
    dataset <- read_htseq(htseq_counts_folder, htseq_count_list[i])
  }
  else {
    temp_dataset <- read_htseq(htseq_counts_folder, htseq_count_list[i])
    dataset <- cbind(dataset, temp_dataset)
    rm(temp_dataset)
  }
}

# remove statistics at the end of each hteq-counts file (last 5 lines)
dataset <- dataset[-seq(nrow(dataset) - 4, nrow(dataset)),]
dataset <- t(dataset)

# change row names
rownames(dataset) <- c("1A", "1B", "1C", "3A", "3B", "3C", "5A", "5B", "5C", "7A", "7B", "7C")


###
# pca
###

png("PCA_htseq_counts.png")
dev.control('enable')
pca <- prcomp(dataset)
pca_graphic <- fviz_pca_ind(pca,
             invisible=c("quali"),
             col.ind=c("#00AFBB","#00AFBB", "#00AFBB", "#00BB22","#00BB22","#00BB22", "#BBBB00","#BBBB00","#BBBB00", "#BB3000","#BB3000","#BB3000"),
             repel=TRUE)
pca_graphic + theme(legend.position="none")
dev.copy(postscript, "PCA_htseq_counts.eps", onefile=TRUE, horizontal=FALSE)
dev.off()
dev.off()
