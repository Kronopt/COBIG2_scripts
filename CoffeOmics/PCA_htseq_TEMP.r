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
           col.names=c("gene", substring(file_name, 9, 11)),
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
rownames(dataset) <- c("10A", "10B", "10C", "11A", "11B", "11C", "12A", "12B", "12C", "1A", "1B", "1C", "2A", "2B", "2C", "3A", "3B", "3C", "4A", "4B", "4C", "5A", "5B", "5C", "6A", "6B", "6C", "7A", "7B", "7C", "8A", "8B", "8C", "9A", "9B", "9C")


###
# pca
###

png("PCA_htseq_counts.png")
dev.control('enable')
pca <- prcomp(dataset)
colours <- c("#bd0000","#bd0000","#bd0000", "#bd5e00","#bd5e00","#bd5e00", "#bdbd00","#bdbd00","#bdbd00", "#5ebd00","#5ebd00","#5ebd00", "#00bd00","#00bd00","#00bd00", "#00bd5e","#00bd5e","#00bd5e", "#00bdbd","#00bdbd","#00bdbd", "#005ebd","#005ebd","#005ebd", "#0000bd","#0000bd","#0000bd", "#5e00bd","#5e00bd","#5e00bd", "#bd00bd","#bd00bd","#bd00bd", "#bd005e","#bd005e","#bd005e")
pca_graphic <- fviz_pca_ind(pca,
             invisible=c("quali"),
             col.ind=colours,
             repel=TRUE)
pca_graphic + theme(legend.position="none") + geom_point(aes(colour = colours))
dev.copy(postscript, "PCA_htseq_counts.eps", onefile=TRUE, horizontal=FALSE)
dev.off()
dev.off()
