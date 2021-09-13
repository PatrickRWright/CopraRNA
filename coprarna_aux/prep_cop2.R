# prepares input for CopraRNA2

args <- commandArgs(trailingOnly = TRUE) 
inputFile <- args[1] # opt_tags.clustered

# read clusters
data <- read.table(inputFile, sep=";", header=TRUE)

prep_lines <- c()

for (i in unique(data$clusternumber)) { 
    curr_cluster<-data[which(data$clusternumber==i),]
    curr_cluster <- curr_cluster[order(curr_cluster$id2),] 
    curr_cluster_string <- paste(c(0,as.character(curr_cluster$id1),as.character('')), collapse=";")
    prep_lines <- c(prep_lines, curr_cluster_string)
}

write.table(file="CopraRNA2_prep.csv", prep_lines, row.names=F, quote=F, col.names=F)

