# get input parameter
## edit 2.0.4 // changing to use of CopraRNA_option_file.txt
options <- read.table("CopraRNA_option_file.txt", sep=":")
relSize <- as.numeric(as.character(options$V2[5]))

# get the maximum cluster size
fasta <- read.table("16s_sequences.fa")
max_cluster_size <- length(grep(">", fasta$V1))

d <- read.table("opt_tags.clustered", sep=";", header=T) ## edit 2.0.4

keep_clusters <- c()
# for 1 to the biggest clusternumber // the statement in the brackets finds the maximum clusternumber    
for (i in 1:(d[length(d$d1),]$clusternumber)) { ## edit 2.0.5.1 // changed to d1 due to IntaRNA 2 output

    if(length(which(d$clusternumber==i))/max_cluster_size >= relSize) { 
        keep_clusters <- c(keep_clusters,i)
    }
}

filtered_tags_clustered <- d[which(d$clusternumber %in% keep_clusters),]

colnames(filtered_tags_clustered)[35]<-"p-value" # ## edit 2.0.5.1 // otherwise it calls it p.value

write.table(file="opt_tags.clustered_rcsize", filtered_tags_clustered, quote=F, row.names=F, sep=";")

