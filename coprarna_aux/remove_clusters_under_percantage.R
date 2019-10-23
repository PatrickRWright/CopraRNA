# get input parameter

options <- read.table("CopraRNA_option_file.txt", sep=":",colClasses = "character")
relSize <- as.numeric(options$V2[grep("relative clustersize", options$V1)])

# get the maximum cluster size
fasta <- read.table("16s_sequences.fa")
max_cluster_size <- length(grep(">", fasta$V1))

d <- read.table("opt_tags.clustered", sep=";", header=T)

keep_clusters <- c()
# for 1 to the biggest clusternumber // the statement in the brackets finds the maximum clusternumber    
for (i in 1:(d[length(d$d1),]$clusternumber)) {

    if(length(which(d$clusternumber==i))/max_cluster_size >= relSize) { 
        keep_clusters <- c(keep_clusters,i)
    }
}

filtered_tags_clustered <- d[which(d$clusternumber %in% keep_clusters),]

# reset "p-value" column header, otherwise it calls it p.value
colnames(filtered_tags_clustered)[grep("p.value",colnames(filtered_tags_clustered),ignore.case=TRUE,perl=TRUE)]<-"p-value" 

write.table(file="opt_tags.clustered_rcsize", filtered_tags_clustered, quote=F, row.names=F, sep=";")

