
#call:
# R --slave -f /home/jens/CopraRNA-git/coprarna_aux/remove_full_hybrids.r 

suppressPackageStartupMessages(require(seqinr))

co<-readLines("CopraRNA_option_file.txt")
hybrid_threshold<-as.numeric(gsub("hybrid_threshold:","",co[grep("hybrid_threshold:", co)]))


check_hybrid<-function(hybrid, srna_length){
	hybrid<-gsub(".*&","",hybrid)
	hybrid<-gsub("\\.","",hybrid)
	return(nchar(hybrid)/srna_length)
}

if(hybrid_threshold<1){
	removed<-c()
	d<-dir()
	sRNAs<-read.fasta("input_sRNA.fa")
	cluster<-read.csv("cluster.tab", sep="\t")
	for(i in 1:length(sRNAs)){
		le<-length(sRNAs[[i]])
		intaRNA<-grep(paste(names(sRNAs)[i],".*.fa.intarna.csv",sep=""),d) 
		intaRNA<-read.csv(d[intaRNA],sep=";")
		hybrids<-unlist(lapply(intaRNA[,"hybridDP"],check_hybrid, srna_length=le))
		long_hybrids<-which(hybrids>=hybrid_threshold)
		if(length(long_hybrids)>0){
			removed<-rbind(removed,intaRNA[long_hybrids,])
			tags<-intaRNA[long_hybrids,1]
			for(j in 1:length(tags)){
				pos_col<-grep(tags[j],cluster)
				pos_row<-grep(tags[j],cluster[,pos_col])
				pos_col<-rep(pos_col,length(pos_row))
				print(c(pos_col,pos_row,tags[j]))
				if(length(pos_row)>0){
					for(jj in 1:length(pos_row)){
						tmp<-cluster[pos_row[jj],pos_col[jj]]
						tmp<-strsplit(tmp," ")[[1]]
						pos<-grep(tags[j],tmp)
						tmp<-tmp[-pos]
						tmp<-paste(tmp,collapse=" ")
						cluster[pos_row[jj],pos_col[jj]]<-tmp
					}	
				}				
			}
		}		
	}
	write.table(cluster,file="cluster.tab", sep="\t", quote=F, row.names=F)
	write.table(removed,file="removed_full_hybrids.txt", sep="\t", quote=F, row.names=F)
} 
