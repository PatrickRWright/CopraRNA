# script by Jens Georg

#dependencies
# clustalo
# mafft

#call:
# R --slave -f ./refine_clustertab.r 
# no arguments

print("re-clustering")
require(seqinr)

# read the organism of interest (ooi) from the ncRNA fasta file. The sRNA of the ooi is considered to be the first sequence.
ooi<-gsub("ncRNA_","",names(read.fasta("./ncrna.fa"))[1])  
print(ooi)


# function for 
count_char <- function(string, char) {
	string_as_vector = unlist(strsplit(string, ""))
	char_counts = table(string_as_vector)
	char_counts[names(char_counts) == char]
}

# wrapper for mafft
mafft<-function(filename="./ncrna.fa", outname="ncrna_aligned.fa"){
	command<-paste("mafft --maxiterate 1000 --retree 1 --localpair --quiet --inputorder ", filename, " > ", outname, sep="" )
	system(command)
} 


# function for density based clustering similar to DBSCAN
## dm = similarity matrix
## eps = minimal similarity for clustering
## minpts = minimal cluster size

cluster<-function(dm,eps,minpts){
	vis<-c()
	clus<-list()
	for(i in 1:ncol(dm)){
		cand<-colnames(dm)[i]
		temp<-c()
		while(is.na(match(cand,vis))){
				vis<-c(vis, cand)
				temp1<-which(dm[,cand]<=eps)
				if(length(temp1)>=minpts){
					temp1<-colnames(dm)[temp1]
					temp<-unique(c(temp,temp1))
					
					cand<-setdiff(temp,vis)[1]
					if(is.na(cand)){
						cand<-vis[1]
					}
				}
			}
		clus[[i]]<-temp
	}
	l<-unlist(lapply(clus, length))
	l<-which(l>0)
	clus<-clus[l]
	clus
}

# read the UTR fasta sequences from all organisms
d<-dir()
fasta<-grep("_upfrom.*_down_.*.fa$", d)
fasta1<-c()
for(i in 1:length(fasta)){
fasta1<-c(fasta1, read.fasta(d[fasta[i]]))

}

#read the original cluster file produced by domclust and write a backup version
clus<-read.csv("cluster.tab", sep="\t")
write.table(clus,file="cluster_backup.tab", sep="\t", quote=F, row.names=F)
clus<-as.matrix(clus)
clus_new<-clus
colnames(clus_new)[1]<-"#id"

# identify the pseudo kegg code assigned to the ooi
ooi_pos<-grep(paste(ooi,"_upfromstartpos_.*fa$",sep=""), dir())
inp<-read.fasta(dir()[ooi_pos])
ooi_pos<-names(inp)[1:100]
ooi2<-c()
for(i in 1:length(ooi_pos)){
	ooi2<-c(ooi2, grep(paste(":",ooi_pos[i],sep=""), clus, ignore.case=T))
}
ooi<-ooi2[1]
ooi<-gsub("\\:.*","",clus[ooi])


# identify numbering for existing clusters 
id1<-gsub("\\|.*","",clus[,1])
id2<-gsub(".*\\|","",clus[,1])
id2<-gsub(" ","",id2)
id1<-as.numeric(id1)
id2<-as.numeric(id2)
id1<-max(id1)
id2<-max(id2)

# check each cluster if re-clustering is neccessary (two or more genes of the ooi in one cluster) and perform re-clustering
for(j in 1:nrow(clus)){
	tmp<-clus[j,2:ncol(clus)]
	# check if the cluster contains a gene from the ooi
	if(tmp[ooi]!=""){ 								
		homologs<-count_char(paste(tmp[ooi],collapse=""),":")
		# sub-clustering only if more than one homologs exist in the ooi
		if(homologs>1){  							
			empty<-which(tmp=="")
			if(length(empty)>0){
			tmp<-tmp[-empty]
			}
			# extract locus tags and write the corresponding UTRs in a temporary fasta file
			locus<-c()
			for(i in 1:length(tmp)){
				loc<-strsplit(as.character(tmp[i])," ")[[1]]
				loc2<-loc
				loc<-gsub(".*\\:","",loc)
				names(loc)<-loc2
				loc<-gsub("\\(.*\\)","",loc)
				locus<-c(locus,loc)
			}
			dup<-which(duplicated(locus))
			if(length(dup)>0){
				locus<-locus[-dup]
			}
			ooi_locus_tags<-names(locus)[grep(ooi,names(locus))]
			pos<-match(locus,names(fasta1))
			write.fasta(fasta1[pos],names=names(locus), file.out="temp.fasta")
			# perform clustalo to get a fast distance matrix (all vs all UTRs in the cluster)
			command<-paste("clustalo -i ", "temp.fasta", " --distmat-out=distmatout.txt --full --percent-id --output-order=input-order --force --max-hmm-iterations=-1", sep="")
			system(command)
			dis<-read.delim("distmatout.txt",sep="",header=F, , skip=1)
			dis<-dis[,2:ncol(dis)]
			unlink("distmatout.txt")
			colnames(dis)<-names(locus)
			rownames(dis)<-names(locus)
			# transform the distance matrix in a similarity matrix and do the densitity based clustering
			clust_tmp<-cluster(1-dis/100,eps=0.40,1)
			# identify clusters that do not contain a gene from the ooi
			## merge the genes from the "no_ooi" clusters in a single vector and remove the clusters from the temporary cluster list
			no_ooi<-setdiff(1:length(clust_tmp),grep(ooi, clust_tmp))
			if(length(no_ooi)>0){
				no_ooi<-unlist(clust_tmp[no_ooi])
				clust_tmp<-clust_tmp[grep(ooi, clust_tmp)]
			} else{
				no_ooi<-c()
			}
			
			# create the new clusters with only one ooi gene per cluster
			final_clus<-list()
			for(jj in 1:length(clust_tmp)){
				tmp_ooi<-grep(ooi, clust_tmp[[jj]])
				
				# if a cluster still contains more than one gene from the ooi (e.g. because the UTRs are too similar), the cluster is multiplied by the number of ooi genes
				if(length(tmp_ooi)>1){
					nam<-clust_tmp[[jj]][tmp_ooi]
					for(jjj in 1:length(tmp_ooi)){
						tmp2<-c(setdiff(clust_tmp[[jj]],nam),nam[jjj])
						final_clus[[length(final_clus)+1]]<-tmp2
					}
				} else {
					final_clus[[length(final_clus)+1]]<-clust_tmp[[jj]]
				}
			}
			

			# the genes of the no_ooi clusters are assigned to all final clusters
			final_clus<-lapply(final_clus, function (x){return(c(x,no_ooi))})
			
			# the cluster file is updated and extended based on the new cluster results
			for(i in 1:length(final_clus)){					
				na1<-unique(gsub("\\:.*","",final_clus[[i]]))
				na<-match(na1, colnames(clus))
				new<-rep("",ncol(clus))
				te<-final_clus[[i]]
				te2<-c()
				for(jj in 1:length(na1)){
					p1<-grep(na1[jj],te)
					te2<-c(te2, paste(te[p1],collapse=" "))
				}
				new[na]<-te2	
				if(i==1){
					ooi_p<-grep(ooi, gsub("\\(.*","",final_clus[[i]]))
					clus_p<-j
					clus_new[clus_p,2:ncol(clus_new)]<-new[2:length(new)]
				}
				if(i>1){
					id1<-id1+1
					id2<-id2+1
					new[1]<-paste(id1,id2,sep="|")
					clus_new<-rbind(clus_new,new)
				}
			}
		}
	}
}
write.table(clus_new,file="cluster.tab", sep="\t", quote=F, row.names=F)