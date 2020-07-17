#!/usr/bin/env Rscript

#dependencies:
# mafft

#call:
# Rscript --slave refine_clustertab.r 
# no arguments

suppressPackageStartupMessages(require(seqinr))
suppressPackageStartupMessages(require(doMC))


# register cores for parallel processing
co<-readLines("CopraRNA_option_file.txt") 
max_cores <- min( detectCores(), as.numeric(gsub("core count:","",co[grep("core count:", co)])), na.rm = TRUE)
registerDoMC(max_cores)


# read the organism of interest (ooi) from the ncRNA fasta file. The sRNA of the ooi is considered to be the first sequence.
ooi<-gsub("ncRNA_","",names(read.fasta("./ncrna.fa"))[1])  
print(paste("re-clustering for ooi",ooi))


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
## dm = distance matrix
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
write.fasta(fasta1, names=names(fasta1),file="utr_seqs.fa")

#read the original cluster file produced by domclust and write a backup version
clus<-read.csv("cluster.tab", sep="\t")
write.table(clus,file="cluster_backup.tab", sep="\t", quote=F, row.names=F)
clus<-as.matrix(clus)
clus_new<-clus
colnames(clus_new)[1]<-"#id"

# identify the pseudo kegg code assigned to the ooi
ooi_pos<-grep(paste(ooi,"_upfrom.*_.*fa$",sep=""), dir())
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

fasta_text<-function(fasta){
	out<-c()
	for(i in 1:length(fasta)){
		out<-c(out, paste(">",names(fasta)[i]))
		out<-c(out, paste(fasta[[i]],collapse=""))
	}
	out
}

parse_fasta<-function(x){
	fasta<-x
	seq_start<-grep(">", fasta)
	seq_name<-fasta[seq_start]
	seq_name<-gsub(">","",seq_name)
	seqs<-vector("list",length(seq_start))
	names(seqs)<-seq_name
	for(i in 1:length(seq_start)){
		if(i<length(seq_start)){
			temp<-fasta[(seq_start[i]+1):(seq_start[i+1]-1)]
			temp<-as.character(temp)
			temp<-gsub(" ","",temp)
			temp<-paste(temp,collapse="")
			temp<-strsplit(temp,"")[[1]]
			seqs[[i]]<-temp
		}
		if(i==length(seq_start)){
			temp<-fasta[(seq_start[i]+1):length(fasta)]
			temp<-as.character(temp)
			temp<-gsub(" ","",temp)
			temp<-paste(temp,collapse="")
			temp<-strsplit(temp,"")[[1]]
			seqs[[i]]<-temp
		}
	}
	seqs
}


dist_hamming<-function(x){
	se<-x
	le<-length(se)
	le_se<-length(se[[1]])
	out<-matrix(,le,le)
	out[]<-0
	colnames(out)<-names(se)
	rownames(out)<-names(se)
	for(i in 1:le){
		for(j in i:le){
			if(i!=j){
				temp<-length(which(se[[i]]!=se[[j]]))/(le_se-length(which(se[[i]]== "-" & se[[j]]=="-" )))
				out[i,j]=out[j,i]=temp
			}
		}
	}
	out
	
}

# check occurence of clusters with multiple ooi entries
mult<-c()
for(j in 1:nrow(clus)){
	tmp<-clus[j,2:ncol(clus)]
	# check if the cluster contains a gene from the ooi
	if(tmp[ooi]!=""){ 								
		homologs<-count_char(paste(tmp[ooi],collapse=""),":")
		# sub-clustering only if more than one homologs exist in the ooi
		if(homologs>1){  
			mult<-c(mult,j)
		}
	}
}

# deal with clusters with multiple ooi entries
if ( length(mult) > 0 ) {

clus_in<-clus
clus2<-clus[mult,]

# divide the data in subsets for parallel processing 
max_cores<-min(max_cores,length(mult))
jobs<-length(mult)%/%max_cores
rest<-length(mult)-max_cores*jobs
jobs<-rep(jobs,max_cores)
if(rest>0){
	jobs[1:rest]<-jobs[1:rest]+1
}

count_vect1<-cumsum(c(1,jobs[1:(length(jobs)-1)]))
count_vect2<-cumsum(jobs)


# generate a temp file per thread to prepare mafft input file
thread2tmpfile = c();
for (i in 1:max_cores) { thread2tmpfile = c(thread2tmpfile, tempfile(pattern="CopraRNA2.refineClusterTab.")); }

# start parallel processing
vari<-foreach(ji=1:max_cores)  %dopar% {
	clus<-clus2[count_vect1[ji]:count_vect2[ji],]
	if(is.vector(clus)){
		clus<-t(as.matrix(clus))
	}
	out_clus<-vector("list",nrow(clus))
	names(out_clus)<-clus[,1]
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
				temp_fasta<-fasta1[pos]
				temp_fasta<-fasta_text(temp_fasta)
				tempf<-thread2tmpfile[ji]
				write.fasta(fasta1[pos],names=names(locus), file.out=tempf)
				ca<-paste("mafft --maxiterate 1000 --localpair --quiet --inputorder ", tempf,"",sep="")
				align<-parse_fasta(system(ca, intern=T))
				dis<-dist_hamming(align)
				#  do the densitity based clustering
				clust_tmp<-cluster(dis,eps=0.40,1)
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
				out_clus[[j]]<-final_clus
			}
		}
	}
	out_clus
}

# cleanup temp files
file.remove(thread2tmpfile); 


final_clus1<-list()
for(i in 1:length(vari)){
	final_clus1<-c(final_clus1,vari[[i]])
}

clus<-clus_in	
# the cluster file is updated and extended based on the new cluster results
for(j in 1:length(final_clus1)){	
	final_clus<-final_clus1[[j]]
	clus_p<-mult[j]
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

} # if cluster with multiple ooi entries

write.table(clus_new,file="cluster.tab", sep="\t", quote=F, row.names=F)