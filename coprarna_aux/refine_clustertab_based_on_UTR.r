

# R --slave -f ../refine_clustertab_based_on_UTR.r 




require(seqinr)


count_char <- function(string, char) {
	string_as_vector = unlist(strsplit(string, ""))
	char_counts = table(string_as_vector)
	char_counts[names(char_counts) == char]
}


mafft<-function(filename="ncrna.fa", outname="ncrna_aligned.fa"){
	command<-paste("mafft --maxiterate 1000 --retree 1 --localpair --quiet --inputorder ", filename, " > ", outname, sep="" )
	system(command)
	#fas<-read.fasta(outname)
	#fas
} 


proc_cdhit<-function(){ #x= pasted sorted cdhit cluster file

	cd<-read.delim("cdhit_res.txt.clstr", header=F, sep="?")
  cd<-as.character(cd[,1])
  cd<-gsub("\t"," ", cd)
  x<-cd
  clustlist<-list()
  numb<-grep(">Cluster", x)
  
  for(i in 1:length(numb)){
    if(i<length(numb)){
      end<-numb[i+1]-1
    }
    if(i==length(numb)){
      end<-length(x)
    }
    
    temp<-x[(numb[i]+1):end]
    temp<-gsub(".*aa, >","",temp)
    temp<-gsub("\\.\\.\\..*","",temp)
    clustlist[[i]]<-temp
    
  }
  clustlist	
}

#ooi<-colnames(clus)[2]


d<-dir()
fasta<-grep("_upfromstartpos_200_down_100.fa$", d)

fasta1<-c()

for(i in 1:length(fasta)){
fasta1<-c(fasta1, read.fasta(d[fasta[i]]))

}


clus<-read.csv("cluster.tab", sep="\t")
clus<-as.matrix(clus)
clus_new<-clus
#temp<-clus[3013,2:ncol(clus)]
ooi<-colnames(clus)[2]
id1<-gsub("\\|.*","",clus[,1])
id2<-gsub(".*\\|","",clus[,1])
id2<-gsub(" ","",id2)
id1<-as.numeric(id1)
id2<-as.numeric(id2)

id1<-max(id1)
id2<-max(id2)
#seq_length<-300
for(j in 1:nrow(clus)){
	print(paste(j,nrow(clus),sep="/"))
	tmp<-clus[j,2:ncol(clus)]
	if(tmp[1]!=""){
		#empty<-length(which(tmp==""))
		homologs<-count_char(paste(tmp[1],collapse=""),":")
		if(homologs>1){  # sub-clustering only if more than one homologs exist in the ooi
			empty<-which(tmp=="")
			if(length(empty)>0){
			tmp<-tmp[-empty]
			}
			
			locus<-c()
			
			name_file1<-gsub("\\:.*","",names(locus))
			#name_file2<-gsub(".*:","",locus)
			n3<-rep(seq_length,length(locus))
			name_file<-cbind(name_file1,locus,n3)
			write.table(name_file, file="names.txt", sep=" ", quote=F, row.names=F, col.names=F)
				for(i in 1:length(tmp)){
				
					loc<-strsplit(as.character(tmp[i])," ")[[1]]
					loc<-gsub(".*\\:","",loc)
					names(loc)<-paste(colnames(clus)[i+1],loc,sep=":")
					loc<-gsub("\\(.*\\)","",loc)
					locus<-c(locus,loc)
					
				}

				pos<-match(locus,names(fasta1))
				
				write.fasta(fasta1[pos],names=names(locus), file.out="temp.fasta")
				#write.fasta(fasta1[pos],names=locus, file.out="temp.fasta")
				fasta2<-fasta1[pos]
				names(fasta2)<-names(locus)
				 mafft(filename="temp.fasta", outname="temp2.fasta")
				  #command<-paste("clustalo -i ", "temp2.fasta", " --distmat-out=distmatout.txt --full --output-order=input-order --percent-id --force --max-hmm-iterations=-1", sep="")
				  #system(command)
				 dat<-read.phyDat("temp2.fasta", format="fasta", type="DNA")
				 dis <- dist.ml(dat, model="F81")
				dis<-as.matrix(dis)
				
				# # out_dis<-matrix(,nrow(dis)*(nrow(dis)-1)*0.5,7)
				# out_dis[,3]<-1
				# out_dis[,4]<-seq_length
				# out_dis[,5]<-1
				# out_dis[,6]<-seq_length
				
				# dis2<-upper.tri(dis)
				# count<-1
				# for(i in 1:nrow(dis)){
					# for(ii in 1:nrow(dis)){
						# if(dis2[i,ii]==T){
							# out_dis[count,1]<-colnames(dis)[i]
							# out_dis[count,2]<-rownames(dis)[ii]
							# out_dis[count,7]<-dis[ii,i]
							# count<-1+count
						# }
					# }
				# }
				# write.table(out_dis, file="dis.txt", sep=" ", quote=F, row.names=F, col.names=F)
				
				
				# system("domclust -S -o5 dis.txt names.txt")
				
				 # dis2<-read.delim("distmatout.txt",sep="",header=F, , skip=1)
				  # na<-dis2[,1]
				 # unlink("distmatout.txt")
				  # unlink("temp_fasta")
				  # dis2<-dis2[,2:ncol(dis2)]
				  # colnames(dis2)<-na
				  # rownames(dis2)<-na
				 # dis3<-as.matrix(dis2)
					
				# knum<-max(2,length(apcluster(dis)))
				# ooi_clus<-grep(ooi, colnames(dis))
				
				### cdhit solution ###########
				system("cd-hit -n 4 -i temp.fasta -o cdhit_res.txt -c 0.60 -d 1000 -g 1")
				cd<-proc_cdhit()
				final_clus<-cd[grep(ooi,cd)]   # select cluster with member from ooi
				
				rest<-unlist(cd)
				pos<-unlist(final_clus)
				rest<-setdiff(rest, pos)
				
				ooi_pos<-grep(ooi, colnames(dis))
				if(length(rest)>0){
					for(i in 1:length(rest)){ 
					#	print(i)					# the remaining UTRS were put into the cluster with the ooi homolog with the lowest distance to the respective UTR
						tmp1<-grep(gsub("\\(.*","",rest[i]), colnames(dis))
						tmp1<-names(which(dis[ooi_pos,tmp1]==min(dis[ooi_pos,tmp1])))
						tmp1<-grep(gsub("\\(.*","",tmp1), final_clus)
						final_clus[[tmp1]]<-c(final_clus[[tmp1]],rest[i])
					}
				}
				# for(i in 1:length(cd)){
					# tmp<-cd[[i]]
					# #tmp<-gsub(".*\\:","",tmp)
					# #tmp<-gsub("\\(.*","",tmp)
					# tmp<-match(tmp, names(fasta2))
					# write.fasta(fasta2[tmp], names=names(fasta2)[tmp], file.out=paste(i, "test.fasta", sep="_"))
				# }
				
				# ##############
				# tmp_clus<-cmeans(dis,knum)
				# tmp_cluster<-tmp_clus$cluster
				# tmp_clus<-tmp_clus$membership

				# diff_clus<-unique(tmp_cluster[which(gsub("\\:.*","",names(tmp_cluster))==ooi)])
				# final_clus<-list()
				# knum<-length(grep(ooi, colnames(dis)))
				
				# if(knum==1){
					# dis2<-sort(dis[,1])
					# t3<-names(dis2)
					# t4<-which(duplicated(t3))
					# if(length(t4)>0){
						# t3<-t3[-t4]
					# }
					# final_clus[[1]]<-t3
				# }
				
				# if(knum>1){
					# if(length(diff_clus)>1){
					# clus1<-c()
					# for(i in 1:knum){
						# t1<-match(colnames(dis)[ooi_clus[i]], rownames(tmp_clus))
						# t1<-which(tmp_clus[t1,]==max(tmp_clus[t1,]))
						# clus1<-c(clus1,t1)
						# t2<-tmp_clus[order(tmp_clus[,t1],decreasing=T),]
						# t3<-rownames(t2)[which(t2[,t1]>=min(0.45,tmp_clus[colnames(dis)[ooi_clus[i]],t1]))]
						# t4<-gsub("\\:.*","",t3)
						# t4<-which(duplicated(t4))
						# if(length(t4)>0){
							# t3<-t3[-t4]
						# }
						# final_clus[[i]]<-t3
					# }
				
				# res<-setdiff(rownames(tmp_clus),unlist(final_clus))
				
				# for(i in 1:length(res)){
					# t1<-match(res[i], rownames(tmp_clus))
					# t1<-which(tmp_clus[t1,clus]==max(tmp_clus[t1,clus]))
					# t2<-final_clus[[t1]]
					# t2<-c(t2,res[i])
					# final_clus[[t1]]<-t2
				
				# }
				
				
				# le<-unlist(lapply(final_clus,length))
				# le<-which(le>=3)
				# final_clus<-final_clus[le]
				# if(length(le>0)){
					# for(i in 1:length(final_clus)){
						# if(i==1){
							# if(length(final_clus[[i]])>=3){
								# na<-gsub("\\:.*","",final_clus[[i]])
								# na<-match(na, colnames(clus))
								# clus_new[j,2:ncol(clus_new)]<-""
								# clus_new[j,na]<-final_clus[[i]]
							# }
						# }
				for(i in 1:length(final_clus)){		
						#if(length(final_clus[[i]])>=3){
							na<-gsub("\\:.*","",final_clus[[i]])
							na<-match(na, colnames(clus))
							
							new<-rep("",ncol(clus))
							new[na]<-final_clus[[i]]
							id1<-id1+1
							id2<-id2+1
							new[1]<-paste(id1,id2,sep="|")
							clus_new<-rbind(clus_new,new)
						#}
						}
			}
		}
	}

write.table(clus_new,file="cluster.tab", sep="\t", quote=F, row.names=F)
write.table(clus,file="cluster_old.tab", sep="\t", quote=F, row.names=F)  


