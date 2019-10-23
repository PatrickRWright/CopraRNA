## edit 1.2.8 new script
## script by Jens Georg

args <- commandArgs(trailingOnly = TRUE) ## edit 2.0.5.1
inputFile_CopraRNA <- args[1] ## edit 2.0.5.1
inputFile_DAVID <- args[2] ## edit 2.0.5.1
output_file <- args[3] ## edit 2.0.5.1


# number of top predictions which should be investigated
co<-readLines("CopraRNA_option_file.txt") 
top<-as.numeric(gsub("top count:","",co[grep("top count:", co)]))



x<-read.delim(inputFile_DAVID, header=FALSE) ## edit 2.0.5.1

y<-grep("Enrichment Score", x[,2])
kk<-c()

for(i in 1:length(y)){
    k<-as.character(x[y[i],2])
    k<-sub("Enrichment Score: ","",k)
    kk<-c(kk,k)
}

kk<-as.numeric(kk)
ks<-which(kk >= 1)
if(length(ks)<1){
	ks<-1
}



ka<-seq(1,length(ks))
ks<-seq(1,length(ks)+1)
scores<-kk[ka]

if(length(y)==(length(ks)-1)){
	y<-c(y,nrow(x)+1)
}

kk<-cbind(paste(as.character(x[y[ks]+2,1]),as.character(x[y[ks]+2,2]),sep=": "),kk[ks])
grouplist<-list()

#round(as.numeric(as.character(x[y[ks]+2,10])),digits=2),

for(i in 1:(length(ks)-1)){

	group<-x[(y[ks[i]]+2):(y[ks[i+1]]-1),]
	#print(group)
	#print(length(group[,1]))
	kegg<-grep("KEGG_PATHWAY",group[,1])
	
	if(length(kegg)>0){
		kegg2<-strsplit(as.character(group[kegg,2]),":")
		kegg2<-do.call(rbind, kegg2)
		keggu<-unique(kegg2[,2])
		keggl<-c()
		for(j in 1:length(keggu)){
			keggl<-c(keggl,which(keggu[j]==kegg2[,2])[1])
			#keggl<-c(keggl,grep(keggu[j],kegg2[,2])[1])
			}
		keggl<-unique(keggl)
		
		kegg<-kegg[-keggl]
		if(length(kegg)>0){
		#print(kegg)
		group<-group[-kegg,]
					}
				}
		
		#print(length(group[,1]))
	grouplist[[i]]<-group
	}
	
	
if(kk[1,1]!= "prediction or create your enrichment manually at the DAVID homepage.: "){	

id<-c()
	for(i in 1:length(grouplist)){
		id<-c(id,as.character(grouplist[[i]][,6]))
		}
	id1<-c()
	for(i in 1:length(id)){
		id1<-c(id1,as.numeric(strsplit(id[i],",")[[1]])) ## edit 2.0.3.1 # removed space in separator
		}
	id1<-unique(id1)
	}

terms<-c()
	for(i in 1:length(grouplist)){
		terms<-c(terms,paste(round(as.numeric(as.character(grouplist[[i]][,10])),digits=2),as.character(grouplist[[i]][,2]), sep="  "))
		}
	
res<-matrix(,length(id1),sum(sapply(grouplist,nrow)))

#colnames(res)<-terms
colnames(res)<-terms
rownames(res)<-id1

res[]<-0





for(i in 1:nrow(res)){
	colcount<-0
	for(j in 1:length(grouplist)){
		h<-grep(id1[i], grouplist[[j]][,6])
		res[i,(h+colcount)]<-1
		colcount<-colcount+nrow(grouplist[[j]])
		}
	
	}
#res<-rbind(terms,res)
	
#----------------	
x<-read.csv(inputFile_CopraRNA, header=TRUE,sep=",") ## edit 2.0.5.1

l<-(strsplit(as.character(x[1:top,4]), "\\|"))

ll<-matrix(,top,8)

colnames(ll)<-c("GeneName", "IntaRNA energy","IntaRNA p-Value","start mRNA", "end mRNA", "start ncRNA", "end ncRNA", "Entrez GeneID")

for(i in 1:nrow(ll)){
 for(j in 1:8){
 ll[i,j]<-l[[i]][j]
 }

}

ll<-cbind(x[1:top,2],ll)

lll<-which(is.na(ll[,2])==TRUE)

if(length(lll) > 0){

ll<-ll[-lll,]

}


ll<-ll[1:15,]


d<-sub("GeneID:","",ll[,9])
d<-sub(")","",d)

d<-as.numeric(d)



id2<-id1






p<-c()
for(i in 1:length(id2)){

pp<-grep(id2[i],x[,4])
p<-c(p,pp[1])

}

res2<-cbind(id2,as.character(x[p,4]),x[p,2])	
	
#---------------------	

dd<-matrix(,nrow(res2),9)
colnames(dd)<-c("Locustag","GeneName", "IntaRNA energy","IntaRNA p-Value","start mRNA", "end mRNA", "start ncRNA", "end ncRNA", "Entrez GeneID")
for(i in 1:nrow(res2)){

da3<-(strsplit(as.character(res2[i,2]), "\\|"))

for(j in 1: length(da3)){
if(length(da3[[j]])==0){
da3[[j]]<-rep(0,8)
}}
da3<-as.data.frame(da3)
colnames(da3)<-seq(1:ncol(da3))
da3<-t(da3)
colnames(da3)<-c("GeneName", "IntaRNA energy","IntaRNA p-Value","start mRNA", "end mRNA", "start ncRNA", "end ncRNA", "Entrez GeneID")
geneID<-sub("GeneID:","",da3[,8])
geneID<-sub("\\)","",geneID)
da3[,8]<-geneID

Locustag<-(strsplit(as.character(da3[,1]), "\\("))
for(j in 1: length(Locustag)){
if(length(Locustag[[j]])!=2){
Locustag[[j]]<-rep("NA",2)
}}

Locustag<-as.data.frame(Locustag)
colnames(Locustag)<-seq(1:ncol(Locustag))
Locustag<-t(Locustag)
colnames(Locustag)<-c("locustag", "Genename")

da4<-c(Locustag, da3[1,2:ncol(da3)])
da4<-t(da4)

#colnames(da4)<-c("Locustag","GeneName", "IntaRNA energy","p-Value","start mRNA", "end mRNA", "start ncRNA", "end ncRNA", "Entrez GeneID")

dd[i,]<-da4[1,]

}
dd<-cbind(res2[,3],dd)
#---------------


dd[,3]<-paste(dd[,2],dd[,3], sep=" - ")

#for(i in 1:nrow(dd)){
#	if(dd[i,3]=="N/A"){
#		dd[i,3]<-dd[i,2]
#		}
#	}
	

res3<-cbind(rownames(res),dd[,1:3],dd[,3],dd[,4:9],res)

meta<-c(nrow(res3),length(grouplist))

for(i in 1:length(grouplist)){
	meta<-c(meta, rep(i, nrow(grouplist[[i]])))
	}

meta<-c(meta,rep("",9))

meta2<-(scores)

meta2<-c(meta2,rep("",length(meta)-length(scores)))


res4<-rbind(meta2,meta, colnames(res3), res3)	
res4[3,1]<-"Entrez GeneID"
res4[3,2]<-"CopraRNA_pValue"

write.table(res4, file=output_file, sep="\t",col.names = FALSE,row.names = FALSE, quote=FALSE)


#./copraRNA.pl ./enrichment_cop1.txt

#./phantomjs ./rasterize ./index_v1_1\ copy.html out.pdf
	
	
