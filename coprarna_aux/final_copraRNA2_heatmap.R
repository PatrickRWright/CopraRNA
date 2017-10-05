
#R --slave -f  final_copraRNA2_heatmap.r --args inputfile=CopraRNA2_final_all_ooi.csv num=25 consensus=overall select=FALSE clustering=sRNA sel=genelist.txt coprarna_reference_file=copra_refseq_positivelist.txt int_p_thres=0.3 prefix=sRNA



#Arguments
# clustering: ordering of columns. "default" by kmeans clustering based on the heatmap table, "ribosomal" by 16S rRNA tree and "sRNA" by sRNA tree. 
# consensus: "overall" or "ooi"
# select: If TRUE the genes supplied in the file specified by the parameter "sel" (default: "genelist.txt", one gene per row, identifier locustag or genename, case sensitive) are visualized. If "FALSE" the top "num" (default: num=25) genes from the CopraRNA_2.0 result specified by the "inputfile" parameter are visualized (default: "CopraRNA2_final_all_ooi.csv")
# coprarna_reference_file: filename of the .txt file specifying the available organisms for CopraRNA prediction
# int_p_thres: IntaRNA p-value threshold for a likely non functional interaction which is colored orange in the heatmap (default: int_p_thres=0.3)
# prefix: Prefix for the name of the output file. <prefix>_conservation_heatmap.pdf

args <- commandArgs(trailingOnly = TRUE)

sel="genelist.txt"
inputfile="CopraRNA2_final_all_ooi.csv"
num=25
consensus="overall"
select=FALSE
clustering="sRNA"
coprarna_reference_file="copra_refseq_positivelist.txt"
int_p_thres=0.3
prefix="sRNA"

 for(i in 1:length(args)){
	temp<-strsplit(args[i],"=")
	temp<-temp[[1]]
	temp1<-temp[1]
	temp2<-temp[2]
	assign(as.character(temp1),temp2)
 }


num<-as.numeric(num)
select<-as.logical(select)
int_p_thres<-as.numeric(int_p_thres)




selected_heatmap<-function(int_p_thres=0.3, sel="genelist.txt", select=T, clustering="ribosomal", consensus="overall", inputfile="CopraRNA2_final_all_ooi.csv", num=25,coprarna_reference_file="copra_refseq_positivelist.txt", prefix="sRNA"){
clus=TRUE
evo_analysis<-read.csv(inputfile,sep=",", header=T) 
selection<-evo_analysis[1:num,"initial_sorting"]
evo_analysis<-evo_analysis[order(evo_analysis[,"initial_sorting"]),]

load("conservation_table.Rdata")
con_table<-conservation_table[[1]]
con_table_sub<-conservation_table[[2]]
if(consensus=="ooi"){
con_table<-conservation_table[[3]]
con_table_sub<-conservation_table[[4]]

}

find_sRNA<-function(x, y){
	e<-grep("Annotation", colnames(y))-1
	temp<-evo_analysis[,3:e]
	genes<-apply(temp, 1, paste, collapse="")
	x<-as.matrix(x)
	out<-unlist(apply((x), 1, find1, genes))
	out
}



find1<-function(x, genes){
	temp<-grep(x, genes, ignore.case=F)
	temp
}

if(select==T){
	genelist<-read.csv(sel, sep="\t", header=F)[,1]
	selection<-unique(find_sRNA(genelist, evo_analysis))
}
evo_analysis<-evo_analysis[selection,]
opt_sub_comp<-function(con_table,con_table_sub){ #for each row in con_table
	opt<-con_table
	subopt<-con_table_sub
	for(i in 1:length(opt)){
		if(is.na(opt[i])==F){
		if(opt[i] == 0 | opt[i] == 2){
			if(is.na(subopt[i])==F){
			if(subopt[i] == 3){
				opt[i]<-5
			}
			if(subopt[i] == 1){
				opt[i]<-4
			}
		}
		}
	}
	}
	opt
	
}
pfun<-function(x){
out<-x[1]
out
}

mRNA_and_sRNA_test<-function(x){
	out<-NA
	if(is.na(x[1])==F){
		if(length(x)==4){
			if(x[2]=="TRUE"){
				if(x[3]=="TRUE" | x[4]=="TRUE"){
					out<-3
				}
				if(x[3]=="FALSE" & x[4]=="FALSE"){
					out<-1
				}
			}
			if(x[2]=="FALSE"){
				
					if(x[3]=="TRUE" | x[4]=="TRUE"){
						out<-2
					}
					if(x[3]=="FALSE" & x[4]=="FALSE"){
						out<-0
					}
				
			}
		}
	}
	out
}

heat_table<-c()
evo_int_pvalue<-c()
con_table<-con_table[selection,]
con_table_sub<-con_table_sub[selection,]
for(i in 1:nrow(con_table)){
	print(i)
	opt<-strsplit(con_table[i,], "\\|")
	subopt<-strsplit(con_table_sub[i,], "\\|")
	p_opt<-as.numeric(unlist(lapply(opt,pfun)))
	p_subopt<-as.numeric(unlist(lapply(subopt,pfun)))
	opt<-unlist(lapply(opt,mRNA_and_sRNA_test))
	subopt<-unlist(lapply(subopt,mRNA_and_sRNA_test))
	opt_subopt<-cbind(opt,subopt)
	
	opt<-opt_sub_comp(opt,subopt)
	heat_table<-rbind(heat_table,opt)
	
	for(j in 1:length(p_opt)){
	#print(j)
		if(is.na(opt[j])==F){
		if(opt[j] == 4 | opt[j] ==5 ){
			p_opt[j]<-p_subopt[j]
		}
		}
	}
	
	evo_int_pvalue<-rbind(evo_int_pvalue,p_opt)
}





ee<-grep("Annotation", colnames(evo_analysis))-1
#e<-grep("Annotation", colnames(evo_analysis))-3
#s<-grep("Annotation", colnames(evo_analysis))+2
#e<-s+e-1
#evo_int_pvalue1<-as.matrix(evo_analysis[,s:e])
colnames(evo_int_pvalue)<-colnames(evo_analysis)[3:ee]
copref<-read.delim(coprarna_reference_file, sep="\t", header=T,comment.char = "#")	
nam<-c()
for(i in 1:ncol(evo_int_pvalue)){
	tnam<-grep(gsub("\\..*","",colnames(evo_int_pvalue)[i]),copref[,1])
	nam<-c(nam,as.character(copref[tnam,2]))
	
	
}


#evo_int_pvalue<-evo_int_pvalue[selection,]




nam2<-c()
for(i in 1:length(nam)){
	temp<-substr(nam[i],1,3)
	temp2<-strsplit(nam[i],"_")[[1]]
	temp<-paste(temp,"_",temp2[2], sep="")
	if(length(temp2)>2){
		temp<-paste(temp, temp2[length(temp2)], sep="_")
	}
	nam2<-c(nam2,temp)
}
nam2<-paste(nam2, colnames(evo_int_pvalue), sep="_")
#colnames(evo_int_pvalue)<-nam2

evo_int_pvalue_scored<-evo_int_pvalue
evo_int_pvalue_scored[]<-0


for(i in 1:ncol(evo_int_pvalue)){
	
	temp<-which(evo_int_pvalue[,i]<=0.001)
	if(length(temp)>0){
		evo_int_pvalue_scored[temp,i]<-1
	}
	temp<-which(evo_int_pvalue[,i]>0.001 & evo_int_pvalue[,i]<=0.01)
	if(length(temp)>0){
		evo_int_pvalue_scored[temp,i]<-2
	}
	temp<-which(evo_int_pvalue[,i]>0.01 & evo_int_pvalue[,i]<=0.1)
	if(length(temp)>0){
		evo_int_pvalue_scored[temp,i]<-3
	}
	temp<-which(evo_int_pvalue[,i]>0.1 & evo_int_pvalue[,i]<=0.2)
	if(length(temp)>0){
		evo_int_pvalue_scored[temp,i]<-4
	}
	temp<-which(evo_int_pvalue[,i]>0.2 & evo_int_pvalue[,i]<=int_p_thres)
	if(length(temp)>0){
		evo_int_pvalue_scored[temp,i]<-5
	}
	temp<-which(evo_int_pvalue[,i]>int_p_thres)
	if(length(temp)>0){
		evo_int_pvalue_scored[temp,i]<-6
	}
	
}


napos<-which(is.na(evo_int_pvalue))
evo_int_pvalue_scored[napos]<-NA


gene_anno<-c()
evo_analysis<-as.matrix(evo_analysis)
for(i in 1:nrow(evo_analysis)){

	e<-grep("Annotation", colnames(evo_analysis))-1
	temp<-evo_analysis[i,3:e]
	temp<-temp[which(temp!="")][1]
	locus_tag<-gsub("\\(.*","",as.character(temp))
	genename<-gsub(".*\\(","",as.character(temp))
	genename<-gsub("\\|.*", "", genename)
	genename<-gsub("\\/", "", genename)
	na<-paste(genename, locus_tag,sep="_")
	gene_anno<-c(gene_anno, na)
}



l<-2
len<-length(which(duplicated(gene_anno)))
while(len>0){
	#print(l)
	gene_anno[which((duplicated(gene_anno)))]<-paste(gene_anno[which((duplicated(gene_anno)))],l,sep="_")
	l<-l+1
	len<-length(which(duplicated(gene_anno)))
}
row.names(evo_int_pvalue_scored)<-gene_anno
row.names(evo_int_pvalue)<-gene_anno
row.names(heat_table)<-gene_anno
#evo_int_pvalue_log<--log10(evo_int_pvalue)

if(clustering=="sRNA"){
command<-paste("clustalo -i ", "input_sRNA.fa", " --distmat-out=distmatout1.txt --full --use-kimura --output-order=input-order --force --max-hmm-iterations=-1", sep="")
system(command)
temp<-read.delim("distmatout1.txt",sep="",header=F, , skip=1)
unlink("distmatout1.txt")
na<-temp[,1]
temp<-temp[,2:ncol(temp)]
colnames(temp)<-na
rownames(temp)<-na
dis<-as.dist(temp) 
clus<-(hclust(dis,method="average"))

ord<-clus$label

ord2<-match(ord, gsub("\\..*","",colnames(evo_int_pvalue_scored)))

evo_int_pvalue_scored<-evo_int_pvalue_scored[,ord2]
heat_table<-heat_table[,ord2]
#nam3<-nam2[ord2]
#clus2<-as.dendrogram(clus)

clus$label<-nam2
}
if(clustering=="ribosomal"){
d<-read.delim("distmat.out", sep="\t", skip=6, header=FALSE)
 dd<-d[2:nrow(d),2:(ncol(d)-2)]
 rownames(dd)<-gsub(" .*","",d[2:nrow(d),ncol(d)])
colnames(dd)<-gsub(" .*","",d[2:nrow(d),ncol(d)])
 dis<-as.dist(t(dd))
 clus<-(hclust(dis,method="average"))

ord<-clus$label

ord2<-match(ord, gsub("\\..*","",colnames(evo_int_pvalue_scored)))

evo_int_pvalue_scored<-evo_int_pvalue_scored[,ord2]
heat_table<-heat_table[,ord2]
#nam3<-nam2[ord2]
#clus2<-as.dendrogram(clus)

clus$label<-nam2
 

}
#require(gplots)
require(pheatmap)
my_palette<-colorRampPalette(c("darkolivegreen1","olivedrab1","olivedrab2","olivedrab3","olivedrab4","orangered1"))(n=6)

nam<-paste(prefix,"conservation_heatmap.pdf", sep="_" )

pdf(nam,paper = "a4r", width = 0, height = 0) 
 heat_table2<-ifelse(heat_table == 3 , "+/+" ,ifelse(heat_table == 2 , "-/+" , ifelse(heat_table == 1 , "+/-" , ifelse(heat_table == 0 , "-/-" , ifelse(heat_table == 4 ,"2+/-" ,ifelse(heat_table == 5 ,"2+/+" ,""))))))
 nan<-is.na(heat_table2)
 heat_table2[nan]<-""
 
 
 
thre1<-paste("<", int_p_thres, sep="")
thre2<-paste(">", int_p_thres, sep="")

int_table2<-evo_int_pvalue_scored
colnames(int_table2)<-nam2
pheatmap(int_table2,legend_breaks=c(1,2,3,4,5,6),main=paste(consensus, "_consensus_based"), legend_labels=c("<0.001","<0.01","<0.1","<0.2",thre1,thre2), scale="none",display_numbers = heat_table2, cluster_rows=F,cluster_cols=clus, fontsize_number = 6 ,col=my_palette, trace="none", colsep=seq(1,ncol(evo_int_pvalue_scored)), rowsep=seq(1,nrow(int_table2)),cexCol=0.8, cexRow=0.9,na.rm=TRUE,na.color=1,srtCol=45,offsetRow = 0,offsetCol = 0,margins=c(12,8),sepcolor="white")

dev.off()


}
selected_heatmap(select=select, clustering=clustering, consensus=consensus, inputfile=inputfile, sel=sel, num=num, coprarna_reference_file=coprarna_reference_file,int_p_thres=int_p_thres, prefix=prefix)


