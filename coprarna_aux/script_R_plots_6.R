# call
# R --slave -f ../script_R_plots_6.R --args CopraRNA1_final.csv 
# script by Jens Georg

args <- commandArgs(trailingOnly = TRUE) ## edit 2.0.2
inputFile <- args[1] ## edit 2.0.2

localMaxima <- function(x) {
  # Use -Inf instead if x is numeric (non-integer)
  y <- diff(c(-.Machine$integer.max, x)) > 0L
  rle(y)$lengths
  y <- cumsum(rle(y)$lengths)
  y <- y[seq.int(1L, length(y), 2L)]
  if (x[[1]] == x[[2]]) {
    y <- y[-1]
  }
  y
}



map<-function(y){
se<-y
dd<-vector("list", length(se))
d<-vector("list", length(se))

d[]<-0
dd[]<-0


x<-length(se)

for(i in 1:x){

A<-grep("a", se[[i]], ignore.case=TRUE)
G<-grep("g", se[[i]], ignore.case=TRUE)
T<-grep("t", se[[i]], ignore.case=TRUE)
U<-grep("u", se[[i]], ignore.case=TRUE)
C<-grep("c", se[[i]], ignore.case=TRUE)
gap<-grep("-", se[[i]], ignore.case=TRUE)
all<-c(A,G,C)
if(length(T) > 0){
all<-c(all, T)}

if(length(U) > 0){
all<-c(all, U)}

all<-sort(all)

if(length(all)<1){
	all<-0
	}
if(length(gap)<1){
	gap<-0
	}	
	
d[[i]]<-all
dd[[i]]<-gap
}
d<-list(d,dd)
d
}


#clustalW aufruf

clustalW<-function (RNAseqs,  input.format = "fasta",      
    clustal.path = NULL, forceDNAtolower = TRUE, forceAAtolower = FALSE) 
{
    seq.nucl <- RNAseqs
    nseqs <- length(seq.nucl)
  
  
        tmp <- tempfile(pattern = "clustal")
        RNAseq.file <- tempfile(pattern = "RNA")
		
        write.fasta(seq.nucl , names = names(seq.nucl), file.out = RNAseq.file)
		kk<-paste(clustal.path, " -outfile=", tmp, " -infile=", 
            RNAseq.file," -pwdnamatrix=IUB -dnamatrix=IUB -pwgapopen=10 -pwgapext=0.1 -gapopen=10 -gapext=0.2 -gapdist=5  -outorder=input -output=FASTA > /dev/null", sep = "") ## edit 2.0.5.1 // moving command line output to /dev/null
			
        system(kk)
       
		read.fasta(tmp)
	
    }
	
require(seqinr)

#daten preprozessierung

pre<-function(up=200,down=100){
da<-read.csv(inputFile) ## edit 2.0.2

danc<-read.fasta("ncrna.fa")
d<-clustalW(danc,clustal.path="clustalw")

d<-map(d)

en<-grep("Annotation", colnames(da))
da2<-da[,3:en-1]
namegenomes<-colnames(da2)
genomes<-list()

lis<-dir()

for(i in 1: length(namegenomes)){

#gen<-grep(namegenomes[i], lis)
#gen<-grep("NC_000913", lis)
#gen2<-grep(paste(".fa","i",sep="="), paste(lis[gen],"i",sep="="))
#lis<-lis[gen[gen2]]
#gen<-grep("interacting",lis)
#lis<-lis[-gen]
#gen<-grep("ncRNA",lis)
#lis<-lis[-gen]

lis<-paste(namegenomes[i],"_upfromstartpos_",up,"_down_",down,".fa",sep="")
#genomes[[i]]<-read.fasta(paste(namegenomes[i],"_upfromstartpos_200_down_50.fa",sep=""))
genomes[[i]]<-read.fasta(lis)
}
names(genomes)<-namegenomes

da4<-list()
da4[[1]]<-da[,1]
for(i in 1:ncol(da2)){

da3<-(strsplit(as.character(da2[,i]), "\\|"))

for(j in 1: length(da3)){
if(length(da3[[j]])==0){
da3[[j]]<-rep(0,8)
}}
da3<-as.data.frame(da3)
colnames(da3)<-seq(1:ncol(da3))
da3<-t(da3)
colnames(da3)<-c("GeneName", "IntaRNA energy","p-Value","start mRNA", "end mRNA", "start ncRNA", "end ncRNA", "Entrez GeneID")
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

da3<-cbind(Locustag, da3[,2:ncol(da3)])

da4[[i+1]]<-da3
}
da4<-list(da4, genomes, d)
da4
}





mRNAalign<-function(x, nu, le) {

pred_results<-x[[1]]
sequences<-x[[2]]


# bn<- matrix(,Anzahl von geplotteten predictions(num), Anzahl von Organismen (le)) Inhalt: Locustag´s 
#   [,1]      [,2]     [,3]        
#1  "mm_2442" "NA"     "mbar_a3378"
#2  "mm_1282" "NA"     "mbar_a1640"
#3  "mm_0491" "ma3600" "mbar_a1818"
#4  "mm_1574" "ma0314" "NA"    

bn<-matrix(,nu,1)
	for(j in 1: le){
	bn<-cbind(bn,pred_results[[j+1]][1:nu,1])
	}
bn<-bn[,2:ncol(bn)]

for(i in 1:ncol(bn)){
	#bn[which(bn[,i]=="NA"),i]<-paste(bn[which(bn[,i]=="NA"),i],bn[1,i],sep="_")
	bn[which(bn[,i]=="NA"),i]<-paste(bn[which(bn[,i]=="NA"),i],i,sep="_")
	}
homolog_table<-find_homologs(bn, sequences)

pos_list<-align_homologs(homolog_table, sequences, bn)
seqs<-pos_list[[1]]
gaps<-pos_list[[2]]
interactions<-interaction_pos(pred_results, seqs, nu, le)

res<-list(seqs,gaps,interactions)
res
}




find_homologs<-function(bn, sequences){

homolog_table<-bn

for(i in 1:nrow(bn)){
	pos_vect<-c()
	for(j in 1:ncol(bn)){
		aaa<-paste("i",bn[i,j],sep="-")
		bbb<-paste("i",names(sequences[[j]]),sep="-")
		pp<-grep(aaa,bbb, ignore.case=TRUE)
		
		if(length(pp)>0){
			pos_vect<-c(pos_vect,pp[1])
			}
		else{
			pos_vect<-c(pos_vect,NA)
			}
		
		}
	homolog_table[i,]<-pos_vect	
	}
homolog_table
}

align_homologs<-function(homolog_table, sequences, bn){


seqs<-vector("list", nrow(homolog_table))
gaps<-vector("list", nrow(homolog_table))


for(i in 1:nrow(homolog_table)){
    seq_list<-vector("list", ncol(bn))
		
	for(j in 1:ncol(homolog_table)){
		if(is.na(as.numeric(homolog_table[i,j]))==FALSE){
			seq_list[[j]]<-sequences[[j]][[as.numeric(homolog_table[i,j])]]
			}
		else{
			seq_list[[j]]<-NA
			}
		}
	
	names(seq_list)<-bn[i,]
	map_homologs<-clustalW(seq_list,clustal.path="clustalw")
	map_homologs<-empty_plots(map_homologs)
	map_homologs<-map(map_homologs)
	seqs[[i]]<-map_homologs[[1]]
	gaps[[i]]<-map_homologs[[2]]
	}
res<-list(seqs,gaps)
res

}

interaction_pos<-function(pred_results, seqs, nu, le, sRNA=FALSE){

a<-5
b<-6

if(sRNA==TRUE){
	a<-7
	b<-8
	}

interactions<-vector("list", nu)

for(i in 1:length(interactions)){
	interactions[[i]]<-vector("list", le)
	}


for(i in 1:nu){
	for(j in 1:le){
		st<-as.numeric(pred_results[[j+1]][i,a])
		en<-as.numeric(pred_results[[j+1]][i,b])
		if(st != 0 & en != 0){
			st<-seqs[[i]][[j]][st]
			en<-seqs[[i]][[j]][en]
			
			
			}
		
		interactions[[i]][[j]]<-seq(st,en)
		
		}
	}
interactions
}

#interactions<-interaction_pos(pred_results, seqs, nu, le)

#homolog_table<-find_homologs(bn, sequences)




empty_plots<-function(map_homologs) {

ls<-grep("NA", names(map_homologs))
if(length(ls)>0){
	for(i in 1:length(ls)){
		map_homologs[[ls[i]]][]<-"-"
		}
	}
map_homologs
}

plot_parameter<-function(position_data, center=200, rescale=TRUE){
seqs<-position_data[[1]]
gaps<-position_data[[2]]
interactions<-position_data[[3]]

target_table<-matrix(,length(seqs)*length(seqs[[1]]),3)
interaction_table<-matrix(,length(seqs)*length(seqs[[1]]),3)
gap_list<-vector("list", length(seqs)*length(seqs[[1]]))


for(i in 1:length(seqs)){
	
	for(j in 1:length(seqs[[1]])){
		su<-0
		if(rescale==TRUE){
			su<-seqs[[i]][[j]][center]
			su<-(center*1.2)-su
			}
		
		target_table[j+(i-1)*length(seqs[[1]]),1]<-min(seqs[[i]][[j]])+su
		target_table[j+(i-1)*length(seqs[[1]]),2]<-max(seqs[[i]][[j]])+su
		target_table[j+(i-1)*length(seqs[[1]]),3]<-j+(i-1)*length(seqs[[1]])+i-1
		interaction_table[j+(i-1)*length(interactions[[1]]),1]<-min(interactions[[i]][[j]])+su
		interaction_table[j+(i-1)*length(interactions[[1]]),2]<-max(interactions[[i]][[j]])+su
		interaction_table[j+(i-1)*length(interactions[[1]]),3]<-j+(i-1)*length(seqs[[1]])+i-1
		gap_list[[j+(i-1)*length(seqs[[1]])]]<-list(gaps[[i]][[j]]+su,j+(i-1)*length(seqs[[1]])+i-1)
		
		
		}
	
	}
	
target_table<-list(target_table,c(length(seqs),length(seqs[[1]])))
interaction_table<-list(interaction_table,c(length(seqs),length(seqs[[1]])))
res<-list(target_table, interaction_table, gap_list)
res
}

#res<-plot_parameter(position_data)

dens_pars<-function(plot_pars){
interaction_table<-plot_pars[[2]][[1]]
gap_list<-plot_pars[[3]]
dens<-c()
for(i in 1:nrow(interaction_table)){
	if(is.na(interaction_table[i,1])==FALSE){
		dens1<-seq(interaction_table[i,1],interaction_table[i,2])
		gap1<-gap_list[[i]][[1]]
		
		overlap<-which(is.element(dens1,gap1))
		
		if(length(overlap>1)){
			dens1<-dens1[-overlap]
			}
		dens<-c(dens,dens1)
		}
	}
dens
}




plot_function<-function(x,num=num, density_pars, plot_pars,center=center, UTR5=FALSE, UTR3=FALSE, fulllength=FALSE,le, sRNA=FALSE){
target_table<-plot_pars[[1]]
interaction_table<-plot_pars[[2]]
gap_list<-plot_pars[[3]]
le<-target_table[[2]][2]

mi<-min(target_table[[1]][,1],na.rm = TRUE)
ma<-max(target_table[[1]][,2],na.rm = TRUE)

collist<-c("steelblue2","orangered1","olivedrab3","darkgoldenrod","darkorchid3","dodgerblue3","chocolate3","chartreuse3")

if(sRNA==FALSE){
	plot(0,0, xlim=c(mi,ma*1.25), ylim=c(0,((num*le+num-1)*2)*1.01), type="n", ylab="", main="Predicted interaction regions mRNAs", xlab="",yaxt="n",bty="n",xaxt="n")
	}
else{
	plot(0,0, xlim=c(mi,ma*1.25), ylim=c(0,((num*le+num-1)*2)*1.01), type="n", ylab="", main="Predicted interaction regions sRNA", xlab="",yaxt="n",bty="n",xaxt="n")
	}

#Interaction plots

count<-0
count2<-1

for(i in 1:(num*le)){
	
	if(count==le){
		count<-0
		count2<-count2+1
		}
	count<-count+1
	xx<-count2/8
	while(xx > 1) xx<-xx-1
	colnumb<-xx*8
	
	lines(c(target_table[[1]][i,1],target_table[[1]][i,2]),c(target_table[[1]][i,3],target_table[[1]][i,3]),col="gray85", lwd=1.5)
	lines(c(interaction_table[[1]][i,1],interaction_table[[1]][i,2]),c(interaction_table[[1]][i,3],interaction_table[[1]][i,3]), col=collist[colnumb], lwd=1.5)
	if(length(gap_list[[i]][[1]])>0){
		for(j in 1:length(gap_list[[i]][[1]])){
			points(c(gap_list[[i]][[1]][j],gap_list[[i]][[1]][j]+1),c(gap_list[[i]][[2]],gap_list[[i]][[2]]), type="l", col="white",pch=20, cex=1, lwd=1.5)
			}
		}
	}
	
#lines between targets and names of targets
for(i in 1:num){	
	k<-i*le+i
	lines(c(mi,ma),c(k,k))
	namv<-c()
	for (j in 1:le){
		namv<-c(namv, as.numeric(x[[1]][[j+1]][i,5]))
		}
	p<-which(namv != 0)	
	
	nam<-paste(x[[1]][[p[1]+1]][i,1],x[[1]][[p[1]+1]][i,2], sep="_")
	
	text(ma*1.15, k-le/2, nam ,cex=0.8 )
	}	
	
	
#Density plot

hhh<-density(density_pars)
maxi<-localMaxima(hhh$y)
mmm<-max(hhh$y)
mmm2<-max(hhh$y)
mmm<-(num*le+num-1)/mmm

hhh$y<-hhh$y*mmm+(num*le+num-1)*1.01+(num*le+num-1)*0.04
points(hhh, type="l", col="black", lwd=2)

#Ablines for local maxima in density plot
for(ii in 1:length(maxi)){
	lines(c(hhh$x[maxi[ii]],hhh$x[maxi[ii]]),c((num*le+num-1)*1.01+(num*le+num-1)*0.04,hhh$y[maxi[ii]]),lwd=2)
	}
#parting line between density plot and interaction plots
abline(h=(num*le+num-1)*1.01+(num*le+num-1)*0.04)

#Axis for density plot

if(sRNA==FALSE){
	ylen<-seq(0,mmm2+0.02, by=0.001)
	ylen1<-ylen
	ylen<-ylen*mmm+(num*le+num-1)*1.01+(num*le+num-1)*0.04
	axis(2,at=ylen,labels=ylen1)
	}

else{
	ylen<-seq(0,mmm2+0.02, by=0.01)
	ylen1<-ylen
	ylen<-ylen*mmm+(num*le+num-1)*1.01+(num*le+num-1)*0.04
	axis(2,at=ylen,labels=ylen1)
	}
	
#Annotations for interaction plot

if(UTR5==TRUE){
lines(c(center*1.2,center*1.2),c(-100,num+num*le), col="gray50")

se3<-c(center*1.2-200,center*1.2-100,center*1.2-50,center*1.2,center*1.2+50,center*1.2+100,center*1.2+200)


axis(1,at=se3,labels=c(-200,-100,-50,"start codon",50,100,200))

}




if(UTR3==TRUE){
lines(c(center*1.2,center*1.2),c(-100,num+num*le), col="gray50")
se3<-c(center*1.2-100,center*1.2,center*1.2+100)


axis(1,at=se3,labels=c(100,"stop codon",100))
}

if(fulllength==TRUE){


	se<-seq(mi,ma,by=100)
	selabels<-seq(0,100*(length(se)-1),by=100)
	axis(1,at=se,labels=selabels)
	title(xlab="nt")
	}

if(sRNA==TRUE){
	xval<-ma
	leng<-seq(0,xval,by=50)
	leng[1]<-1
		if((xval-max(leng))>20){
		leng<-c(leng, xval)}
		else{
		leng[length(leng)]<-xval
		}
	axis(1,at=leng,labels=leng)
	title(xlab="nt")
	}
}

#plot_function(x,5, density_pars, plot_pars,center, UTR5=TRUE, UTR3=FALSE, fulllength=FALSE)
sRNA_align<-function(x, nu=nu, le=le){

pred_results<-x[[1]]
sRNA<-x[[3]]
seqs<-vector("list", nu)
for(i in 1:nu){
	seqs[[i]]<-sRNA[[1]]
	}
gaps<-vector("list", nu)
for(i in 1:nu){
	gaps[[i]]<-sRNA[[2]]
	}	

interactions<-interaction_pos(pred_results, seqs, nu, le, sRNA=TRUE)

res<-list(seqs,gaps,interactions)
res
}


CopraRNA_plot<-function(thres=0.01, num=15,name="ncRNA" ){
require(seqinr)
options <- read.table("CopraRNA_option_file.txt", sep=":") ## edit 2.0.4 // removed coprarna_pars.txt
up <- as.numeric(as.character(options$V2[2])) ## edit 2.0.4 // removed coprarna_pars.txt
down <- as.numeric(as.character(options$V2[3])) ## edit 2.0.4 // removed coprarna_pars.txt
type <- as.character(options$V2[4]) ## edit 2.0.4 // removed coprarna_pars.txt
if(type=="CDS"){
	fulllength<-TRUE
	}

	if(type=="5utr"){
	UTR5<-TRUE
	center<-up
	}

if(type=="3utr"){
	center<-up
	UTR3<-TRUE
	}


x<-pre(up=up,down=down)

nu<-length(which(x[[1]][[1]]<=thres))
nu2<-nu
if(num>nu){
	#nu2<-num
	num<-nu
	}

le<-length(x[[1]])-1
position_data<-mRNAalign(x, nu2, le)
plot_pars<-plot_parameter(position_data, center=center, rescale=TRUE)
density_pars<-dens_pars(plot_pars)

na<-"mRNA_regions.pdf"
pdf(file=na, paper="a4r", width=0, height=0)
plot_function(x,num, density_pars, plot_pars,center,UTR5=TRUE)
dev.off()

postscript(file="mRNA_regions.ps")
plot_function(x,num, density_pars, plot_pars,center,UTR5=TRUE)
dev.off()

sRNA_position_data<-sRNA_align(x, nu=nu2, le=le)
sRNA_plot_pars<-plot_parameter(sRNA_position_data, center=center, rescale=FALSE)
sRNA_density_pars<-dens_pars(sRNA_plot_pars)

na<-"sRNA_regions.pdf"
pdf(file=na, paper="a4r", width=0, height=0)
plot_function(x,num, sRNA_density_pars, sRNA_plot_pars,center, sRNA=TRUE)
dev.off()

postscript(file="sRNA_regions.ps")
plot_function(x,num, sRNA_density_pars, sRNA_plot_pars,center, sRNA=TRUE)
dev.off()

}


CopraRNA_plot(thres=0.01, num=20,name="test" )
