
# Parameter:
# ooi as refseq id
# script by Jens Georg

# R --slave -f ./join_pvals_coprarna2.R --args NC_000913 

args <- commandArgs(trailingOnly = TRUE) 
ooi <- args[1] 

da<-read.csv("CopraRNA2_prep_anno_addhomologs_padj_amountsamp.csv") ## edit prw changed file name 
options <- read.table("CopraRNA_option_file.txt", sep=":") 
root<-as.numeric(as.character(options[14,2]))

en<-grep("Annotation", colnames(da))
da2<-da[,3:(en-1)]
namegenomes<-colnames(da2)
genomes<-list()

lis<-dir()

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

z.transform<-function(p,weight){
        z <- qnorm(p, lower.tail = FALSE)
        cp <- pnorm(sum(weight * z)/sqrt(sum(weight^2)), lower.tail = FALSE)
		cp
		}

z.transform.rho<-function(p,weight,rho){
        z <- qnorm(p, lower.tail = FALSE)
        cp <- pnorm(sum(weight * z)/sqrt((1-rho)*sum(weight^2)+rho*sum(weight)^2), lower.tail = FALSE)
		cp
		}
		
dat<-da4
le<-length(dat)-1
ro<-nrow(dat[[2]])
weight<-read.csv("zscore.weight", header=F, sep=";")

weight<-weight[ match(colnames(da2),toupper(weight[,1])),]


if(root==0){
	weight[,2]<-1
}

if(root>0){
	weight[,2]<-weight[,2]^(1/root)
}

out_rho<-as.matrix(da[,3:(ncol(da)-1)])
out_rho_ooi<-as.matrix(da[,3:(ncol(da)-1)])

rho<-c()
rho.names<-c()
rho_w<-c()

pv3<-c()
pv4<-c()

selfun<-function(pos, dat){
	temp<-c()
	for(i in 1:length(pos)){
		temp<-cbind(temp, as.numeric(dat[[pos[i]]][,4]))
	}
	temp
}

qtrans<-function(dat2){
	for(i in 1:ncol(dat2)){
		dat2[,i]<-qnorm(dat2[,i])
	}
	dat2
}
minorg<-max(3)

for(i in 1:ro){
        #print(i)
	p<-c()
	w<-c()
	na<-c()
	for(j in 1:le){
		
		if((dat[[j+1]][i,1])!="NA"){
			ptemp<-dat[[j+1]][i,4]
			
			w<-c(w,weight[j,2])
			p<-c(p, as.numeric(ptemp))
			na<-c(na, colnames(da)[(j+2)])
		}
	}
	if(length(p)>(minorg-1)){
	names(p)<-na
	names(w)<-na
	p<-sort(p)
	
	w<-w[match(names(w),names(p))]

	pv4temp<-c()
	pv3temp<-c()
	
	for(j in 1:(length(p)-(minorg-1))){
	
		ptemp<-p[1:(minorg-1+j)]
		wtemp<-w[1:(minorg-1+j)]
		tempname1<-names(ptemp)
		ooi_exist<-which(tempname1==ooi)
		
		tempname<-paste(sort(names(ptemp)), collapse="_")
		exist<-match(tempname, rho.names)
		if(is.na(exist)==F){
			rhotemp<-rho[exist]	
					
		}
		if(is.na(exist)==T){
			te<-match(tempname1, colnames(da))-1
			dat2<-selfun(te, dat)
			dat2[which(dat2==0)]<-NA
			dat2<-na.omit(dat2)
			dat2<-qtrans(dat2)
			
			a<-rowSums(t(t(dat2)*wtemp))
			b<-sum(wtemp)^2
			d<-sum(wtemp^2)
						
			l<-length(a)
			k<-1/l
			y<-k * seq(1,l)
			dat3<-cbind(y,a,d,b)
			dat3<-as.data.frame(dat3)
			rho1<-0
			
			rhotemp<-0
			
			rhotemp<-coef(nls(y~sort(pnorm(a/sqrt((1-rho)*d+rho*b))), data=dat3, start=list(rho=0), control=list(minFactor = 1/128, tol = 1e-05, warnOnly = TRUE)))
			rho<-c(rho, rhotemp)
			
			rho.names<-c(rho.names,tempname)
		}
		p3<-z.transform.rho(p[1:(2+j)], w[1:(2+j)],rhotemp)
		
		pv3temp<-c(pv3temp,p3)
		pv4temp<-c(pv4temp,p3)
		
		if(length(ooi_exist)==0){
			pv4temp[length(pv4temp)]<-1
			
		}
	}
	
	pv3pos<-names(p[1:(tail(which(pv3temp==min(pv3temp)),n=1)+minorg-1)])
	pv3pos<-match(pv3pos,colnames(out_rho))
	out_rho[i,pv3pos]<-paste("!",as.character(out_rho[i,pv3pos]), sep="")
	out_rho[i,pv3pos]
	
	pv4pos<-names(p[1:(tail(which(pv4temp==min(pv4temp)),n=1)+minorg-1)])
	pv4pos<-match(pv4pos,colnames(out_rho_ooi))
	out_rho_ooi[i,pv4pos]<-paste("!",as.character(out_rho_ooi[i,pv4pos]), sep="")
	out_rho_ooi[i,pv4pos]
	
	pv3temp<-min(pv3temp)
	pv3<-c(pv3,pv3temp)
	pv4temp<-min(pv4temp)
	pv4<-c(pv4,pv4temp)
	}
	if(length(p)<=(minorg-1)){
		pv3<-c(pv3,NA)
		pv4<-c(pv4,NA)
	}
}

pv3_fdr<-p.adjust(pv3, method="BH")
pv4_fdr<-p.adjust(pv4, method="BH")

pv3<-cbind(pv3_fdr,pv3)
pv4<-cbind(pv4_fdr,pv4)

colnames(pv3)<-c("fdr","p-value") ## edit prw changed spelling
colnames(pv4)<-c("fdr","p-value") ## edit prw changed spelling

out_rho<-cbind(pv3, out_rho)
out_rho_ooi<-cbind(pv4, out_rho_ooi)

out_rho<-out_rho[order(as.numeric(out_rho[,2])),]
out_rho_ooi<-out_rho_ooi[order(as.numeric(out_rho_ooi[,2])),]

write.table(out_rho, file="CopraRNA2_final_all_evo.csv",sep=",", quote=F, row.names=F) ## edit prw changed file
write.table(out_rho_ooi, file="CopraRNA2_final_all_ooi.csv",sep=",", quote=F, row.names=F) ## edit prw changed file
