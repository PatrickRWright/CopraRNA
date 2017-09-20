# Parameter:
# ooi as refseq id
# script by Jens Georg

# R --slave -f ../join_pvals_coprarna2.R --args NC_000913 

args <- commandArgs(trailingOnly = TRUE) 
ooi <- args[1] 


da<-read.csv("CopraRNA2_prep_anno_addhomologs_padj_amountsamp.csv") ## edit prw changed file name 
options <- read.table("CopraRNA_option_file.txt", sep=":") 
root<-as.numeric(as.character(options[14,2]))

all_dat<-read.csv("opt_tags.clustered", sep=";", header=T)

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
da3[[j]]<-rep(NA,8)
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

 int_table<-c()
 for(i in 1:le){
	 int_table<-cbind(int_table,as.numeric(dat[[i+1]][,4]))
	
 }
 colnames(int_table)<-colnames(da2)

out_rho<-as.matrix(da[,3:(ncol(da)-1)])
out_rho_ooi<-as.matrix(da[,3:(ncol(da)-1)])

out_rho_evo<-as.matrix(da[,3:(ncol(da)-3)])
out_rho_evo[]<-NA

out_intaRNA_pvalue<-out_rho_evo

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
        print(paste(i,"/",ro,sep=""))
	p<-c()
	w<-c()
	na<-c()
	for(j in 1:le){
		
		if((dat[[j+1]][i,1])!="NA"){
			ptemp<-dat[[j+1]][i,4]
			
			w<-c(w,weight[j,2])
			p<-c(p, as.numeric(ptemp))
			na<-c(na, colnames(da)[(j+2)])
			out_intaRNA_pvalue[i,j]<-as.numeric(ptemp)
		}
	}
	if(length(p)>(minorg-1)){
	names(p)<-na
	names(w)<-na
	p<-sort(p)
	
	w<-w[match(names(p),names(w))]

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
	
	for(jj in 1:length(pv3temp)){
		natemp<-match(names(p)[1:(2+jj)], colnames(out_rho_evo))
		if(jj>3){
			natemp<-natemp[length(natemp)]
		}
		out_rho_evo[i,natemp]<-pv3temp[jj]
	}
	
	
	for(jj in 1:ncol(out_rho_evo)){
		natemp<-match(colnames(out_rho_evo)[jj], names(p))
		if(is.na(natemp)==F){
		if(natemp<=3){
			p_temp<-min(pv3temp)
		}
		if(natemp>3){
			p_temp<-min(pv3temp[(natemp-2):length(pv3temp)])
		}
		
		out_rho_evo[i,jj]<-p_temp
		
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

initial_sorting<-seq(1,nrow(pv3))

colnames(pv3)<-c("fdr","p-value") ## edit prw changed spelling
colnames(pv4)<-c("fdr","p-value") ## edit prw changed spelling

out_rho<-cbind(pv3, out_rho,initial_sorting)
out_rho_ooi<-cbind(pv4, out_rho_ooi,initial_sorting)

out_rho<-out_rho[order(as.numeric(out_rho[,2])),]
out_rho_ooi<-out_rho_ooi[order(as.numeric(out_rho_ooi[,2])),]

write.table(out_rho, file="CopraRNA2_final_all_balanced.csv",sep=",", quote=F, row.names=F) ## edit prw changed file
write.table(out_rho_ooi, file="CopraRNA2_final_all_ooi.csv",sep=",", quote=F, row.names=F) ## edit prw changed file

out_rho_evo<-as.matrix(out_rho_evo)
out_evo_rank<-out_rho_evo
top<-c()
rank_thres<-200

anno<-da

 for(i in 1:ncol(out_evo_rank)){
	 out_evo_rank[,i]<-rank(as.numeric(out_rho_evo[,i]), na.last="keep" ,ties.method="first")
	 temp<-which(as.numeric(out_evo_rank[,i])<=rank_thres)
	 top<-c(top, temp)
 }

 top<-unique(top)

initial_sorting<-initial_sorting[top]
 evo_ooi_pvalue<-out_rho_evo[top,]
 evo_ooi_rank<-out_evo_rank[top,]
 evo_int_pvalue<-out_intaRNA_pvalue[top,]
 evo_anno<-anno[top,3:ncol(anno)]


#### rank sum ordering , NA = 200

rank_list<-list()
for(i in 1:ncol(evo_ooi_rank)){
	temp<-evo_ooi_rank[,i]
	names(temp)<-seq(1,nrow(evo_ooi_rank))
	temp<-names(sort(temp,na.last=NA))
	print(length(temp))
	rank_list[[i]]<-temp

}

require(RobustRankAggreg)
rank_list2<-aggregateRanks(rank_list,method = "RRA")
rank_list3<-as.numeric(as.character(rank_list2[,"Name"]))

initial_sorting<-initial_sorting[rank_list3]
evo_ooi_pvalue<-evo_ooi_pvalue[rank_list3,]
evo_ooi_rank<-evo_ooi_rank[rank_list3,]
evo_int_pvalue<-evo_int_pvalue[rank_list3,]
evo_anno<-evo_anno[rank_list3,]

evo_ooi_pvalue_scored<-evo_ooi_pvalue
evo_ooi_rank_scored<-evo_ooi_rank
evo_int_pvalue_scored<-evo_int_pvalue

evo_ooi_pvalue_scored[]<-0
evo_ooi_rank_scored[]<-0
evo_int_pvalue_scored[]<-0

evo_rank_thres<-250
int_p_thres<-0.3
ooi_p_tres<-0.001

for(i in 1:ncol(evo_ooi_pvalue)){
	temp<-which(as.numeric(evo_ooi_pvalue[,i])<=ooi_p_tres)
	evo_ooi_pvalue_scored[temp,i]<-1
	
	temp<-which(as.numeric(evo_ooi_rank[,i])<=evo_rank_thres)
	evo_ooi_rank_scored[temp,i]<-1
	
	temp<-which(as.numeric(evo_int_pvalue[,i])<=int_p_thres)
	evo_int_pvalue_scored[temp,i]<-1
}

evo_ooi_p_result<-c()
evo_ooi_rank_result<-c()
evo_int_p_result<-c()

for(i in 1:nrow(evo_ooi_pvalue)){
	#print(i)
	temp<-sum(as.numeric(evo_ooi_pvalue_scored[i,]))
	temp<-temp/ncol(evo_ooi_pvalue_scored)
	evo_ooi_p_result<-c(evo_ooi_p_result, temp)
	
	temp<-sum(as.numeric(evo_ooi_rank_scored[i,]))
	temp<-temp/ncol(evo_ooi_rank_scored)
	evo_ooi_rank_result<-c(evo_ooi_rank_result, temp)
	
	temp<-sum(as.numeric(evo_int_pvalue_scored[i,]))
	temp<-temp/ncol(evo_int_pvalue_scored)
	evo_int_p_result<-c(evo_int_p_result, temp)	
}
rho_out<-list(rho.names,rho)

save(rho_out, file="rho_out.Rdata")

Rank_p_value<-rank_list2
Rank_p_value<-Rank_p_value[,2]
evo_analysis<-cbind(Rank_p_value,evo_ooi_rank_result,evo_ooi_p_result,evo_int_p_result,evo_anno,evo_ooi_rank,evo_ooi_pvalue,evo_int_pvalue,initial_sorting)
write.table(evo_analysis, file="CopraRNA2_final_all_evo.csv", sep="\t", row.names=F)

evo_int_pvalue<-matrix(as.numeric(evo_int_pvalue),nrow(evo_int_pvalue),ncol(evo_int_pvalue))
evo_int_pvalue_scored<-matrix(as.numeric(evo_int_pvalue_scored),nrow(evo_int_pvalue_scored),ncol(evo_int_pvalue_scored))

napos<-which(is.na(evo_int_pvalue))
evo_int_pvalue_scored[napos]<--1

colnames(evo_int_pvalue_scored)<-colnames(da2)
write.table(evo_int_pvalue_scored, file="IntaRNA_heatmap_table", sep="\t", row.names=F)


