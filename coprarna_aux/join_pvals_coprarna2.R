# Parameter:
# ooi as refseq id
# ooi_consensus : calculates p_Value based on the predicted interaction regions in the ooi. Organisms with non matching predictions are excluded from the calculation.
# overall_consensus : calculates p_Value based on the consensus interaction regions based on all organisms. Organisms with non matching predictions are excluded from the calculation.
# sRNA : only the interaction regions in the mRNA are considered (default). If sRNA is an argument also the interaction regions in the sRNA are considered
# script by Jens Georg

# R --slave -f ../fast_pvalue5.r --args NC_000913 ooi_consensus overall_consensus

ooi_consensus<-FALSE
overall_consensus<-FALSE
mRNA_only<-TRUE

args <- commandArgs(trailingOnly = TRUE) 
ooi <- args[1] 

ooi_consensus1<-na.omit(match("ooi_consensus", args))
if(length(ooi_consensus1)>0){
	ooi_consensus<-TRUE
}

overall_consensus1<-na.omit(match("overall_consensus", args))
if(length(overall_consensus1)>0){
	overall_consensus<-TRUE
}

srna<-na.omit(match("sRNA", args))
if(length(srna)>0){
	mRNA_only<-FALSE
}

mnum=40000		
	
da<-read.csv("CopraRNA2_prep_anno_addhomologs_padj_amountsamp.csv") ## edit prw changed file name 
options <- read.table("CopraRNA_option_file.txt", sep=":") 
root<-as.numeric(as.character(options[14,2]))



en<-grep("Annotation", colnames(da))
da2<-da[,3:(en-1)]
namegenomes<-colnames(da2)
genomes<-list()

myfun<-function(da2){
	int_table<-gsub("^([^\\|]*\\|[^\\|]*\\|)", "",da2)
	int_table<-gsub("\\|.*","",int_table)
	int_table
}


int_table<-apply(da2, 2, myfun)


int_table<-matrix(as.numeric(int_table),nrow(int_table),ncol(int_table))
colnames(int_table)<-colnames(da2)
rownames(int_table)<-seq(1,nrow(int_table))

weight<-read.csv("zscore.weight", header=F, sep=";")

weight<-weight[ match(colnames(da2),toupper(weight[,1])),]


if(root==0){
	weight[,2]<-1
}

if(root>0){
	weight[,2]<-weight[,2]^(1/root)
}


qtrans<-function(dat2){
		if(is(dat2)[1]=="numeric"){
			dat2<-matrix(dat2,1,length(dat2))
		}
        for(i in 1:ncol(dat2)){
                dat2[,i]<-qnorm(dat2[,i], lower.tail=FALSE)
        }
        dat2
}


posvec<-rep(NA,10^7)
vari<-rep(NA,10^7)
count<-0
for(i in 1:nrow(int_table)){
	temp<-order(int_table[i,], na.last=NA)
	#print(i)
	if(length(temp)>2){               
		for(j in 1:(length(temp)-2)){
			count<-count+1
			vari[count]<-paste(sort(temp[1:(j+2)]), collapse="_")
			posvec[count]<-i
		}
	}
}

vari<-na.omit(vari)
posvec<-na.omit(posvec)
dups<-which(duplicated(vari))
first_occurence<-match(vari[dups], vari)
posvec2<-posvec

for(i in 1:length(dups)){
	posvec2[first_occurence[i]]<-paste(posvec2[first_occurence[i]],posvec[dups[i]], sep="_")

}
posvec2<-posvec2[-dups]
vari<-unique(vari)


spa<-floor(length(vari)/mnum)

rest<-length(vari)-(mnum*spa)
spal<-rep(mnum,spa)
if(rest>0){
spa<-spa+1
spal<-c(spal,rest)
}

out2<-rep((list(rep(NA,nrow(int_table)))), spa)
#ooilist<-rep(rep((list(rep(NA,nrow(int_table)))),spa), ncol(int_table))

ooilist<-rep(list(rep((list(rep(NA,nrow(int_table)))), spa)),ncol(int_table))

na<-paste("_",vari,"_",sep="")
count<-1
for(jj in 1:spa){

out<-rep((list(rep(NA,nrow(int_table)))), spal[jj])

for(i in count:(count+spal[jj]-1)){
	
	temp1<-as.numeric(strsplit(vari[i],"_")[[1]])
	temp<-na.omit(int_table[,temp1])
	wtemp<-weight[temp1,2]
	position<-as.numeric(strsplit(posvec2[i],"_")[[1]])
	print(i)
	dat2<-qtrans(temp)
    a<-rowSums(t(t(dat2)*wtemp))
    b<-sum(wtemp)^2
    d<-sum(wtemp^2)
	l<-length(a)
    k<-1/l
	y<-k * seq(1,l)
    dat3<-data.frame(y,a,d,b)
    rhotemp<-coef(nls(y~sort(pnorm(a/sqrt((1-rho)*d+rho*b), lower.tail=F)), data=dat3, start=list(rho=0), control=list(minFactor = 1/128, tol = 1e-05, warnOnly = TRUE)))
	out[[(i-count+1)]][position]<-pnorm(a[match(as.character(position),rownames(temp))]/sqrt((1-rhotemp)*d+rhotemp*b),lower.tail=F)
	
	
}

for(i in 1:length(ooilist)){
	print(i)
	temp<-grep(paste("_",i,"_",sep=""), na[count:(count+spal[jj]-1)])
	#print(temp)
	if(length(temp)>0){
		ooilist[[i]][[jj]]<-do.call(pmin, c(out[temp],list(na.rm=T)))
		#out_ooi[,i]<-do.call(pmin, c(out[temp],list(na.rm=T)))
	}
}

count<-count+spal[jj]
#print(count)
out2[[jj]]<-do.call(pmin, c(out,list(na.rm=T)))

}
out<-out2

out_evo<-do.call(pmin, c(out,list(na.rm=T)))
out_fdr<-p.adjust(out_evo, method="BH")
out_evo<-cbind(out_fdr,out_evo)
colnames(out_evo)<-c("fdr","p-value")
out_evo<-cbind(out_evo, da[,3:ncol(da)])

initial_sorting<-seq(1,nrow(int_table))
out_evo<-cbind(out_evo, initial_sorting)
out_evo<-out_evo[order(as.numeric(out_evo[,2])),]

write.table(out_evo, file="CopraRNA2_final_all_balanced.csv",sep=",", quote=F, row.names=F)


out_ooi<-int_table
out_ooi[]<-NA

for(i in 1:ncol(out_ooi)){
	#print(i)
	#temp<-grep(paste("_",i,"_",sep=""), na)
	#out_ooi[,i]<-do.call(pmin, c(out[temp],list(na.rm=T)))
	out_ooi[,i]<-do.call(pmin, c(ooilist[[i]],list(na.rm=T)))
}

ooi_pos<-grep(ooi, colnames(int_table))




out_ooi_res<-out_ooi[,ooi_pos]
out_ooi_fdr<-p.adjust(out_ooi_res, method="BH")
out_ooi_res<-cbind(out_ooi_fdr,out_ooi_res)
colnames(out_ooi_res)<-c("fdr","p-value")
out_ooi_res<-cbind(out_ooi_res, da[,3:ncol(da)])

initial_sorting<-seq(1,nrow(int_table))
out_ooi_res<-cbind(out_ooi_res, initial_sorting)
out_ooi_res<-out_ooi_res[order(as.numeric(out_ooi_res[,2])),]



write.table(out_ooi_res, file="CopraRNA2_final_all_ooi.csv",sep=",", quote=F, row.names=F)



out_rho_evo<-out_ooi
out_evo_rank<-out_rho_evo
top<-c()
rank_thres<-200

anno<-da

 for(i in 1:ncol(out_evo_rank)){
# print(i)
	 out_evo_rank[,i]<-rank(as.numeric(out_rho_evo[,i]), na.last="keep" ,ties.method="first")
	 temp<-which(as.numeric(out_evo_rank[,i])<=rank_thres)
	 top<-c(top, temp)
 }

 top<-unique(top)

initial_sorting<-initial_sorting[top]
 evo_ooi_pvalue<-out_rho_evo[top,]
 evo_ooi_rank<-out_evo_rank[top,]
 evo_int_pvalue<-int_table[top,]
 evo_anno<-anno[top,3:ncol(anno)]


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


Rank_p_value<-rank_list2
Rank_p_value<-Rank_p_value[,2]
evo_analysis<-cbind(Rank_p_value,evo_ooi_rank_result,evo_ooi_p_result,evo_int_p_result,evo_anno,evo_ooi_rank,evo_ooi_pvalue,evo_int_pvalue,initial_sorting)
write.table(evo_analysis, file="CopraRNA2_final_all_evo_table2.csv", sep="\t", row.names=F)

evo_int_pvalue<-matrix(as.numeric(evo_int_pvalue),nrow(evo_int_pvalue),ncol(evo_int_pvalue))
evo_int_pvalue_scored<-matrix(as.numeric(evo_int_pvalue_scored),nrow(evo_int_pvalue_scored),ncol(evo_int_pvalue_scored))



################## ooi based consensus dependent pValue combination ################

if(overall_consensus==TRUE){

load("conservation_table.Rdata")


con_table<-conservation_table[[1]]
con_table_sub<-conservation_table[[2]]

p_table<-con_table
p_table[]<-NA

mRNA_test<-function(x){
	out<-NA
	if(is.na(x[1])==F){
		if(length(x)==4){
			if(x[2]=="TRUE"){
				out<-as.numeric(x[1])
			}
		}
	}
	out
}
mRNA_and_sRNA_test<-function(x){
	out<-NA
	if(is.na(x[1])==F){
	if(length(x)==4){
		if(x[2]=="TRUE"){
			if(x[3]=="TRUE" | x[4]=="TRUE"){
				out<-as.numeric(x[1])
			}
		}
	}
	}
	out
}


for(i in 1:nrow(p_table)){
	
	if(mRNA_only==TRUE){
		con_temp<-con_table[i,]
		con_temp<-strsplit(con_temp, "\\|")
		con_temp_sub<-con_table_sub[i,]
		con_temp_sub<-strsplit(con_temp_sub, "\\|")
		con_temp<-unlist(lapply(con_temp,mRNA_test))
		con_temp_sub<-unlist(lapply(con_temp_sub,mRNA_test))
		na<-which(is.na(con_temp))
		if(length(na)>0){
			con_temp[na]<-con_temp_sub[na]
		}
		
		p_table[i,]<-con_temp
	}
	
	if(mRNA_only==FALSE){
		con_temp<-con_table[i,]
		con_temp<-strsplit(con_temp, "\\|")
		con_temp_sub<-con_table_sub[i,]
		con_temp_sub<-strsplit(con_temp_sub, "\\|")
		con_temp<-unlist(lapply(con_temp,mRNA_and_sRNA_test))
		con_temp_sub<-unlist(lapply(con_temp_sub,mRNA_and_sRNA_test))
		na<-which(is.na(con_temp))
		if(length(na)>0){
			con_temp[na]<-con_temp_sub[na]
		}
		
		p_table[i,]<-con_temp
	}
}
nan<-colnames(p_table)
p_table<-matrix(as.numeric(p_table),nrow(p_table),ncol(p_table))
colnames(p_table)<-nan

rownames(p_table)<-seq(1,nrow(p_table))

posvec<-rep(NA,10^7)
vari<-rep(NA,10^7)
count<-0
for(i in 1:nrow(p_table)){
	temp<-order(p_table[i,], na.last=NA)
	#print(i)
	if(length(temp)>2){               ##############<----hier 2 auf 3
		for(j in 1:(length(temp)-2)){
			count<-count+1
			vari[count]<-paste(sort(temp[1:(j+2)]), collapse="_")
			posvec[count]<-i
		}
	}
}

vari<-na.omit(vari)
posvec<-na.omit(posvec)
dups<-which(duplicated(vari))
first_occurence<-match(vari[dups], vari)
posvec2<-posvec

for(i in 1:length(dups)){
	posvec2[first_occurence[i]]<-paste(posvec2[first_occurence[i]],posvec[dups[i]], sep="_")

}
posvec2<-posvec2[-dups]
vari<-unique(vari)

# out<-rep((list(rep(NA,nrow(p_table)))), length(vari))

# for(i in 1:length(vari)){
	
	# temp1<-as.numeric(strsplit(vari[i],"_")[[1]])
	# temp<-na.omit(int_table[,temp1])
	# wtemp<-weight[temp1,2]
	# position<-as.numeric(strsplit(posvec2[i],"_")[[1]])
	# print(i)
	# dat2<-qtrans(temp)
    # a<-rowSums(t(t(dat2)*wtemp))
    # b<-sum(wtemp)^2
    # d<-sum(wtemp^2)
	# l<-length(a)
    # k<-1/l
	# y<-k * seq(1,l)
    # dat3<-data.frame(y,a,d,b)
    
	# rhotemp<-coef(nls(y~sort(pnorm(a/sqrt((1-rho)*d+rho*b), lower.tail=F)), data=dat3, start=list(rho=0), control=list(minFactor = 1/128, tol = 1e-05, warnOnly = TRUE)))
	
	# temp<-na.omit(p_table[position,temp1])
	# wtemp<-weight[temp1,2]
	# position<-as.numeric(strsplit(posvec2[i],"_")[[1]])
	
	# dat2<-qtrans(temp)
    # a<-rowSums(t(t(dat2)*wtemp))
    # b<-sum(wtemp)^2
    # d<-sum(wtemp^2)
	
	# out[[i]][position]<-pnorm(a/sqrt((1-rhotemp)*d+rhotemp*b),lower.tail=F)

# }

################################

spa<-floor(length(vari)/mnum)

rest<-length(vari)-(mnum*spa)
spal<-rep(mnum,spa)
if(rest>0){
spa<-spa+1
spal<-c(spal,rest)
}

out2<-rep((list(rep(NA,nrow(int_table)))), spa)
#ooilist<-rep(rep((list(rep(NA,nrow(int_table)))),spa), ncol(int_table))

ooilist<-rep((list(rep(NA,nrow(int_table)))), spa)

na<-paste("_",vari,"_",sep="")
count<-1
ooi_pos<-grep(ooi, colnames(p_table))
for(jj in 1:spa){

out<-rep((list(rep(NA,nrow(int_table)))), spal[jj])

for(i in count:(count+spal[jj]-1)){
	
	temp1<-as.numeric(strsplit(vari[i],"_")[[1]])
	temp<-na.omit(int_table[,temp1])
	wtemp<-weight[temp1,2]
	position<-as.numeric(strsplit(posvec2[i],"_")[[1]])
	print(i)
	dat2<-qtrans(temp)
    a<-rowSums(t(t(dat2)*wtemp))
    b<-sum(wtemp)^2
    d<-sum(wtemp^2)
	l<-length(a)
    k<-1/l
	y<-k * seq(1,l)
    dat3<-data.frame(y,a,d,b)
    
	rhotemp<-coef(nls(y~sort(pnorm(a/sqrt((1-rho)*d+rho*b), lower.tail=F)), data=dat3, start=list(rho=0), control=list(minFactor = 1/128, tol = 1e-05, warnOnly = TRUE)))
	
	temp<-na.omit(p_table[position,temp1])
	wtemp<-weight[temp1,2]
	position<-as.numeric(strsplit(posvec2[i],"_")[[1]])
	
	dat2<-qtrans(temp)
    a<-rowSums(t(t(dat2)*wtemp))
    b<-sum(wtemp)^2
    d<-sum(wtemp^2)
	
	out[[i-count+1]][position]<-pnorm(a/sqrt((1-rhotemp)*d+rhotemp*b),lower.tail=F)
	
	
}

#for(i in 1:length(ooilist)){
	#print(i)
	
	temp<-grep(paste("_",ooi_pos,"_",sep=""), na[count:(count+spal[jj]-1)])
	#print(temp)
	if(length(temp)>0){
	ooilist[[jj]]<-do.call(pmin, c(out[temp],list(na.rm=T)))
	#out_ooi[,i]<-do.call(pmin, c(out[temp],list(na.rm=T)))
	}
#}

count<-count+spal[jj]
#print(count)
out2[[jj]]<-do.call(pmin, c(out,list(na.rm=T)))

}
out<-out2

################################




out_evo<-do.call(pmin, c(out,list(na.rm=T)))
out_fdr<-p.adjust(out_evo, method="BH")
out_evo<-cbind(out_fdr,out_evo)
colnames(out_evo)<-c("fdr","p-value")
out_evo<-cbind(out_evo, da[,3:ncol(da)])

initial_sorting<-seq(1,nrow(p_table))
out_evo<-cbind(out_evo, initial_sorting)
out_evo<-out_evo[order(as.numeric(out_evo[,2])),]



#ooi_pos<-grep(ooi, colnames(p_table))
#na<-paste("_",vari,"_",sep="")
#temp<-grep(paste("_",ooi_pos,"_",sep=""), na)
out_ooi<-do.call(pmin, c(ooilist,list(na.rm=T)))
out_ooi_fdr<-p.adjust(out_ooi, method="BH")
out_ooi<-cbind(out_ooi_fdr,out_ooi)
colnames(out_ooi)<-c("fdr","p-value")
out_ooi<-cbind(out_ooi, da[,3:ncol(da)])

initial_sorting<-seq(1,nrow(p_table))
out_ooi<-cbind(out_ooi, initial_sorting)
out_ooi<-out_ooi[order(as.numeric(out_ooi[,2])),]



write.table(out_evo, file="CopraRNA2_final_all_balanced_consensus.csv",sep=",", quote=F, row.names=F) ## edit prw changed file
write.table(out_ooi, file="CopraRNA2_final_all_ooi_consensus.csv",sep=",", quote=F, row.names=F) ## edit prw changed file

}


################## ooi consensus dependent pValue combination ################

if(ooi_consensus==TRUE){

load("conservation_table.Rdata")



con_table<-conservation_table[[3]]
con_table_sub<-conservation_table[[4]]

p_table<-con_table
p_table[]<-NA

mRNA_test<-function(x){
	out<-NA
	if(is.na(x[1])==F){
		if(length(x)==4){
			if(x[2]=="TRUE"){
				out<-as.numeric(x[1])
			}
		}
	}
	out
}
mRNA_and_sRNA_test<-function(x){
	out<-NA
	if(is.na(x[1])==F){
	if(length(x)==4){
		if(x[2]=="TRUE"){
			if(x[3]=="TRUE" | x[4]=="TRUE"){
				out<-as.numeric(x[1])
			}
		}
	}
	}
	out
}


for(i in 1:nrow(p_table)){
	
	if(mRNA_only==TRUE){
		con_temp<-con_table[i,]
		con_temp<-strsplit(con_temp, "\\|")
		con_temp_sub<-con_table_sub[i,]
		con_temp_sub<-strsplit(con_temp_sub, "\\|")
		con_temp<-unlist(lapply(con_temp,mRNA_test))
		con_temp_sub<-unlist(lapply(con_temp_sub,mRNA_test))
		na<-which(is.na(con_temp))
		if(length(na)>0){
			con_temp[na]<-con_temp_sub[na]
		}
		
		p_table[i,]<-con_temp
	}
	
	if(mRNA_only==FALSE){
		con_temp<-con_table[i,]
		con_temp<-strsplit(con_temp, "\\|")
		con_temp_sub<-con_table_sub[i,]
		con_temp_sub<-strsplit(con_temp_sub, "\\|")
		con_temp<-unlist(lapply(con_temp,mRNA_and_sRNA_test))
		con_temp_sub<-unlist(lapply(con_temp_sub,mRNA_and_sRNA_test))
		na<-which(is.na(con_temp))
		if(length(na)>0){
			con_temp[na]<-con_temp_sub[na]
		}
		
		p_table[i,]<-con_temp
	}
}
nan<-colnames(p_table)
p_table<-matrix(as.numeric(p_table),nrow(p_table),ncol(p_table))
colnames(p_table)<-nan


posvec<-rep(NA,10^7)
vari<-rep(NA,10^7)
count<-0
for(i in 1:nrow(p_table)){
	temp<-order(p_table[i,], na.last=NA)
	
	if(length(temp)>2){              
		for(j in 1:(length(temp)-2)){
			count<-count+1
			vari[count]<-paste(sort(temp[1:(j+2)]), collapse="_")
			posvec[count]<-i
		}
	}
}

vari<-na.omit(vari)
posvec<-na.omit(posvec)
dups<-which(duplicated(vari))
first_occurence<-match(vari[dups], vari)
posvec2<-posvec

for(i in 1:length(dups)){
	posvec2[first_occurence[i]]<-paste(posvec2[first_occurence[i]],posvec[dups[i]], sep="_")

}
posvec2<-posvec2[-dups]
vari<-unique(vari)


spa<-floor(length(vari)/mnum)

rest<-length(vari)-(mnum*spa)
spal<-rep(mnum,spa)
if(rest>0){
spa<-spa+1
spal<-c(spal,rest)
}

#out2<-rep((list(rep(NA,nrow(int_table)))), spa)
#ooilist<-rep(rep((list(rep(NA,nrow(int_table)))),spa), ncol(int_table))

ooilist<-rep((list(rep(NA,nrow(int_table)))), spa)

na<-paste("_",vari,"_",sep="")
count<-1
ooi_pos<-grep(ooi, colnames(p_table))
for(jj in 1:spa){

out<-rep((list(rep(NA,nrow(int_table)))), spal[jj])

for(i in count:(count+spal[jj]-1)){
	
	temp1<-as.numeric(strsplit(vari[i],"_")[[1]])
	temp<-na.omit(int_table[,temp1])
	wtemp<-weight[temp1,2]
	position<-as.numeric(strsplit(posvec2[i],"_")[[1]])
	print(i)
	dat2<-qtrans(temp)
    a<-rowSums(t(t(dat2)*wtemp))
    b<-sum(wtemp)^2
    d<-sum(wtemp^2)
	l<-length(a)
    k<-1/l
	y<-k * seq(1,l)
    dat3<-data.frame(y,a,d,b)
    
	rhotemp<-coef(nls(y~sort(pnorm(a/sqrt((1-rho)*d+rho*b), lower.tail=F)), data=dat3, start=list(rho=0), control=list(minFactor = 1/128, tol = 1e-05, warnOnly = TRUE)))
	
	temp<-na.omit(p_table[position,temp1])
	wtemp<-weight[temp1,2]
	position<-as.numeric(strsplit(posvec2[i],"_")[[1]])
	
	dat2<-qtrans(temp)
    a<-rowSums(t(t(dat2)*wtemp))
    b<-sum(wtemp)^2
    d<-sum(wtemp^2)
	
	out[[i-count+1]][position]<-pnorm(a/sqrt((1-rhotemp)*d+rhotemp*b),lower.tail=F)
	
	
}

#for(i in 1:length(ooilist)){
	#print(i)
	
	temp<-grep(paste("_",ooi_pos,"_",sep=""), na[count:(count+spal[jj]-1)])
	#print(temp)
	if(length(temp)>0){
	ooilist[[jj]]<-do.call(pmin, c(out[temp],list(na.rm=T)))
	#out_ooi[,i]<-do.call(pmin, c(out[temp],list(na.rm=T)))
	}
#}

count<-count+spal[jj]
#print(count)
#out2[[jj]]<-do.call(pmin, c(out,list(na.rm=T)))

}
#out<-out2

################################




# out_evo<-do.call(pmin, c(out,list(na.rm=T)))
# out_fdr<-p.adjust(out_evo, method="BH")
# out_evo<-cbind(out_fdr,out_evo)
# colnames(out_evo)<-c("fdr","p-value")
# out_evo<-cbind(out_evo, da[,3:ncol(da)])

# initial_sorting<-seq(1,nrow(p_table))
# out_evo<-cbind(out_evo, initial_sorting)
# out_evo<-out_evo[order(as.numeric(out_evo[,2])),]



#ooi_pos<-grep(ooi, colnames(p_table))
#na<-paste("_",vari,"_",sep="")
#temp<-grep(paste("_",ooi_pos,"_",sep=""), na)
out_ooi<-do.call(pmin, c(ooilist,list(na.rm=T)))

# out<-rep((list(rep(NA,nrow(p_table)))), length(vari))



# for(i in 1:length(vari)){
	
	# temp1<-as.numeric(strsplit(vari[i],"_")[[1]])
	# temp<-na.omit(int_table[,temp1])
	# wtemp<-weight[temp1,2]
	# position<-as.numeric(strsplit(posvec2[i],"_")[[1]])
	# print(i)
	# dat2<-qtrans(temp)
    # a<-rowSums(t(t(dat2)*wtemp))
    # b<-sum(wtemp)^2
    # d<-sum(wtemp^2)
	# l<-length(a)
    # k<-1/l
	# y<-k * seq(1,l)
    # dat3<-data.frame(y,a,d,b)
    
	# rhotemp<-coef(nls(y~sort(pnorm(a/sqrt((1-rho)*d+rho*b), lower.tail=F)), data=dat3, start=list(rho=0), control=list(minFactor = 1/128, tol = 1e-05, warnOnly = TRUE)))
	
	# temp<-na.omit(p_table[position,temp1])
	
	# wtemp<-weight[temp1,2]
	# position<-as.numeric(strsplit(posvec2[i],"_")[[1]])
	
	# dat2<-qtrans(temp)
    # a<-rowSums(t(t(dat2)*wtemp))
    # b<-sum(wtemp)^2
    # d<-sum(wtemp^2)
	
	# out[[i]][position]<-pnorm(a/sqrt((1-rhotemp)*d+rhotemp*b),lower.tail=F)

# }

# ooi_pos<-grep(ooi, colnames(p_table))
# na<-paste("_",vari,"_",sep="")
# temp<-grep(paste("_",ooi_pos,"_",sep=""), na)
# out_ooi<-do.call(pmin, c(out[temp],list(na.rm=T)))
out_ooi_fdr<-p.adjust(out_ooi, method="BH")
out_ooi<-cbind(out_ooi_fdr,out_ooi)
colnames(out_ooi)<-c("fdr","p-value")
out_ooi<-cbind(out_ooi, da[,3:ncol(da)])

initial_sorting<-seq(1,nrow(p_table))
out_ooi<-cbind(out_ooi, initial_sorting)
out_ooi<-out_ooi[order(as.numeric(out_ooi[,2])),]

write.table(out_ooi, file="CopraRNA2_final_all_ooi_ooiconsensus.csv",sep=",", quote=F, row.names=F)

}

