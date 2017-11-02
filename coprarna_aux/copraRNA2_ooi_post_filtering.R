#R --slave -f ../copraRNA2_ooi_post_filtering.r --args ooi=NC_000913 thres=0.45 

# post-filter CopraRNA results for predictions with a p-value >= a given threshold (default: 0.45) in the ooi and/or not hitting the consensus interaction site in the mRNA

#thres<-0.45
#ooi<-"NC_000913"
inputfile<-"CopraRNA2_final_all_ooi.csv"
inputfile2<-"CopraRNA2_final_all_ooi_ooiconsensus.csv"
inputfile3<-"CopraRNA2_final_all_ooi_consensus.csv"
args <- commandArgs(trailingOnly = TRUE) 
 for(i in 1:length(args)){
	temp<-strsplit(args[i],"=")
	temp<-temp[[1]]
	temp1<-temp[1]
	temp2<-temp[2]
	assign(as.character(temp1),temp2)
 }

thres<-as.numeric(thres)

### filter for p-value ##############
dat<-read.csv(inputfile,sep=",", header=T)
dat<-dat[order(dat[,"initial_sorting"]),]
pos<-grep(ooi, colnames(dat))
dat2<-dat[,pos]
dat2<-as.matrix(dat2)
fun<-function(x){
	out<-as.numeric(strsplit(x, "\\|")[[1]][3])
	out
}
dat3<-unlist(lapply(dat2,fun))
negthres<-which(dat3>= thres)



dat_2<-read.csv(inputfile2,sep=",", header=T)
dat_2<-dat_2[order(dat_2[,"initial_sorting"]),]
pos<-grep(ooi, colnames(dat_2))
dat_22<-dat_2[,pos]
dat_22<-as.matrix(dat_22)
fun<-function(x){
	out<-as.numeric(strsplit(x, "\\|")[[1]][3])
	out
}
dat_3<-unlist(lapply(dat_22,fun))
negthres2<-which(dat_3>= thres)


dat_3<-read.csv(inputfile3,sep=",", header=T)
dat_3<-dat_3[order(dat_3[,"initial_sorting"]),]
pos<-grep(ooi, colnames(dat_3))
dat_222<-dat_3[,pos]
dat_222<-as.matrix(dat_222)
fun<-function(x){
	out<-as.numeric(strsplit(x, "\\|")[[1]][3])
	out
}
dat_4<-unlist(lapply(dat_222,fun))
negthres3<-which(dat_4>= thres)



#### filter for consensensus intersction site ###########
load("conservation_table.Rdata")
con_table<-conservation_table[[1]]
con_table_sub<-conservation_table[[2]]

pos<-grep(ooi, colnames(con_table))

opt<-con_table[,pos]
subopt<-con_table_sub[,pos]

negpos<-c()
for(i in 1:length(opt)){
	if(is.na(opt[i])==F){
		temp_opt<-strsplit(opt[i], "\\|")[[1]]
		temp_sub<-strsplit(subopt[i], "\\|")[[1]]
		if(temp_opt[2]=="FALSE"){
			if(is.na(subopt[i])==F){
				if(temp_sub[2]=="FALSE"){
					negpos<-c(negpos,i)
				}
				if(temp_sub[2]=="TRUE" & as.numeric(temp_sub[1])>=thres){
					negpos<-c(negpos,i)
				}
			}
			if(is.na(subopt[i])==T){
				negpos<-c(negpos,i)
			}
		}
	}
}

####################################


if(length(negthres)>0){
	dat1<-dat[-negthres,]
}
if(length(negthres2)>0){
	dat_2<-dat_2[-negthres2,]
}

if(length(negthres3)>0){
	dat_3<-dat_3[-negthres3,]
}
neg<-union(negthres,negpos)

if(length(neg)>0){
	dat2<-dat[-neg,]
}

dat1<-dat1[order(dat1[,"p.value"]),]
dat2<-dat2[order(dat2[,"p.value"]),]
dat_2<-dat_2[order(dat_2[,"p.value"]),]
dat_3<-dat_3[order(dat_3[,"p.value"]),]

write.table(dat1, file="CopraRNA2_final_all_ooi_filtered1.csv",sep=",", quote=F, row.names=F)
write.table(dat2, file="CopraRNA2_final_all_ooi_filtered2.csv",sep=",", quote=F, row.names=F)
write.table(dat_2, file="CopraRNA2_final_all_ooi_ooi_consensus_filtered.csv",sep=",", quote=F, row.names=F)
write.table(dat_3, file="CopraRNA2_final_all_ooi_consensus_filtered.csv",sep=",", quote=F, row.names=F)
