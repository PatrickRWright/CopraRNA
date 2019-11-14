# script by Jens Georg

# dependencies:

## Tools: 
### Mafft
### dialign-tx

## Files:
### CopraRNA_available_organisms.txt
### path to dialign config file
### jalview_props.txt

#call:
# R --slave -f ./copraRNA2_find_conserved_sites.r 


suppressPackageStartupMessages(require(phangorn))
suppressPackageStartupMessages(require(seqinr))
suppressPackageStartupMessages(require(doMC))

# get absolute path
initial.options <- commandArgs(trailingOnly = FALSE)
path<-initial.options [4]
path<-sub("copraRNA2_find_conserved_sites.r","",path)
#print(path)

# preset path to required files, path can also be specified as argument
copref_path<-paste(path,"CopraRNA_available_organisms.txt",sep="")
dialign_conf<-paste(path,"dialign_conf/",sep="")

# read the organism of interest (ooi) from the ncRNA fasta file. The sRNA of the ooi is considered to be the first sequence.
ooi<-gsub("ncRNA_","",names(read.fasta("ncrna.fa"))[1])

# register cores for parallel processing
co<-readLines("CopraRNA_option_file.txt") 
max_cores<-as.numeric(gsub("core count:","",co[grep("core count:", co)]))
registerDoMC(max_cores)

# number of top predictions which should be investigated
top<-as.numeric(gsub("top count:","",co[grep("top count:", co)]))

# IntaRNA parameters
winsize<-as.numeric(gsub("win size:","",co[grep("win size:", co)]))

maxbpdist<-as.numeric(gsub("max bp dist:","",co[grep("max bp dist:", co)]))

maxbpdist<-as.numeric(gsub("max bp dist:","",co[grep("max bp dist:", co)]))

temperature<-as.numeric(gsub("temperature:","",co[grep("temperature:", co)]))

# method to calculate pyhlogentic weights from the 16S alignment. "clustal" = ClustalW method, "copra" = CopraRNA_1 method
weight_method="clustal"

# path to ribosomal RNA fasta
ribosomal_rna="16s_sequences.fa"


# path to CopraRNA result file (in order to investigate the top predictions)
copra_result<-"CopraRNA_result_all.csv"




# perform Mafft alignment if not already present (default alignments are present)
own_alignment=FALSE


conservation_oois=ooi

# transforming arguments in valid variables 
args <- commandArgs(trailingOnly = TRUE) 
if(length(args)>0){
	for(i in 1:length(args)){
		temp<-strsplit(args[i],"=")
		temp<-temp[[1]]
		temp1<-temp[1]
		temp2<-temp[2]
		assign(as.character(temp1),temp2)
	}
}
top<-as.numeric(top)

# function to create alignment annotation files for alignment drawing/visualization with jalview
jalview_anno<-function(tab_aligned, tabsub_aligned,peaks1,test, ooi=conservation_oois, nam2,alignment,sRNA_alignment2 ,all_orgs,name="srna", type="png"){
	peaks_all<-test[[3]]
	if(length(peaks_all)>0){
		peaks<-unlist(lapply(peaks1, function(x){return(x[1])}))
		peaks<-sort(peaks)
		tab_aligned[,"cluster_id"]<-unlist(lapply(tab_aligned[,"cluster_id"] ,function(x){return(strsplit(as.character(x), split="\\|")[[1]][1])}))
		tabsub_aligned[,"cluster_id"]<-unlist(lapply(tabsub_aligned[,"cluster_id"] ,function(x){return(strsplit(as.character(x), split="\\|")[[1]][1])}))
		colo2<-c("#caff70","#7aa9df","#cf97bb" ,"#fede7e","#fc9187","#d3648f","#4c9998","#988377","#576e81","#9fe1e5")
		colo<-c("#799943","#5680b0","#9c6488","#e4c054","#e26a5f","#b43768","#006362","#614736","#294258","#77b1b5")
		colo3<-c()
		for(i in 1:length(peaks1)){
			fc <- colorRampPalette(c(colo[i], colo2[i]))
			tmp<-fc(length(peaks1[[i]]))
			names(tmp)<-peaks1[[i]]
			colo3<-c(colo3,tmp)
		}
		colo3<-c(colo3,"#b0b0b0")
		names(colo3)[length(colo3)]<-"not_conserved"
		peaks<-names(colo3)
	} else {
		peaks<-"not_conserved"
		colo3<-"#b0b0b0"
		names(colo3)<-"not_conserved"
	}
	ooi_pos<-grep(ooi[1], tab_aligned[,"orgs"])
	if(length(ooi_pos)>0){
		ooi_int<-as.character(tab_aligned[ooi_pos,"cluster_id"])
		peaks<-unlist(unique(c(ooi_int,peaks)))
		colo4<-colo3
		colo3<-gsub("#","",unique(c(colo3[ooi_int],colo3)))
	}
	na<-which(is.na(peaks))
	if(length(na)>0){
		peaks<-peaks[-na]
		colo3<-colo3[-na]
	}
	anno_mRNA<-paste(peaks, "\t", colo3, sep="")
	anno_sRNA<-anno_mRNA
	tab<-rbind(tab_aligned,tabsub_aligned)
	an<-match(tab[,"orgs"],all_orgs)
	an<-nam2[an]
	ntab<-paste(an,tab[,"name2"],tab[,"orgs"],sep="/")
	na<-which(is.na(tab[,"cluster_id"]))
	tab2<-NA
	tab3<-tab[order(unlist(tab[,"cluster_id"]),na.last = TRUE),]
	ord<-match(unique(tab3[,"name2"]),names(alignment))
	if(length(na)>0){
		tab2 <-tab[na,]
		if(is.matrix(tab2)==F){
			tab2<-t(as.matrix(tab2))
		}
		ntab2<-ntab[na]
		tab<-tab[-na,]
		ntab<-ntab[-na]
		if(is.matrix(tab)==F){
			tab<-t(as.matrix(tab))
		}
	}
	an<-match(tab_aligned[,"orgs"],all_orgs)
	an<-nam2[an]
	names(alignment)<-paste(an,names(alignment),tab_aligned[,"orgs"],sep="/")
	names(sRNA_alignment2)<-names(alignment)
	jal_anno_m<-"JALVIEW_ANNOTATION"
	jal_anno_s<-"JALVIEW_ANNOTATION"
	gr<-names(alignment)[match(ooi, gsub(".*/","",names(alignment)))]
	jal_anno_m<-c(jal_anno_m, paste("SEQUENCE_GROUP","oois", 1,length(alignment[[1]]),-1,paste(gr, collapse="\t"),sep="\t"  ))
	jal_anno_s<-c(jal_anno_s, paste("SEQUENCE_GROUP","oois", 1,length(sRNA_alignment2[[1]]),-1,paste(gr, collapse="\t"),sep="\t"  ))
	jal_anno_m<-c(jal_anno_m,paste("PROPERTIES","oois","colour=None","outlineColour=black","displayBoxes=true","displayText=true","idColour=FFF8DC", sep="\t"))
	jal_anno_s<-c(jal_anno_s,paste("PROPERTIES","oois","colour=None","outlineColour=black","displayBoxes=true","displayText=true","idColour=FFF8DC", sep="\t"))
	if(nrow(tab)>0){
		for(i in 1:nrow(tab)){
			pos<-as.numeric(strsplit(as.character(tab[i,"int_mrna2"]),",")[[1]])
			pos1<-diff(pos)
			s<-c(0,which(pos1>1),length(pos))
			for(j in 1:(length(s)-1)){
				tmp<-paste(as.numeric(tab[i,"Energy"])*-1, "\t", ntab[i], "\t-1\t", pos[s[j]+1], "\t", pos[s[j+1]],"\t",tab[i,"cluster_id"], sep="")
				anno_mRNA<-c(anno_mRNA,tmp)
			}
		}
	}
	if(is.matrix(tab2)==T){
		for(i in 1:nrow(tab2)){
			pos<-as.numeric(strsplit(as.character(tab2[i,"int_mrna2"]),",")[[1]])
			pos1<-diff(pos)
			s<-c(0,which(pos1>1),length(pos))
			for(j in 1:(length(s)-1)){
				tmp<-paste(as.numeric(tab2[i,"Energy"])*-1, "\t", ntab2[i], "\t-1\t", pos[s[j]+1], "\t", pos[s[j+1]],"\t","not_conserved", sep="")
				anno_mRNA<-c(anno_mRNA,tmp)
			}
		}
	}
	alignment<-alignment[ord]
	write.table(anno_mRNA, file=paste(name,"_mRNA_features.txt", sep=""), quote=F, row.names=F, col.names=F)
	write.fasta(alignment, names=names(alignment), file.out=paste(name,"_mRNA_alignment.fasta", sep=""))
	if(nrow(tab)>0){
		for(i in 1:nrow(tab)){
			pos<-as.numeric(strsplit(as.character(tab[i,"int_srna2"]),",")[[1]])
			pos1<-diff(pos)
			s<-c(0,which(pos1>1),length(pos))
			for(j in 1:(length(s)-1)){
				tmp<-paste(as.numeric(tab[i,"Energy"])*-1, "\t",  ntab[i], "\t-1\t", pos[s[j]+1], "\t", pos[s[j+1]],"\t",tab[i,"cluster_id"], sep="")
				anno_sRNA<-c(anno_sRNA,tmp)
			}
		}
	}
	if(is.matrix(tab2)==T){
		if(nrow(tab)>0){
			for(i in 1:nrow(tab2)){
				pos<-as.numeric(strsplit(as.character(tab2[i,"int_srna2"]),",")[[1]])
				pos1<-diff(pos)
				s<-c(0,which(pos1>1),length(pos))
				for(j in 1:(length(s)-1)){
					tmp<-paste(as.numeric(tab2[i,"Energy"])*-1, "\t",  ntab2[i], "\t-1\t", pos[s[j]+1], "\t", pos[s[j+1]],"\t","not_conserved", sep="")
					anno_sRNA<-c(anno_sRNA,tmp)
				}
			}
		}
	}
	sRNA_alignment2<-sRNA_alignment2[ord]
	write.table(anno_sRNA, file=paste(name,"_sRNA_features.txt", sep=""), quote=F, row.names=F, col.names=F)
	write.fasta(sRNA_alignment2, names=names(sRNA_alignment2), file.out=paste(name,"_sRNA_alignment.fasta", sep=""))
	write.table(jal_anno_m, file=paste(name,"_mRNA_annotation.txt", sep=""), quote=F, row.names=F, col.names=F)
	write.table(jal_anno_s, file=paste(name,"_sRNA_annotation.txt", sep=""), quote=F, row.names=F, col.names=F)
}

# function for density based clustering (dbscan-like) based on a similarity matrix
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

# combine all UTR sequences in one fasta file (needed if alignmenta are not available)
fast<-function(){
	falist1<-dir()
	int<-grep("intarna", falist1)
	falist1<-falist1[-int]
	falist<-grep("_([0-9]){1,}_down_([0-9]){1,}\\.fa$",falist1)
	
	
	for(i in 1:length(falist)){
		temp<-read.fasta(falist1[falist[i]], strip.desc = TRUE)
		write.fasta(temp, names=names(temp),file="utr_seqs.fa", open="a")
		}
	}

# extract locus_tags from all UTRs
fast2<-function(){
	falist1<-dir()
	int<-grep("intarna", falist1)
	falist1<-falist1[-int]
	falist<-grep("_([0-9]){1,}_down_([0-9]){1,}\\.fa$",falist1)
	namelist<-vector("list",length(falist))
	name<-falist1[falist]
	name<-strsplit(name,"_")
	name<-unlist(lapply(name, function(x){paste(x[1],x[2],sep="_")}))
	names(namelist)<-name
	for(i in 1:length(falist)){
		tmp<-readLines(falist1[falist[i]])
		header<-grep(">", tmp)
		header<-gsub(">","",tmp[header])
		namelist[[i]]<-header
		}
	namelist	
	}

seqle <- function(x,incr=1) { 
  if(!is.numeric(x)) x <- as.numeric(x) 
  n <- length(x)  
  y <- x[-1L] != x[-n] + incr 
  i <- c(which(y|is.na(y)),n) 
  list(lengths = diff(c(0L,i)),
       values = x[head(c(0L,i)+1L,-1L)]) 
} 

# extrat sub-alignment based on given alignment positions fro alignment in fasta format
cutalign<-function(x,s,e){
	x<-x[s:e]
	x
}

# function to remove gaps from a alignment in fasta_format		
remove_gaps<-function(x){
			y<-grep("-",x)
			if(length(y)>0){
				x<-x[-y]
			}
			x
		}
		
# function to remove nucleotides that are within the borders of the predicted interaction site, but do not interact
exclude_nonbound<-function(tab_aligned){
	outm<-c()
	outs<-c()
	for(i in 1:nrow(tab_aligned)){
		temp<-strsplit(as.character(tab_aligned[i,"hybridDP"]),"&")[[1]]
		m<-paste(which(strsplit(temp[1],"")[[1]]=="(")+as.numeric(tab_aligned[i,1])-1,collapse=",")
		s<-paste(which((strsplit(temp[2],"")[[1]])==")")+as.numeric(tab_aligned[i,5])-1,collapse=",")
		outm<-c(outm,m)
		outs<-c(outs,s)
	}
	out<-list(outm,outs)
	out
}
	
# function to modify the borders of conserved interaction sites based interactions assigned to the site
# uses longest and shortes interaction to define borders		
refine_site<-function(out_comb, tab_comb){
	sites1<-out_comb
	sites<-na.omit(unique(unlist(sites1)))
	na<-grep("NA",sites)
	if(length(na)>0){
		sites<-sites[-na]
	}
	for(i in 1:length(sites)){
		memb<-grep(sites[i], sites1)
		pos<-as.numeric(strsplit(sites[i],"_")[[1]])
		if(min(as.numeric(tab_comb[memb,1]))>pos[1]){
			pos[1]<-min(as.numeric(tab_comb[memb,1]))
		}
		if(max(as.numeric(tab_comb[memb,2]))<pos[2]){
			pos[2]<-max(as.numeric(tab_comb[memb,2]))
		}
		
		if(min(as.numeric(tab_comb[memb,5]))>pos[3]){
			pos[3]<-min(as.numeric(tab_comb[memb,5]))
		}
		if(max(as.numeric(tab_comb[memb,6]))<pos[4]){
			pos[4]<-max(as.numeric(tab_comb[memb,6]))
		}
		for(j in 1:length(memb)){
			tmp<-which(sites1[[memb[j]]]==sites[i])
			sites1[[memb[j]]][tmp]<-paste(pos,collapse="_")
		}
		
	}
	sites1
}	

# function to modify the borders of conserved interaction sites based interactions assigned to the site
# uses a consensus to define borders		
refine_consensus_site<-function(out_comb, tab_comb,posprobs=posprobs, thres_m=0.2, thres_s=0.5){
	sites1<-out_comb
	sites<-na.omit(unique(unlist(sites1)))
	na<-grep("NA",sites)
	if(length(na)>0){
		sites<-sites[-na]
	}
	for(i in 1:length(sites)){
		pos_new<-c()
		memb<-grep(sites[i], sites1)
		pos<-as.numeric(strsplit(sites[i],"_")[[1]])
		org<-tab_comb[memb,"orgs"]
		a_table<-posprobs[[org[1]]]
		for(j in 2:length(memb)){
			a_table<-a_table+posprobs[[org[j]]]
		}
		tmp<-a_table[pos[1]:pos[2],pos[3]:pos[4]]
		tmp<-tmp/max(tmp)
		tmp2<-rowSums(tmp)
		tmp2<-tmp2/max(tmp2)
		pos2<-which(tmp2>=thres_m)+pos[1]-1
		pos_new[1]<-min(pos2)
		pos_new[2]<-max(pos2)
		tmp2<-colSums(tmp)
		tmp2<-tmp2/max(tmp2)
		pos2<-which(tmp2>=thres_s)+pos[3]-1
		pos_new[3]<-min(pos2)
		pos_new[4]<-max(pos2)
		for(j in 1:length(memb)){
			tmp<-which(sites1[[memb[j]]]==sites[i])
			sites1[[memb[j]]][tmp]<-paste(pos_new,collapse="_")
		}
	}
	sites1
}	
		
# condense peaks that largely overlap into one peak
condense_peak<-function(peaks, string=F, ov=10){
	if(string==T){
		tmp<-strsplit(peaks,"_")
		peaks<-(lapply(tmp, function(x){return(as.numeric(x[1:4]))}))
	}
	dup<-which(duplicated(unlist(lapply(peaks,paste, collapse="_"))))
	if(length(dup)>0){
		peaks<-peaks[-dup]
	}
	peaks_l<-order(unlist(lapply(peaks, function(x){return(x[2]-x[1])})), decreasing=T)
	peaks<-peaks[peaks_l]
	peaks2<-peaks
	i<-1
	while(i<=(length(peaks2)-1)){
		m_i<-peaks2[[i]][1]:peaks2[[i]][2]
		s_i<-peaks2[[i]][3]:peaks2[[i]][4]
		tmp<-c()
		for(j in (i+1):length(peaks2)){
			m_i2<-peaks2[[j]][1]:peaks2[[j]][2]
			s_i2<-peaks2[[j]][3]:peaks2[[j]][4]
			i_m<-max(length(intersect(m_i,m_i2)),0,na.rm=T)/min(length(m_i),length(m_i2))
			i_s<-max(length(intersect(s_i,s_i2)),0,na.rm=T)/min(length(s_i),length(s_i2))
			i_m2<-max(length(m_i),length(m_i2))-max(length(intersect(m_i,m_i2)),0,na.rm=T)
			i_s2<-max(length(s_i),length(s_i2))-max(length(intersect(s_i,s_i2)),0,na.rm=T)
			if(min(i_m,i_s)>0.8 & i_m2 < ov){
				tmp<-c(tmp, paste(peaks2[[j]],collapse="_"))
			}
		}
		m<-match(unique(tmp),unlist(lapply(peaks2, paste,collapse="_")))
		if(length(m)>0){
			peaks2<-peaks2[-m]
		}
		i<-i+1
	}
	if(string==T){
		peaks2<-lapply(peaks2, paste, collapse="_")
	}
	peaks2
}
	
# combine overlapping peaks to master peak and sub peak	
condense_peak2<-function(test){
	peaks<-test[[3]]
	tmp<-strsplit(peaks,"_")
	peaks<-(lapply(tmp, function(x){return(as.numeric(x[1:4]))}))
	peaks3<-list()
	com<-rbind(test[[1]],test[[2]])
	peaks_l<-order(unlist(lapply(peaks, function(x){return(x[2]-x[1])})), decreasing=T)
	peaks<-peaks[peaks_l]
	peaks2<-peaks
	i<-1
	while(i<=(length(peaks2))-1){
		m_i<-peaks2[[i]][1]:peaks2[[i]][2]
		s_i<-peaks2[[i]][3]:peaks2[[i]][4]
		int<-grep(paste(peaks2[[i]],collapse="_"),com[,"cluster_id"])
		tmp<-c()
		peaks2[[i]]<-paste(peaks2[[i]],collapse="_")
		for(j in (i+1):length(peaks2)){
			m_i2<-peaks2[[j]][1]:peaks2[[j]][2]
			s_i2<-peaks2[[j]][3]:peaks2[[j]][4]
			i_m<-max(length(intersect(m_i,m_i2)),0,na.rm=T)/min(length(m_i),length(m_i2))
			i_s<-max(length(intersect(s_i,s_i2)),0,na.rm=T)/min(length(s_i),length(s_i2))
			
			i_m2<-max(length(m_i),length(m_i2))-max(length(intersect(m_i,m_i2)),0,na.rm=T)
			i_s2<-max(length(s_i),length(s_i2))-max(length(intersect(s_i,s_i2)),0,na.rm=T)
			int2<-grep(paste(peaks2[[j]],collapse="_"),com[,"cluster_id"])
			int2<-length(intersect(int,int2))
			if(min(i_m,i_s)>0.65 & int2 >0){
				peaks2[[i]]<-c(peaks2[[i]],paste(peaks2[[j]],collapse="_"))
				tmp<-c(tmp, paste(peaks2[[j]],collapse="_"))
			}
		}
		m<-match(unique(tmp),unlist(lapply(peaks2, paste,collapse="_")))
		if(length(m)>0){
			peaks2<-peaks2[-m]
		}
		i<-i+1
	}
	if(is.numeric(peaks2[[length(peaks2)]][1])){
		peaks2[[length(peaks2)]]<-paste(peaks2[[length(peaks2)]],collapse="_")
	}
	peaks2
}


check_length<-function(summ2, min_length=6){
	tmp<-lapply(summ2, function(x){return(as.numeric(strsplit(x,"_")[[1]]))})
	tmp<-unlist(lapply(tmp, function(x){return(x[2]-x[1])}))
	out<-summ2[which(tmp>=min_length)]
	out
}

# function to identify peaks and assign interactions to peaks
peak_find2D<-function(align_table, tab_aligned, tabsub_aligned, alignment,sRNA_alignment2, min_thres=2, perc=0.1, min_length=6,eps=0.4,ooi=ooi, ooi_force=TRUE, posprobs=posprobs)	{
	fun3<-function(x, se="|"){
		if(is.null(x[1])==FALSE & is.na(x[1])==F){
			y<-strsplit(x,"_")
			y<-lapply(y, as.numeric)
			y<-unlist(lapply(y,function(x){x[2]-x[1]}))
			x<-paste(x[order(y,decreasing=T)], collapse="|")
		}
		else{
			x<-NA
		}
		x
	}
	o<-nrow(tab_aligned)
	thres<-max(2, nrow(tab_aligned)*0.1)
	peaks<-c()
	summ2<-c()
	tab_aligned2<-tab_aligned
	tabsub_aligned2<-tabsub_aligned
	tab_aligned2[,"name"]<-paste(tab_aligned[,"name"],"opt",sep="_")
	
	tabsub_aligned2[,"name"]<-paste(tabsub_aligned[,"name"],"sub",sep="_")
	con<-contourLines(x=1:nrow(align_table),y=1:ncol(align_table),z=align_table)#, levels=perc)
	if(length(con)>0){
		peaks<-list()
		for(i in 1:length(con)){
			m<-c(min(con[[i]][[2]]),max(con[[i]][[2]]))
			s<-c(min(con[[i]][[3]]),max(con[[i]][[3]]))
			ms<-round(c(m,s))
			if((ms[2]-ms[1]+1)>=min_length & (ms[4]-ms[3]+1)>=min_length ){
				peaks[[(length(peaks)+1)]]<-round(c(m,s))
			}
		}
		if(length(peaks)>0){
			peaks<-condense_peak(peaks,ov=20)
		
			out_opt<-assign_int3(peaks,tab_aligned, alignment,sRNA_alignment2, ov_thres=0.3 )
			out_sub<-list()
			if(nrow(tabsub_aligned)>0){
				out_sub<-assign_int3(peaks,tabsub_aligned, alignment,sRNA_alignment2, ov_thres=0.3 )
			}
			out<-c(unlist(out_opt),unlist(out_sub))
			summ<-table(out)
			summ2<-names(summ)[which(summ>=thres)]
			if(length(summ2)>0){
				peaks<-lapply(summ2,function(x){return(as.numeric(strsplit(x,"_")[[1]]))})
				out_opt<-assign_int3(peaks,tab_aligned, alignment,sRNA_alignment2, ov_thres=0.3 )
				out_sub<-list()
				if(nrow(tabsub_aligned)>0){
					out_sub<-assign_int3(peaks,tabsub_aligned, alignment,sRNA_alignment2, ov_thres=0.3 )
				}
				out<-c(unlist(out_opt),unlist(out_sub))
				summ<-table(out)
				summ2<-names(summ)[which(summ>=thres)]
				if(length(summ2)>0){				
					out_opt<-lapply(out_opt, function(x){
							tmp<-na.omit(match(summ2,x ))
							if(length(tmp)>0){
								return (x[sort(tmp)])
							} else{
								return (NA)
							}
							})
					out_sub<-lapply(out_sub, function(x){
							tmp<-na.omit(match(summ2,x ))
							if(length(tmp)>0){
								return (x[sort(tmp)])
							} else{
								return (NA)
							}
							})
					names(out_opt)<-rep("opt",length(out_opt))
					if(nrow(tabsub_aligned)>0){
						names(out_sub)<-rep("sub",length(out_sub))
					}
					out_comb<-c(out_opt,out_sub)
					tab_comb<-rbind(tab_aligned2,tabsub_aligned2)
					if(length(out_comb)>0){
						out_comb<-refine_site(out_comb, tab_comb)
					}
					summ2<-na.omit(unique(unlist(out_comb)))
					summ2<-check_length(summ2, min_length=min_length)
					out_comb_old<-1
					count<-1
					while(identical(out_comb_old,out_comb)==F & length(summ2)>0){
						out_comb_old<-out_comb
						out_comb<-comb_peaks3(tab_comb,summ2,alignment,sRNA_alignment2,out_comb, minpts=2, eps=eps)
						nu<-which(unlist(lapply(out_comb,is.null)))
						if(length(nu)<length(out_comb)){
							out_comb[nu]<-NA
							out_comb<-refine_consensus_site(out_comb, tab_comb, posprobs=posprobs, thres_m=0.2, thres_s=0.5)
						}	
						summ2<-table(unlist(out_comb))
						summ2<-names(summ2)[which(summ2>=thres)]
						summ2<-check_length(summ2, min_length=min_length)
						if(length(summ2)>0){
							out_comb<-lapply(out_comb, function(x){
							tmp<-na.omit(match(summ2,x ))
							if(length(tmp)>0){
								return (x[sort(tmp)])
							} else{
								return (NA)
							}
							})
							# check if sites are unique, otherwise condense them
							for(i in 1:length(out_comb)){
								if(is.na(out_comb[i])[1]==F){
									out_comb[[i]]<-unlist(condense_peak(out_comb[[i]], string=T, ov=min_length))
								}
							}
							summ2<-table(unlist(out_comb))
							summ2<-names(summ2)
							summ2<-check_length(summ2, min_length=min_length)
							if(length(summ2)>0){
								out_comb<-assign_int_fin(summ2, out_comb,tab_comb, alignment,sRNA_alignment2, ov_thres=0.60, ov_thres2=0.45)
								summ2<-names(table(unlist(out_comb)))
								summ2<-check_length(summ2, min_length=min_length)
							}
						}
					}
					if(length(summ2)>0){
						cluster_id<-lapply(out_comb,fun3)
						tab_comb<-cbind(tab_comb,cluster_id)
						tab_aligned<-tab_comb[1:o,]
						cluster_id<-c()
						tabsub_aligned<-cbind(tabsub_aligned,cluster_id)
						if(nrow(tabsub_aligned)>0){
							tabsub_aligned<-tab_comb[((o+1):nrow(tab_comb)),]
							if(is.matrix(tabsub_aligned)==F){
								tabsub_aligned<-t(as.matrix(tabsub_aligned))
							}
						}
					}
				}
			}
		}
	}
	if(length(con)==0 | length(peaks)==0 | length(summ2)==0){
		cluster_id<-rep(NA,nrow(tab_aligned))
		tab_aligned<-cbind(tab_aligned,cluster_id)
		cluster_id<-c()
		if(nrow(tabsub_aligned)>0){
			cluster_id<-rep(NA,nrow(tabsub_aligned))
		}
		tabsub_aligned<-cbind(tabsub_aligned,cluster_id)
	}
	if(ooi_force==TRUE){
		for(jjj in 1:length(ooi)){
			ooipos<-grep(ooi[jjj], tab_aligned[,"orgs"])
			if(length(ooipos)==1){
				if(is.na(tab_aligned[ooipos, "cluster_id"])){
					peaks<-list(c(as.numeric(tab_aligned[ooipos,1]),as.numeric(tab_aligned[ooipos,2]),as.numeric(tab_aligned[ooipos,5]),as.numeric(tab_aligned[ooipos,6])))
					#peaks<-list(paste(peaks, collapse="_"))
					out_opt<-assign_int3(peaks,tab_aligned, alignment,sRNA_alignment2, ov_thres=0.3 )
					out_sub<-list()
					if(nrow(tabsub_aligned)>0){
						out_sub<-assign_int3(peaks,tabsub_aligned, alignment,sRNA_alignment2, ov_thres=0.3 )
					}
					
					out_opt<-lapply(out_opt,fun3)
					out_sub<-lapply(out_sub,fun3)
					out<-c(unlist(out_opt),unlist(out_sub))
					summ<-table(out)
					summ2<-names(summ)
					out_opt<-unlist(out_opt)
					out_sub<-unlist(out_sub)
					names(out_opt)<-rep("opt",length(out_opt))
					if(nrow(tabsub_aligned)>0){
						names(out_sub)<-rep("sub",length(out_sub))
					}
					out_comb<-c(out_opt,out_sub)
					tab_comb<-rbind(tab_aligned2,tabsub_aligned2)
					if(length(na.omit(out_comb))>1){
						summ2<-out_comb[ooipos]
						out_comb<-comb_peaks3(tab_comb,summ2,alignment,sRNA_alignment2,out_comb,minpts=1)
						pos<-which(out_comb==unlist(out_comb[ooipos]))
						if(length(pos)>0){
							for(jj in 1:length(pos)){
								if(pos[jj]<=o){
									empt<-which(is.na(tab_aligned[pos[jj],"cluster_id"]))
									tab_aligned[pos[jj][empt],"cluster_id"]<-summ2
								}
								if(pos[jj]>o){
									empt<-which(is.na(tabsub_aligned[pos[jj]-o,"cluster_id"]))
									tabsub_aligned[pos[jj][empt]-o,"cluster_id"]<-summ2
								}
							}
						}
					}
					if(length(na.omit(out_comb))==1){
						tab_aligned[ooipos,"cluster_id"]<-summ2
					}
				}
			}
			
		}
	}
	out<-list(tab_aligned,tabsub_aligned,summ2)
	out
}

# compare interaction sites that are assigned to a peak in the site probability landscape
comb_peaks3<-function(tab_comb,summ2,alignment,sRNA_alignment2,out_comb, minpts=2, eps=0.55){
	out<-vector("list", nrow(tab_comb))
	peak_consensus<-vector("list",length(summ2))
	# investigate each peak
	for(ii in 1:length(summ2)){
		# which interaction sites are assigned to peak ii
		cands<-grep(summ2[ii], out_comb)
		
		m<-matrix(,length(cands),length(cands))
		m[]<-0
		colnames(m)<-tab_comb[cands,"name"]
		rownames(m)<-tab_comb[cands,"name"]
		m_id<-m

		# cut the peak position from the alignment and re-align the sequences
		# with the consistency based aligner dialign-tx
		coo1<-as.numeric(strsplit(summ2[[ii]],"_")[[1]][1:4])
		fasta_temp<-alignment[as.character(unique(tab_comb[cands,"name2"]))]
		s=max(1,coo1[1]-2)
		e=min(coo1[2]+2,length(fasta_temp[[1]]))
		fasta_temp<-lapply(fasta_temp, cutalign, s=s, e=e)
		fasta_temp2<-lapply(fasta_temp, remove_gaps)
		short<-which(unlist(lapply(fasta_temp2,length))<2)
		fasta_s<-sRNA_alignment2[as.character(unique(tab_comb[cands,"orgs"]))]
		if(length(short)>0){
			fasta_temp<-fasta_temp[-short]
			fasta_temp2<-fasta_temp2[-short]
			short_cands<-match(names(fasta_temp2)[short], tab_comb[cands,"name"])
			cands<-cands[-short_cands]
			fasta_s<-fasta_s[-short]
		}
		if(length(fasta_temp2)>1){
			tempfi<-tempfile(pattern="CopraRNA2.findConservedSites.")
			tempfi2<-tempfile(pattern="CopraRNA2.findConservedSites.")
			write.fasta(fasta_temp2, names=tab_comb[cands,"name"], file.out=tempfi)
			# avoid massive information output of dialign-tx (> /dev/null)
			command<-paste("dialign-tx -D ", dialign_conf," ",tempfi, " ", tempfi2, " > /dev/null 2>> CopraRNA2_subprocess.oe")
			system(command)
			fasta_temp2<-read.fasta(tempfi2)
			s2=max(1,coo1[3]-2)
			e2=min(coo1[4]+2,length(fasta_s[[1]]))
			fasta_s<-lapply(fasta_s, cutalign, s=s2, e=e2)
			fasta_s2<-lapply(fasta_s, remove_gaps)
			short2<-which(unlist(lapply(fasta_s2,length))<2)
			if(length(short2)>0){
				fasta_s<-fasta_s[-short2]
				fasta_s2<-fasta_s2[-short2]
				short_cands2<-match(names(fasta_s2)[short], tab_comb[cands,"name"])
				cands<-cands[-short_cands2]
				fasta_temp<-fasta_temp[-short2]
				fasta_temp2<-fasta_temp2[-short2]
			}
			if(length(fasta_s2)>1){
				write.fasta(fasta_s2, names=tab_comb[cands,"name"], file.out=tempfi)
				# avoid massive information output of dialign-tx (> /dev/null)
				command<-paste("dialign-tx -D ", dialign_conf," ",tempfi, " ", tempfi2, " > /dev/null 2>> CopraRNA2_subprocess.oe")
				system(command)
				fasta_s2<-read.fasta(tempfi2)
				unlink(tempfi)
				unlink(tempfi2)
				
				# after re-alignment, calculate pairwise positional overlaps and identities of the interactions assigned to peak ii
				for(i in 1:length(cands)){
					a_i<-match(tab_comb[cands[i],"name2"],names(fasta_temp))
					for(j in i:length(cands)){
						if(i!=j){
							a_j<-match(tab_comb[cands[j],"name2"],names(fasta_temp))
							a<-as.numeric(strsplit(tab_comb[cands[i],"int_mrna"],",")[[1]])
							a<-a-s+1
							as<-as.numeric(strsplit(tab_comb[cands[i],"int_srna"],",")[[1]])
							as<-as-s2+1
							ga<-which(fasta_temp[[a_i]]!="-")
							a<-match(a,ga)
							ga2<-which(fasta_temp2[[a_i]]!="-")
							a<-ga2[a]
							gas<-which(fasta_s[[a_i]]!="-")
							as<-match(as,gas)
							gas2<-which(fasta_s2[[a_i]]!="-")
							as<-gas2[as]
							as<-rev(as)
							temp<-paste(a,as)
							na<-grep("NA",temp)
							if(length(na)>0){
								temp<-temp[-na]
							}
							b<-a
							a<-as.numeric(strsplit(tab_comb[cands[j],"int_mrna"],",")[[1]])
							a<-a-s+1
							as<-as.numeric(strsplit(tab_comb[cands[j],"int_srna"],",")[[1]])
							as<-as-s2+1
							ga<-which(fasta_temp[[a_j]]!="-")
							a<-match(a,ga)
							ga2<-which(fasta_temp2[[a_j]]!="-")
							a<-ga2[a]
							gas<-which(fasta_s[[a_j]]!="-")
							as<-match(as,gas)
							gas2<-which(fasta_s2[[a_j]]!="-")
							as<-gas2[as]
							as<-rev(as)
							temps<-paste(a,as)
							na<-grep("NA",temps)
							if(length(na)>0){
								temps<-temps[-na]
							}
							ov<-length(intersect(temp,temps))/length(union(temp,temps))
							m[i,j]=m[j,i]=ov
							mi<-min(c(a,b),na.rm = T)
							ma<-max(c(a,b),na.rm = T)
							if(is.finite(mi) & is.finite(ma)){
								mima<-mi:ma
								gap<-which(fasta_temp2[[a_i]][mima]=="-" | fasta_temp2[[a_j]][mima]=="-")
								if(length(gap)>0){
									mima<-mima[-gap]
								}
								ov2<-which(fasta_temp2[[a_i]][mima]==fasta_temp2[[a_j]][mima])
								ov2<-length(ov2)/length(mima)
								m_id[i,j]=m_id[j,i]=ov2
							}
							if(is.infinite(mi) | is.infinite(ma)){
								m_id[i,j]=m_id[j,i]=0
								m[i,j]=m[j,i]=0
							}
						}
					}
				}
				m<-(m*m_id)**0.5  # geometric mean between positional overlap and sequence identity
				m<-1-m # similarity to distance
				diag(m)<-0
				clus<-cluster(m, eps=eps,minpts=minpts)
				if(length(clus)>0){
					for(jj in 1:length(clus)){
						na<-paste(summ2[[ii]],jj,sep="_")
						pos<-match(clus[[jj]],tab_comb[,"name"])
						for(jjj in 1:length(pos)){
							out[[pos[jjj]]]<-c(out[[pos[jjj]]],na)
						
						}
						
						
						
					}
				}
			}
		}
	}
	out
}

# pre-assign interaction sites to peaks defined by the site probability landscape
assign_int3<-function(peaks, tab_aligned, alignment,sRNA_alignment2, ov_thres=0.3){

	colid_m<-"name2"
	colid_s<-"orgs"
	out_opt<-rep(list(NA),nrow(tab_aligned))
	# investigate each mRNA sRNA interaction pair
	for(i in 1:nrow(tab_aligned)){

		# readout start and end of the predicted mRNA site mapped to the alignment
		# and build a vector from start:end
		m<-seq(as.numeric(tab_aligned[i,1]),as.numeric(tab_aligned[i,2]))
		# readout start and end of the predicted sRNA site mapped to the alignment
		# and build a vector from start:end
		s<-seq(as.numeric(tab_aligned[i,5]),as.numeric(tab_aligned[i,6]))
		
		# take the length of the mRNA and sRNA interaction vectors and substract the gap positions in the alignment
		le_m<-length(m)-length(which(alignment[[as.character(tab_aligned[i,colid_m])]][as.numeric(tab_aligned[i,1]):as.numeric(tab_aligned[i,2])]=="-"))
		le_s<-length(s)-length(which(sRNA_alignment2[[as.character(tab_aligned[i,colid_s])]][as.numeric(tab_aligned[i,5]):as.numeric(tab_aligned[i,6])]=="-"))
		
	
		
		# compare the overlap of interaction pair i with all detected peaks in the site probability landscape
		
		temp<-rep(0, length(peaks))		# stores ratio of overlap length and length of union from interaction site and peak for the mRNA
		temp_l<-rep(0, length(peaks))	# stores length of overlap mRNA/peak
		temp_s<-rep(0, length(peaks))	# stores ratio of overlap length and length of union from interaction site and peak for the sRNA
		temp_ls<-rep(0, length(peaks))	# stores length of overlap sRNA/peak
		
		for(j in 1:length(peaks)){
			# built vectors from start to end of x (mRNA) and y (sRNA) position of the peaks 
			m_p<-peaks[[j]][1]:peaks[[j]][2]
			s_p<-peaks[[j]][3]:peaks[[j]][4]
			
			# substract the gaps from the vector length
			le_m_p<-length(m_p)-length(which(alignment[[as.character(tab_aligned[i,colid_m])]][m_p]=="-"))
			le_s_p<-length(s_p)-length(which(sRNA_alignment2[[as.character(tab_aligned[i,colid_s])]][s_p]=="-"))
			
			# calculate length of overlap between predicted interaction vectors and peak vectors without gaps
			inter_m<-length(intersect(m, m_p))-length(which(alignment[[as.character(tab_aligned[i,colid_m])]][intersect(m, m_p)]=="-"))
			inter_s<-length(intersect(s, s_p))-length(which(sRNA_alignment2[[as.character(tab_aligned[i,colid_s])]][intersect(s, s_p)]=="-"))
			
			# if and overlap between peaks exist for mRNA AND sRNA
			# and the each overlap is bigger than the "ov_thres" threshold for the interaction site OR the peak
			# calculate percentage 
			if(inter_m>0 & inter_s>0){
				if(inter_m/ le_m >= ov_thres | inter_m/ le_m_p >= ov_thres){
					if(inter_s/ le_s >= ov_thres | inter_s/ le_s_p >= ov_thres){
						temp[j]<-length(intersect(m, m_p))/length(union(m, m_p))  # Ratio of overlap length and length of union from interaction site and peak for the mRNA
						temp_l[j]<-inter_m # length of overlap
						temp_s[j]<-length(intersect(s, s_p))/length(union(s, s_p))  # Ratio of overlap length and length of union from interaction site and peak for the sRNA
						temp_ls[j]<-inter_s # length of overlap
					}
				}
			}
			
		}
		
		# actual assignment of a predicted interaction to a peak
		if(max(temp)>0 & max(temp_s)>0){
				# normalize the ratio and overlap vectors and calculate the geometric mean
				temp2<-(temp/max(temp)*temp_l/max(temp_l)*temp_s/max(temp_s)*temp_ls/max(temp_ls))^0.25
				names(temp2)<-unlist(lapply(peaks,function(x){
					if(is.null(x)==FALSE & length(x)){
						x<-paste(x, collapse="_")
					}
					else{
						x<-NA
					}	
					x
					}))
				names(temp)<-names(temp2)
				temp<-temp[order(temp2,decreasing = T)]
				out_opt[[i]]<-names(temp)[which(temp>0)]#[1]
				
		}
	}
	out_opt
}


# finally interaction sites to peaks defined by the site probability landscape 
assign_int_fin<-function(peaks, out_comb, tab_comb, alignment,sRNA_alignment2, ov_thres=0.6, ov_thres2=0.5){
	tab_aligned<-tab_comb
	colid_m<-"name2"
	colid_s<-"orgs"
	out_opt<-rep(list(NA),nrow(tab_aligned))
	# investigate each mRNA sRNA interaction pair
	for(i in 1:nrow(tab_aligned)){
		if(is.na(out_comb[i])==F){
			# readout start and end of the predicted mRNA site mapped to the alignment
			# and build a vector from start:end
			m<-seq(as.numeric(tab_aligned[i,1]),as.numeric(tab_aligned[i,2]))
			# readout start and end of the predicted sRNA site mapped to the alignment
			# and build a vector from start:end
			s<-seq(as.numeric(tab_aligned[i,5]),as.numeric(tab_aligned[i,6]))
			
			# take the length of the mRNA and sRNA interaction vectors and substract the gap positions in the alignment
			le_m<-length(m)-length(which(alignment[[as.character(tab_aligned[i,colid_m])]][as.numeric(tab_aligned[i,1]):as.numeric(tab_aligned[i,2])]=="-"))
			le_s<-length(s)-length(which(sRNA_alignment2[[as.character(tab_aligned[i,colid_s])]][as.numeric(tab_aligned[i,5]):as.numeric(tab_aligned[i,6])]=="-"))
			
		
			
			# compare the overlap of interaction pair i with all detected peaks in the site probability landscape
			
			temp<-rep(0, length(peaks))		# stores ratio of overlap length and length of union from interaction site and peak for the mRNA
			#temp_l<-rep(0, length(peaks))	# stores length of overlap mRNA/peak
			temp_s<-rep(0, length(peaks))	# stores ratio of overlap length and length of union from interaction site and peak for the sRNA
			#temp_ls<-rep(0, length(peaks))	# stores length of overlap sRNA/peak
			
			peaks<-lapply(out_comb[[i]], function(x){return(as.numeric(strsplit(x,"_")[[1]]))})
			
			for(j in 1:length(peaks)){
				# built vectors from start to end of x (mRNA) and y (sRNA) position of the peaks 
				m_p<-peaks[[j]][1]:peaks[[j]][2]
				s_p<-peaks[[j]][3]:peaks[[j]][4]
				
				# substract the gaps from the vector length
				le_m_p<-length(m_p)-length(which(alignment[[as.character(tab_aligned[i,colid_m])]][m_p]=="-"))
				le_s_p<-length(s_p)-length(which(sRNA_alignment2[[as.character(tab_aligned[i,colid_s])]][s_p]=="-"))
				
				# calculate length of overlap between predicted interaction vectors and peak vectors without gaps
				inter_m<-length(intersect(m, m_p))-length(which(alignment[[as.character(tab_aligned[i,colid_m])]][intersect(m, m_p)]=="-"))
				inter_s<-length(intersect(s, s_p))-length(which(sRNA_alignment2[[as.character(tab_aligned[i,colid_s])]][intersect(s, s_p)]=="-"))
				
				# if and overlap between peaks exist for mRNA AND sRNA
				# and the each overlap is bigger than the "ov_thres" threshold for the peak
				# calculate percentage 
				temp[j]<-length(intersect(m, m_p))/length(m_p)  # Ratio of overlap length and length of the peak for the mRNA
				temp_s[j]<-length(intersect(s, s_p))/length(s_p)  # Ratio of overlap length and length of thed peak for the sRNA
		
			}
		
			# actual assignment of a predicted interaction to a peak
			if(max(temp)>0 & max(temp_s)>0){
					# normalize the ratio and overlap vectors and calculate the geometric mean
					temp2<-(temp*temp_s)^0.5
					names(temp2)<-unlist(lapply(peaks,function(x){
								x<-paste(x, collapse="_")
								x
						}))
					names(temp)<-names(temp2)
					temp<-temp[order(temp2,decreasing = T)]
					tmp3<-names(temp)[which(temp>ov_thres)]
					if(length(tmp3)>0){
						out_opt[[i]]<-tmp3
					} else {
						tmp3<-names(temp)[which(temp>ov_thres2)][1]
						if(length(tmp3)>0){
							out_opt[[i]]<-tmp3
						}
					}
					
					
			}
		}
	}
	out_opt
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

#mafft("16s_sequences.fa", outname="16s_sequences_mafft_align.fas", mode="accurate")
ribo<-read.phyDat("16s_sequences_mafft_align.fas", format="fasta", type="DNA")
dm <- dist.ml(ribo, model="F81")
tree_rib<- NJ(dm)
weight1<-read.csv("weights.txt",skip=1, header=F, sep="\t")
weight<-weight1[,2]
names(weight)<-weight1[,1]


# function to map interaction positions to the alignments and wrapper for the other functions
build_anno<-function(ooi="NC_000911", conservation_oois=ooi){
	dir.create("evo_alignments2")
	fastutr<-read.fasta("utr_seqs.fa") # for own alignments
	fastnames<-tolower(names(fastutr)) # for own alignments
	
	# create alignment of sRNAs
	sRNA_alignment<-parse_fasta(system("mafft --maxiterate 1000 --localpair --quiet input_sRNA.fa ", intern=T))

	wd<-getwd()
	
	# reference for homolog clusters
	dat<-read.csv("CopraRNA2_prep_anno_addhomologs_padj_amountsamp.csv")
	
	dat<-as.matrix(dat)
	e<-grep("Annotation", colnames(dat))
	all_orgs<-colnames(dat)[3:(e-1)]
	
	# read optimal IntaRNA results from all organisms
	d2<-dir()
	f<-grep("tags.clustered$", d2)
	opt<-read.csv(d2[f], sep=";")
	subopt1<-grep("_subopt.intarna.csv",d2)
	subopt<-c()
	for(i in 1:length(subopt1)){
		temp<-read.csv(d2[subopt1[i]], sep=";")
		subopt<-rbind(subopt,temp)
	}
	subopt[,1]<-tolower(subopt[,1])
	opt<-as.matrix(opt)
	subopt<-as.matrix(subopt)
	
	# read availbale organism file for annotation purposes
	copref<-read.delim(copref_path, sep="\t", header=T,comment.char = "#")	
	e<-grep("Annotation", colnames(dat))
	temp_dat<-dat[,3:(e-1)]
	s_names<-colnames(dat)[3:(e-1)]
	names(sRNA_alignment)<-s_names
	nam<-c()
	for(i in 1:ncol(temp_dat)){
		tnam<-grep(gsub("\\..*","",colnames(temp_dat)[i]),copref[,1])
		nam<-c(nam,as.character(copref[tnam,2]))
	}
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
	tree_rib[[3]]<-paste(nam2,tree_rib[[3]])
	
	
	#select top predictions for interaction site analysis
	copra<-read.csv(copra_result, sep=",")
	ooi_pos<-grep(ooi, colnames(copra))
	ooi_pos2<-grep(ooi, colnames(dat))
	copra<-copra[1:top,ooi_pos]
	copra<-gsub("\\(.*","",copra)
	ex<-na.omit(match(copra, gsub("\\(.*","",dat[,ooi_pos2])))
	align_pos<-1:nrow(dat)
	dat_old<-dat
	dat<-dat[ex,]
	align_pos<-ex
	dat_all<-dat
	
	# divide the data in subsets for parallel processing 
	max_cores<-min(top, max_cores)
	jobs<-nrow(dat)%/%max_cores
	rest<-nrow(dat)-max_cores*jobs
	jobs<-rep(jobs,max_cores)
	if(rest>0){
		jobs[1:rest]<-jobs[1:rest]+1
	}
	count_vect1<-cumsum(c(1,jobs[1:(length(jobs)-1)]))
	count_vect2<-cumsum(jobs)

	# generate a temp file per thread to prepare mafft input file
	thread2tmpfile = c();
	for (i in 1:max_cores) { thread2tmpfile = c(thread2tmpfile, tempfile(pattern="CopraRNA2.findConservedSites.")); }

	# start parallel processing of site conservation
	vari<-foreach(ji=1:max_cores)  %dopar% {
		dat<-dat_all[count_vect1[ji]:count_vect2[ji],]
		if(is.matrix(dat)==F){
			dat<-t(as.matrix(dat))
		}
		all_orgs<-colnames(dat)[3:(e-1)]
		int_sites<-vector("list", nrow(dat))
		peak_list<-vector("list", nrow(dat))
		for(i in  1:nrow(dat)){
			e<-grep("Annotation", colnames(dat))
			temp<-dat[i,3:(e-1)]
			temp<-temp[which(temp!="")]
			orgs<-names(temp)
			if(length(temp)>1){
				poslist<-c()
				query<-gsub("\\(.*","",as.character(temp))
				genename<-gsub(".*\\(","",as.character(temp))
				genename<-gsub("\\|.*", "", genename)
				genename<-gsub("\\/", "", genename)
		
				# read information from the optimal and suboptimal IntaRNA predictions for all UTRs of the respective homolog cluster
				pos_opt<-(match(query, opt[,1]))
				tab<-opt[pos_opt,c("start1","end1","seedStart1","seedEnd1","start2","end2","seedStart2","seedEnd2")]
				tab<-cbind(tab,paste(genename, "_", query,sep=""),query,opt[pos_opt,"p.value"],orgs,opt[pos_opt,"hybridDP"],opt[pos_opt,"E"],opt[pos_opt,"seq1"],opt[pos_opt,"seq2"])
				colnames(tab)<-c("start","end","start_seed", "end_seed","start_sRNA","end_sRNA","start_seed_sRNA", "end_seed_sRNA","name" ,"name2","pvalue","orgs","hybridDP","Energy","seq_mrna","seq_srna")
				
				pos_sub<-(match(query, (subopt[,1])))
				exist_sub<-which(is.na(pos_sub)==F)
				pos_sub<-na.omit(pos_sub)
				tabsub<-matrix(,0,16)
				if(length(pos_sub)>0){
					tabsub<-subopt[pos_sub,c("start1","end1","seedStart1","seedEnd1","start2","end2","seedStart2","seedEnd2")]
					if(is.matrix(tabsub)==F){
						tabsub<-t(as.matrix(tabsub))
					}
					tabsub<-cbind(tabsub,paste(genename[exist_sub], "_", query[exist_sub],sep=""),query[exist_sub],subopt[pos_sub,"p.value"],orgs[exist_sub],subopt[pos_sub,"hybridDP"],subopt[pos_sub,"E"],subopt[pos_sub,"seq1"],subopt[pos_sub,"seq2"])
					
				}
				if(nrow(tabsub)==0){
					tabsub<-matrix(,1,16)
				}
				colnames(tabsub)<-c("start","end","start_seed", "end_seed","start_sRNA","end_sRNA","start_seed_sRNA", "end_seed_sRNA","name" ,"name2","pvalue","orgs","hybridDP","Energy","seq_mrna","seq_srna")
				tabsub<-na.omit(tabsub)
				pos<-as.numeric(names(sort(table(poslist), decreasing=T))[1])
							
				# create alignments 
				tempf2<-thread2tmpfile[ji];
				temp2<-na.omit(match(tolower(tab[,10]), fastnames))
				write.fasta(fastutr[temp2],file=tempf2 ,names=tab[,"name2"])
				ca<-paste("mafft --maxiterate 1000 --localpair --quiet --inputorder ", tempf2,"",sep="")
				alignment<-parse_fasta(system(ca, intern=T))
				m_names<-na.omit( match(orgs, names(dat[i,3:(e-1)])))
				s_names2<-na.omit(match(orgs, s_names))
				sRNA_alignment2<-sRNA_alignment[s_names2]
			
				ta<-exclude_nonbound(tab)
				tab<-cbind(tab,ta[[1]],ta[[2]],ta[[1]],ta[[2]])
				colnames(tab)[c(17,18,19,20)]<-c("int_mrna","int_srna","int_mrna2","int_srna2")
				tab_aligned<-as.matrix(tab)
				gaps_mRNA=gaps_sRNA=rep(NA,nrow(tab_aligned))
				tab_aligned<-cbind(tab_aligned,gaps_mRNA,gaps_sRNA)
				
				if(nrow(tabsub)>0){
					tasu<-exclude_nonbound(tabsub)
					tabsub<-cbind(tabsub,tasu[[1]],tasu[[2]],tasu[[1]],tasu[[2]])
				}
				if(nrow(tabsub)==0){
					tabsub<-cbind(tabsub,c(),c(),c(),c())
				}
				colnames(tabsub)[c(17,18,19,20)]<-c("int_mrna","int_srna","int_mrna2","int_srna2")
				tabsub_aligned<-as.matrix(tabsub)
				gaps_mRNA=gaps_sRNA=rep(NA,nrow(tabsub_aligned))
				tabsub_aligned<-cbind(tabsub_aligned,gaps_mRNA,gaps_sRNA)
				
				# map the position information of the IntaRNA predictions to the mRNA and sRNA alignments
				if(nrow(tabsub)>0){
					for(j in 1:nrow(tabsub)){
						align_num_sRNA<-which(names(sRNA_alignment2)==tabsub[j,"orgs"])
						align_num<-which(names(alignment)==tabsub[j,"name2"])[1]
						nogaps<-which(alignment[[align_num]]!="-")
						nogaps2<-which(sRNA_alignment2[[align_num_sRNA]]!="-")
						tabsub_aligned[j,1]<-nogaps[as.numeric(tabsub[j,1])]
						tabsub_aligned[j,2]<-nogaps[as.numeric(tabsub[j,2])]
						tabsub_aligned[j,5]<-nogaps2[as.numeric(tabsub[j,5])]
						tabsub_aligned[j,6]<-nogaps2[as.numeric(tabsub[j,6])]
						mpos<-as.numeric(strsplit(as.character(tabsub_aligned[j,"int_mrna"]),",")[[1]])
						tabsub_aligned[j,"int_mrna"]<-paste(nogaps[mpos],collapse=",")
						spos<-as.numeric(strsplit(as.character(tabsub_aligned[j,"int_srna"]),",")[[1]])
						tabsub_aligned[j,"int_srna"]<-paste(nogaps2[spos],collapse=",")
						}
				}
				 
				for(j in 1:nrow(tab)){
					align_num_sRNA<-which(names(sRNA_alignment2)==tab[j,"orgs"])
					align_num<-which(names(alignment)==tab[j,"name2"])[1]
					nogaps<-which(alignment[[align_num]]!="-")
					nogaps2<-which(sRNA_alignment2[[align_num_sRNA]]!="-")
					gaps_mRNA<-paste(nogaps,collapse=",")
					gaps_sRNA<-paste(nogaps2,collapse=",")
					tab_aligned[j,1]<-nogaps[as.numeric(tab[j,1])]
					tab_aligned[j,2]<-nogaps[as.numeric(tab[j,2])]
					tab_aligned[j,5]<-nogaps2[as.numeric(tab[j,5])]
					tab_aligned[j,6]<-nogaps2[as.numeric(tab[j,6])]
					mpos<-as.numeric(strsplit(as.character(tab_aligned[j,"int_mrna"]),",")[[1]])
					tab_aligned[j,"int_mrna"]<-paste(nogaps[mpos],collapse=",")
					spos<-as.numeric(strsplit(as.character(tab_aligned[j,"int_srna"]),",")[[1]])
					tab_aligned[j,"int_srna"]<-paste(nogaps2[spos],collapse=",")
					tab_aligned[j,"gaps_mRNA"]<-gaps_mRNA
					tab_aligned[j,"gaps_sRNA"]<-gaps_sRNA
				}
				
				
				# map individual spot probability matrices to the mRNA and sRNA alignments 
				# weight the spot probabilities based on the phylogenetic calculated organism weights
				# combine all spot probabilities and normalize to values between 0 and 1
				tab<-as.matrix(tab)
				align_table_int_pos<-matrix(,length(alignment[[1]]),length(sRNA_alignment[[1]]))
				align_table_int_pos[]<-0
				temp_align_table<-matrix(,length(alignment[[1]]),length(sRNA_alignment[[1]]))
				temp_align_table[]<-0
				
				
				
				posprobs<-vector("list",nrow(tab))
				names(posprobs)<-tab[,"orgs"]
				for(j in 1:nrow(tab_aligned)){
					temp_org<-tab_aligned[j,"orgs"]
					temp_locus<-tab_aligned[j,"name2"]
					
					# run IntaRNA for spot probabilities
					temp_table<-paste("IntaRNA  --target ",tab[j,"seq_mrna"] , " --tAccW " ,winsize, " --tAccL ",maxbpdist, " --query ",tab[j,"seq_srna"]," --qAccW ", winsize, " --qAccL ", maxbpdist, " --temperature ", temperature, " --out /dev/null --out=SpotProb:STDOUT", sep="")
					temp_table<-as.matrix(read.csv(textConnection(system(temp_table,intern=T)),sep=";", row.names=1,comment.char = "#"))
					temp_align_table2<-temp_align_table
					mRNA_no_gaps<-eval(parse( text=paste("c(",tab_aligned[j,"gaps_mRNA"],")",sep="") ))
					sRNA_no_gaps<-eval(parse( text=paste("c(",tab_aligned[j,"gaps_sRNA"],")",sep="") ))
					mRNA_gaps<-seq(1:nrow(temp_align_table2))[-mRNA_no_gaps]
					sRNA_gaps<-seq(1:ncol(temp_align_table2))[-sRNA_no_gaps]
					for(jjj in 1:length(mRNA_no_gaps)){
						temp_align_table2[mRNA_no_gaps[jjj],sRNA_no_gaps]<-temp_table[jjj,]
					}
					if(length(mRNA_gaps)>0){
						for(jjj in 1:length(mRNA_no_gaps)){
							gaps<-which(mRNA_gaps< mRNA_no_gaps[jjj] & mRNA_gaps> max(0,mRNA_no_gaps[jjj-1]))
							if(length(gaps)>0){
								for(jj in 1:ncol(temp_align_table2)){
									val<-min(temp_align_table2[mRNA_no_gaps[jjj],jj],temp_align_table2[mRNA_no_gaps[jjj-1],jj])
									temp_align_table2[mRNA_gaps[gaps],jj]<-val
								}
							}
						}
					}
					if(length(sRNA_gaps)>0){
						for(jjj in 1:length(sRNA_no_gaps)){
							gaps<-which(sRNA_gaps<sRNA_no_gaps[jjj] & sRNA_gaps> max(0,sRNA_no_gaps[jjj-1]))
							if(length(gaps)>0){
								for(jj in 1:nrow(temp_align_table2)){
									val<-min(temp_align_table2[jj,sRNA_no_gaps[jjj]],temp_align_table2[jj,sRNA_no_gaps[jjj-1]])
									
									temp_align_table2[jj,sRNA_gaps[gaps]]<-val
								}
							}
						}
					}
					posprobs[[j]]<-temp_align_table2*weight[temp_org]
					align_table_int_pos<-align_table_int_pos+temp_align_table2*weight[temp_org]
				}
				align_table<-align_table_int_pos
				align_table_int_pos<-align_table_int_pos/max(align_table_int_pos)
				
				# identify conserved interaction regions and map the organism specific predicted sites to these regions 
				test<-peak_find2D(align_table_int_pos, tab_aligned, tabsub_aligned,alignment,sRNA_alignment2, min_thres=0.2, perc=0.1, min_length=6,eps=0.55,ooi=conservation_oois, ooi_force=FALSE, posprobs=posprobs)
				na3<-paste("./evo_alignments2/",tab[1,"name"],sep="")
				dir.create(na3)
				
				#plot Spot Probabilities heatmap based on normalized multiple alignment
				x<-align_table
				w<-2000
				h<-ncol(x)/nrow(x)*w
				png(paste(na3,"/spot_probabilities.png",sep=""), width = w, height = h)	
					filled.contour(1:nrow(x),1:ncol(x),x,col=topo.colors(20), xlab="mRNA_alignment [nt]",ylab="sRNA_alignment [nt]",main=paste(tab[1,"name"]," combined interaction spot-probabilities",sep=""),zlim=c(0,1), 
						plot.axes = {contour(1:nrow(x),1:ncol(x),x, nlevels = 20, drawlabels = TRUE, axes = FALSE, rame.plot = FALSE, add = TRUE);axis(1); axis(2)}
					)
				dev.off()
								
				peaks1<-list()
				if(length(test[[3]])>0){
					# if a peak overlaps with another peak and both peaks are assigned to the same interaction site, the smaller peak is assigned a subpeak of the larger peak
					peaks1<-condense_peak2(test)
					ooi_peak<-test[[1]][1,"cluster_id"]
					ooi_peak<-gsub("\\|.*","",ooi_peak)
					if(is.na(ooi_peak)==F){
						ooi_p2<-grep(ooi_peak,peaks1)
						peaks1<-c(peaks1[ooi_p2],peaks1[-ooi_p2])
						
					}
					
				}
				tab_aligned<-as.matrix(test[[1]])
				tabsub_aligned<-as.matrix(test[[2]])
				if(is.matrix(tabsub_aligned)==F){
					tabsub_aligned<-t(as.matrix(tabsub_aligned))
				}
				if(is.matrix(tab_aligned)==F){
					tab_aligned<-t(as.matrix(tab_aligned))
				}
				tab_aligned[,"name"]<-gsub("_opt","",tab_aligned[,"name"])
				tabsub_aligned[,"name"]<-gsub("_sub","",tabsub_aligned[,"name"])
				jalview_anno(tab_aligned, tabsub_aligned,peaks1, test,ooi, nam2=nam2,alignment,sRNA_alignment2 ,name=paste(na3,"/",tab[1,"name"],sep=""),all_orgs=all_orgs)
				if(length(test[[3]])>0){
					op<-which(is.na(tab_aligned[,"cluster_id"])==F)
					su<-which(is.na(tabsub_aligned[,"cluster_id"])==F)
					sites1<-rep(NA, length(all_orgs))
					sites1<-rbind(sites1,sites1)
					colnames(sites1)<-all_orgs
					sites1[1,unlist(tab_aligned[op,"orgs"])]<-unlist(tab_aligned[op,"cluster_id"])
					sites1[2,unlist(tabsub_aligned[su,"orgs"])]<-unlist(tabsub_aligned[su,"cluster_id"])
					int_sites[[i]]<-sites1
					peak_list[[i]]<-peaks1
				}
			}
		}
		names(int_sites)<-dat[,ooi_pos2]
		names(peak_list)<-dat[,ooi_pos2]
		temp_out<-list(int_sites,peak_list)
		temp_out
	}

	# cleanup temp files
	file.remove(thread2tmpfile); 

	dat<-dat_old
	int_sites<-vector("list", nrow(dat))
	peak_list<-vector("list", nrow(dat))
	for(i in 1:length(vari)){
		tmp<-match(names(vari[[i]][[1]]),dat[,3])
		int_sites[tmp]<-vari[[i]][[1]]
		peak_list[tmp]<-vari[[i]][[2]]
		
	}
	save(int_sites, file="int_sites.Rdata")
	save(peak_list, file="peak_list.Rdata")
}
build_anno(ooi=ooi, conservation_oois=ooi)
