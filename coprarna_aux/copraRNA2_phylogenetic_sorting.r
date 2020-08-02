#!/usr/bin/env Rscript# script by Jens Georg

#call:
# R --slave -f ~/copraRNA2_phylogenetic_sorting.r 

suppressPackageStartupMessages(require(phangorn))
suppressPackageStartupMessages(require(seqinr))
suppressPackageStartupMessages(require(doMC))

# get absolute path
initial.options <- commandArgs(trailingOnly = FALSE)
path<-initial.options [4]
path<-sub("copraRNA2_phylogenetic_sorting.r","",path)
print(path)

# preset path to required files, path can also be specified as argument
copref_path<-paste(path,"CopraRNA_available_organisms.txt",sep="")


# read the organism of interest (ooi) from the ncRNA fasta file. The sRNA of the ooi is considered to be the first sequence.
ooi<-gsub("ncRNA_","",names(read.fasta("ncrna.fa"))[1])

# register cores for parallel processing
co<-readLines("CopraRNA_option_file.txt") 
max_cores<-as.numeric(gsub("core count=","",co[grep("core count=", co)]))
registerDoMC(max_cores)


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

# wrapper function for mafft
mafft<-function(filename="ncrna.fa", outname="ncrna_aligned.fa", mode="fast"){
	if(mode=="accurate"){
		command<-paste("mafft --maxiterate 1000 --localpair --quiet --inputorder ", filename, " > ", outname, sep="" )
	}
	if(mode=="fast"){
		command<-paste("mafft --retree 2 --maxiterate 0 --quiet --inputorder ", filename, " > ", outname, sep="" )
	}
	if(mode=="very_fast"){
		command<-paste("mafft --retree 1 --maxiterate 0 --quiet --inputorder ", filename, " > ", outname, sep="" )
	}
	system(command)
} 


# function that calculates phylogentic weights from the 16S ribosomal RNA alignment
tree_weights<-function(tree, method="clustal"){
	tip<-Ntip(tree)
	node<-Nnode(tree)
	
	di<-dist.nodes(tree)
	#root<-tip+1
	root<-setdiff(tree[[1]][,1],tree[[1]][,2])
	out<-vector("list",tip)
	names(out)<-1:tip
	for(i in 1:tip){
		temp<-numeric(0)
		t1<-i
		while(t1!=(tip+1)){
			t1<-tree[[1]][which(tree[[1]][,2]==t1),1]
			temp<-c(temp,t1)
		}
		out[[i]]<-temp
	}
	count<-table(unlist(out))
	le<-0
	for(i in 1:length(count)){
		t1<-tree[[1]][which(tree[[1]][,2]==as.numeric(names(count)[i])),1]
		le<-c(le,di[t1,as.numeric(names(count)[i])])
	}
	if(method=="clustal"){
		count<-le/count
	}
	if(method=="copra"){
		names(le)<-names(count)
		count<-le/2
		
	}
	
	
	weight<-numeric()
	
	for(i in 1:tip){
		t1<-di[as.numeric(names(out)[i]),out[[i]][1]]
		po<-match(out[[i]],names(count))
		weight<-c(weight,sum(count[po],t1))
		
	}
	names(weight)<-tree$tip.label
	weight<-weight/sum(weight)
	weight
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


build_anno<-function(ooi="NC_000911"){
	fastutr<-read.fasta("utr_seqs.fa") # for own alignments
	fastnames<-tolower(names(fastutr)) # for own alignments
	
	# create alignment of sRNAs
	sRNA_alignment<-parse_fasta(system("mafft --maxiterate 1000 --localpair --quiet input_sRNA.fa ", intern=T))
	
	wd<-getwd()
	
	# read optimal IntaRNA results from all organisms
	d2<-dir()
	f<-grep("tags.clustered$", d2)
	opt<-read.csv(d2[f], sep=";")
	opt<-as.matrix(opt)

	
	# reference for homolog clusters
	dat<-read.csv("CopraRNA2_prep_anno_addhomologs_padj_amountsamp.csv")
	dat<-as.matrix(dat)
	e<-grep("Annotation", colnames(dat))
	
	# read which Refseq IDs are present
	all_orgs<-colnames(dat)[3:(e-1)]

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
	
	
	ooi_pos<-grep(ooi, colnames(dat))
	align_pos<-1:nrow(dat)
	dat_old<-dat
	dat_all<-cbind(dat,align_pos)
		
	
	# divide the data in subsets for parallel processing 
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
	for (i in 1:max_cores) { thread2tmpfile = c(thread2tmpfile, tempfile(pattern="CopraRNA2.phyloSort.")); }

	
	# start parallel processing of phylogenetic order
	vari<-foreach(ji=1:max_cores)  %dopar% {	
		dat<-dat_all[count_vect1[ji]:count_vect2[ji],]
		all_orgs<-colnames(dat)[3:(e-1)]
		order_table<-dat[,3:(e-1)]
		colnames(order_table)<-colnames(dat)[3:(e-1)]
		order_table[]<-NA
		positions<-dat[,"align_pos"]
		order_table_list<-rep(list(order_table),length(all_orgs))
		names(order_table_list)<-all_orgs
		for(i in  1:nrow(dat)){
			e<-grep("Annotation", colnames(dat))
			temp<-dat[i,3:(e-1)]
			temp<-temp[which(temp!="")]
			orgs<-names(temp)
			if(length(temp)>1){
			
				query<-gsub("\\(.*","",as.character(temp))
				genename<-gsub(".*\\(","",as.character(temp))
				genename<-gsub("\\|.*", "", genename)
				genename<-gsub("\\/", "", genename)
				
				pos_opt<-(match(query, opt[,1]))
				
				tab<-opt[pos_opt,c("start1","end1","seedStart1","seedEnd1","start2","end2","seedStart2","seedEnd2")]
				tab<-cbind(tab,paste(genename, "_", query,sep=""),query,opt[pos_opt,"p.value"],orgs)
				colnames(tab)<-c("start","end","start_seed", "end_seed","start_sRNA","end_sRNA","start_seed_sRNA", "end_seed_sRNA","name" ,"name2","pvalue","orgs")
				
				tempf<-thread2tmpfile[ji]
				temp2<-na.omit(match(tolower(tab[,10]), fastnames))
				write.fasta(fastutr[temp2],file=tempf ,names=tab[,"orgs"])
				ca<-paste("mafft --maxiterate 1000 --localpair --quiet --inputorder ", tempf,"",sep="")
				alignment<-parse_fasta(system(ca, intern=T))
				datp<-as.DNAbin(alignment)
				dis<-dist.dna(datp, model = "F81",pairwise.deletion = T)
				dis<-as.matrix(dis)
				
				for(jj in 1:length(orgs)){
					ooip<-match(orgs[jj], colnames(dis))
					
					if(is.na(ooip)==F){
						ord<-sort(dis[,ooip],na.last=T)
						ord[which(is.na(ord))]<-1
						ord<-ord[2:length(ord)]
						dup<-which(duplicated(ord))
						
						order_table_list[[orgs[jj]]][i,match(names(sort(dis[,ooip],na.last=T)), colnames(order_table_list[[jj]]))]<-seq(1,nrow(dis))
						order_table_list[[orgs[jj]]][i,orgs[jj]]<-1
						while(length(dup)>0){
							
							tmp<-which(ord==ord[dup[1]])
							tmpd<-which(ord[dup]==ord[dup[1]])
							tmp2<-as.numeric(tab[match(names(tmp), tab[,"orgs"]),"pvalue"])
							names(tmp2)<-names(tmp)
							tmp2<-sort(tmp2)
							st<-min(as.numeric(order_table_list[[orgs[jj]]][i,match(names(tmp2),colnames(order_table_list[[orgs[jj]]]))]), na.rm=T)
							tmp3<-match(names(tmp2), colnames(order_table_list[[jj]]))
							order_table_list[[orgs[jj]]][i,tmp3]<-seq(st,length(tmp3)+st-1)
							
							dup<-dup[-tmpd]
						}
					}
					
				}
			}
		}
		
		for(jj in 1:length(order_table_list)){
			row.names(order_table_list[[jj]])<-positions
		}
		temp_out<-order_table_list
		temp_out
	}

	# cleanup temp files
	file.remove(thread2tmpfile); 

	dat<-dat_old
	dat<-cbind(dat,align_pos)
	order_table<-dat[,3:(e-1)]
	colnames(order_table)<-colnames(dat)[3:(e-1)]
	order_table[]<-NA
	order_table_list<-rep(list(order_table),length(all_orgs))
	for(i in 1:length(vari)){
		for(ii in 1:length(vari[[i]])){
			tmp<-match(rownames(vari[[i]][[ii]]),dat[,"align_pos"])
			tmp2<-match(colnames(vari[[i]][[ii]]), colnames(order_table))
			order_table_list[[ii]][tmp,tmp2]<-vari[[i]][[ii]]
		}
	}
	order_table<-order_table_list
	names(order_table)<-colnames(dat_old)[3:(e-1)]
	order_table_all_orgs<-order_table
	save(order_table_all_orgs, file="order_table_all_orgs.Rdata")
	
	
}

build_anno(ooi=ooi)
