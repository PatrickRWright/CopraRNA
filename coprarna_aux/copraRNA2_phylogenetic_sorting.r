# script by Jens Georg

# dependencies: 
## CopraRNA_available_organisms.txt
## Mafft

#call:
# R --slave -f ./copraRNA2_phylogenetic_sorting.r --args copref_path=/home/jens/For_CopraRNA2.0/CopraRNA_available_organisms.txt dialign_conf=/home/jens/For_CopraRNA2.0/dialign_conf/ weight_method=clustal ribosomal_rna=16s_sequences.fa

require(phangorn)
require(seqinr)
require(doMC)

# get absolute path
initial.options <- commandArgs(trailingOnly = FALSE)
path<-initial.options [4]
path<-sub("copraRNA2_phylogenetic_sorting.r","",path)
print(path)

# preset path to required files, path can also be specified as argument
copref_path<-paste(path,"CopraRNA_available_organisms.txt",sep="")
#copref_path<-"/home/jens/For_CopraRNA2.0/CopraRNA_available_organisms.txt"


# read the organism of interest (ooi) from the ncRNA fasta file. The sRNA of the ooi is considered to be the first sequence.
ooi<-gsub("ncRNA_","",names(read.fasta("ncrna.fa"))[1])

# register cores for parallel processing
co<-readLines("CopraRNA_option_file.txt") 
co2<-grep("core count:", co)
max_cores<-as.numeric(gsub("core count:","",co[co2]))
registerDoMC(max_cores)



# method to calculate pyhlogentic weights from the 16S alignment. "clustal" = ClustalW method, "copra" = CopraRNA_1 method
weight_method="clustal"

# path to ribosomal RNA fasta
ribosomal_rna="16s_sequences.fa"

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
# combine all UTR sequences in one fasta file (needed if alignments are not available)
fast<-function(){
	require(seqinr)
	falist1<-dir()
	int<-grep("intarna", falist1)
	falist1<-falist1[-int]
	falist<-grep("_([0-9]){1,}_down_([0-9]){1,}\\.fa$",falist1)
	
	
	for(i in 1:length(falist)){
		temp<-read.fasta(falist1[falist[i]], strip.desc = TRUE)
		write.fasta(temp, names=names(temp),file="utr_seqs.fa", open="a")
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
		out[[i]]<-temp#[1:(length(temp)-1)]
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

# re-align the 16S sequences with an accurate mafft setting and calculate weights
mafft(ribosomal_rna, outname="16s_sequences_mafft_align.fas", mode="accurate")
ribo<-read.phyDat("16s_sequences_mafft_align.fas", format="fasta", type="DNA")
dm <- dist.ml(ribo, model="F81")
fitJC<- upgma(dm)
weight<-tree_weights(fitJC, method=weight_method)

build_anno<-function(ooi="NC_000911",own_alignment=F){
	if(own_alignment==TRUE){
		dir.create("evo_alignments")	
		fast()
		fastutr<-read.fasta("utr_seqs.fa") # for own alignments
		fastnames<-tolower(names(fastutr)) # for own alignments
	}
	
	# create alignment of sRNAs
	input<-paste("mafft --maxiterate 1000 --localpair --quiet", " --localpair", " --quiet input_sRNA.fa >", "aligned_sRNA.fa", sep="")
	system(input)
	sRNA_alignment<-read.fasta("aligned_sRNA.fa")
	
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
	
	#identify the column postion of the organism of interest (ooi) and exclude rows without a homolog from the ooi 
	ooi_pos<-grep(ooi, colnames(dat))
	ex<-which(dat[,ooi]!="")
	align_pos<-1:nrow(dat)
	dat_old<-dat
	dat<-dat[ex,]
	align_pos<-ex
	dat_all<-dat
	
	
	
	# divide the data in subsets for parallel processing 
	jobs<-nrow(dat)%/%max_cores
	rest<-nrow(dat)-max_cores*jobs
	jobs<-rep(jobs,max_cores)
	jobs[1]<-jobs[1]+rest
	count_vect1<-cumsum(c(1,jobs[1:(length(jobs)-1)]))
	count_vect2<-cumsum(jobs)
	

	# start parallel processing of phylogenetic order
	vari<-foreach(ji=1:max_cores)  %dopar% {	
		dat<-dat_all[count_vect1[ji]:count_vect2[ji],]
		all_orgs<-colnames(dat)[3:(e-1)]
		order_table<-dat[,3:(e-1)]
		colnames(order_table)<-colnames(dat)[3:(e-1)]
		order_table[]<-NA
			
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
				
				if(own_alignment==TRUE){
					temp2<-na.omit(match(tolower(tab[,10]), fastnames))
					write.fasta(fastutr[temp2],file="test.fa" ,names=tab[,"name2"])
					input<-paste("mafft --maxiterate 1000 --localpair --quiet"," --localpair",  " --quiet test.fa > test2.fa", sep="")
					system(input)
					alignment<-read.fasta("test2.fa")
					unlink("test2.fa")
					unlink("test.fa")
				}

				if(own_alignment==FALSE){
					alignment<-read.fasta(paste(wd,"/target_alignments/",align_pos[count_vect1[ji]+i-1], ".aln", sep=""))
					dup<-which(duplicated(names(alignment)))
					if(length(dup)>0){
						alignment<-alignment[-dup]
					}
					unlink("test2.fa")
					unlink("test.fa")
				}
				
				write.fasta(alignment, names=tab[,"orgs"], file.out=paste(align_pos[count_vect1[ji]+i-1],"_" ,tab[1,9],"_ooi_cons_pos_mRNA.fasta", sep=""))
				datp<-read.dna(paste(align_pos[count_vect1[ji]+i-1],"_" ,tab[1,9],"_ooi_cons_pos_mRNA.fasta", sep=""), format = "fasta")
				unlink(paste(align_pos[count_vect1[ji]+i-1],"_" ,tab[1,9],"_ooi_cons_pos_mRNA.fasta", sep=""))
				dis<-dist.dna(datp, model = "F81",pairwise.deletion = T)
				dis<-as.matrix(dis)
				ooip<-match(ooi, colnames(dis))
				
				if(is.na(ooip)==F){
					ord<-sort(dis[,ooip],na.last=T)
					ord[which(is.na(ord))]<-1
					ord<-ord[2:length(ord)]
					dup<-which(duplicated(ord))
					
					order_table[i,match(names(sort(dis[,ooip],na.last=T)), colnames(order_table))]<-seq(1,nrow(dis))
					order_table[i,ooi]<-1
					while(length(dup)>0){
						
						tmp<-which(ord==ord[dup[1]])
						tmpd<-which(ord[dup]==ord[dup[1]])
						tmp2<-as.numeric(tab[match(names(tmp), tab[,"orgs"]),"pvalue"])
						names(tmp2)<-names(tmp)
						tmp2<-sort(tmp2)
						st<-min(as.numeric(order_table[i,match(names(tmp2),colnames(order_table))]))
						tmp3<-match(names(tmp2), colnames(order_table))
						order_table[i,tmp3]<-seq(st,length(tmp3)+st-1)
						
						dup<-dup[-tmpd]
					}
				}
			}
		}
		
		row.names(order_table)<-dat[,ooi_pos]
		temp_out<-order_table
		temp_out
	}
	dat<-dat_old
	order_table<-dat[,3:(e-1)]
	colnames(order_table)<-colnames(dat)[3:(e-1)]
	order_table[]<-NA
	for(i in 1:length(vari)){
		tmp<-match(rownames(vari[[i]]),dat[,3])
		order_table[tmp,]<-vari[[i]]
	}
	save(order_table, file="order_table_all.Rdata")
}

build_anno(ooi=ooi,own_alignment=F)
