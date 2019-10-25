

#Call:
# R --slave -f ./CopraRNA2html.r

suppressPackageStartupMessages(require(seqinr))
suppressPackageStartupMessages(require(doMC))


# get absolute path
initial.options <- commandArgs(trailingOnly = FALSE)
path<-initial.options [4]
path<-sub("CopraRNA2html.r","",path)

# register cores for parallel processing
co<-readLines("CopraRNA_option_file.txt") 
max_cores<-as.numeric(gsub("core count:","",co[grep("core count:", co)]))
registerDoMC(max_cores)

# jalview properties file
jalprops<-paste(path,"jalview_props.txt",sep="")

# CopraRNA result file 
inputfile="CopraRNA_result_all.csv"
dat<-read.csv(inputfile, sep=",")
	
# number of top predictions which should be investigated
co<-readLines("CopraRNA_option_file.txt") 
numMax<-as.numeric(gsub("top count:","",co[grep("top count:", co)]))
# ensure number does not exceed available data
num <- min(numMax, nrow(dat))


# markdown template path
markdown<-paste(path,"mardown_template.Rmd",sep="")

# preset path to required files, path can also be specified as argument
copref_path<-paste(path,"CopraRNA_available_organisms.txt",sep="")

# read the organism of interest (ooi) from the ncRNA fasta file. The sRNA of the ooi is considered to be the first sequence.
ooi<-gsub("ncRNA_","",names(read.fasta("ncrna.fa"))[1])
oois<-ooi


interactions<-function(ooi1=ooi, oois1=oois, inpfile=inputfile, number=num, outname=""){

	ooi<-ooi1
	dat<-read.csv(inpfile, sep=",")
	dat<-as.matrix(dat)
	e<-grep("Annotation", colnames(dat))
	oois<-colnames(dat)[4:(e-1)]


	int<-list()
	for(i in 1:length(oois)){
		tmp<-grep(paste(oois[i],"_upfrom.*_down_.*.fa.intarna.sorted.csv$",sep=""),dir())
		tmp<-read.csv(dir()[tmp],sep=";")
		int[[length(int)+1]]<-as.matrix(tmp)
	}

	names(int)<-oois

	selection_name1<-gsub("\\(.*","",dat[1:number,match(ooi, colnames(dat))])
	selection_name2<-gsub("\\|.*","",dat[1:number,match(ooi, colnames(dat))])
	selection_name2<-gsub(".*\\(","",selection_name2)
	selection_name<-paste(selection_name2,selection_name1,sep="_")
	
	for(ji in 1:number){
		temp<-dat[ji,4:(e-1)]
		empty<-which(temp!="")
		loc<-gsub("\\(.*","",as.character(temp)[empty])
		loc2<-names(temp)[empty]
		out<-c()
		for(ij in 1:length(loc)){
			org<-loc2[ij]
			tmp<-which(toupper(int[org][[1]][,1])==toupper(loc[ij]))
			if(length(tmp)>0){
				out<-c(out, org,"")
				for(j in 1:length(tmp)){
					temp<-as.character(int[org][[1]][tmp[j],])
					names(temp)<-colnames(int[org][[1]])
					if(j==1){
						out<-c(out,paste(loc[ij],"_optimal_interaction", " ",temp["E"] ,sep=""))
					}
					if(j==2){
						out<-c(out,paste(loc[ij],"_suboptimal_interaction"," ",temp["E"],sep=""))
					}
					inter<-strsplit(temp["hybridDP"],"&")
					m<-strsplit(inter[[1]][1],"")[[1]]
					s<-rev(strsplit(inter[[1]][2],"")[[1]])
					m_seq<-strsplit(temp["subseq1"],"")[[1]]
					s_seq<-rev(strsplit(temp["subseq2"],"")[[1]])
					
					m_seq_un<-rep(" ",length(m_seq))
					s_seq_un<-rep(" ",length(s_seq))
					
					m_bulge<-which(m==".")
					s_bulge<-which(s==".")
					names(m_bulge)<-rep("m",length(m_bulge))
					names(s_bulge)<-rep("s",length(s_bulge))
					
					bulge<-sort(c(m_bulge,s_bulge))
					
					if(length(bulge)>0){
						for(jj in 1:length(bulge)){
							if(names(bulge)[jj]=="m"){
								m_seq_un[bulge[jj]]<-m_seq[bulge[jj]]
								m_seq[bulge[jj]]<-" "
								
							}
							if(names(bulge)[jj]=="s"){
								s_seq_un[bulge[jj]]<-s_seq[bulge[jj]]
								s_seq[bulge[jj]]<-" "
								
							}
							
						}
					}
					count<-1
					while(count<=length(m_seq) & count<=length(s_seq)){
						if(m_seq[count]==" " & s_seq[count] != " "){
							s_seq<-append(s_seq," ",count-1)
							s_seq_un<-append(s_seq_un," ",count-1)
						}
						if(s_seq[count]==" " & m_seq[count] != " "){
							m_seq<-append(m_seq," ",count-1)
							m_seq_un<-append(m_seq_un," ",count-1)
						}
						count<-count+1
					}
					bind<-rep(" ", length(m_seq))
					bind[which(m_seq != " " & s_seq != " ")]<-"|"
					out<-c(out,paste(m_seq_un, collapse=""),paste(m_seq, collapse=""),paste(bind, collapse=""),paste(s_seq, collapse=""),paste(s_seq_un, collapse=""),"")	
				}
			}	
		}
		out2<-c("<style>",
		"p.small {",
		  "line-height: 0.6;",
		  "font-family: courier;",
		  "white-space: pre;",
		"  white-space: post;",
		"}",
		"</style>",
		"<p class='small'>"
		)
		out<-paste(out,"<br>", sep="")
		out<-c(out2,out,"</p>")
		if(file.exists(paste("./evo_alignments2/",selection_name[ji],sep=""))){
			write.table(out, file=paste("./evo_alignments2/",selection_name[ji],"/interactions.html",sep=""), sep="\t", row.names=F, col.names=F, quote=F)
		}
	}
}


jalview<-function(inpfile=inputfile, number=num, align_folder="./evo_alignments2"){
	dat<-read.csv(inpfile, sep=",")
	selection<-gsub("\\(.*","",dat[1:number,match(ooi, colnames(dat))])
	d<-dir(align_folder)
	
	# divide the data in subsets for parallel processing 
	max_cores<-min(number, max_cores)
	jobs<-nrow(dat)%/%max_cores
	rest<-nrow(dat)-max_cores*jobs
	jobs<-rep(jobs,max_cores)
	jobs[1:rest]<-jobs[1:rest]+1
	count_vect1<-cumsum(c(1,jobs[1:(length(jobs)-1)]))
	count_vect2<-cumsum(jobs)

	foreach(ji=1:max_cores)  %dopar% {
		for(i in count_vect1[ji]:count_vect2[ji]){
			pos<-grep(selection[i],d)
			na<-d[pos]
			fol<-paste(align_folder,"/",na,"/",sep="")
			if(file.exists(fol)){
				jal<-paste("jalview -nodisplay  -props ",jalprops, " -open ", paste(fol,na,"_mRNA_alignment.fasta", sep=""), " -features " ,paste(fol,na,"_mRNA_features.txt", sep=""), " -annotations " ,paste(fol, na,"_mRNA_annotation.txt", sep=""),  " -png  ",  paste(fol,na,"_mRNA.PNG", sep="") , sep="" )
				system(jal)
				jal<-paste("jalview -nodisplay -props ",jalprops, " -open ", paste(fol,na,"_sRNA_alignment.fasta", sep=""), " -features " ,paste(fol,na,"_sRNA_features.txt", sep=""), " -annotations " ,paste(fol, na,"_sRNA_annotation.txt", sep=""),  " -png ",  paste(fol,na,"_sRNA.PNG", sep="") , sep="" )
				system(jal)
			}
		}
	}
}




html_table<-function(ooi1=ooi, oois1=oois, inpfile=inputfile, number=num){

	ma<-readLines(markdown)
	tab<-read.csv(inpfile, sep=",")
	selection_name1<-gsub("\\(.*","",dat[1:number,match(ooi, colnames(dat))])
	selection_name2<-gsub("\\|.*","",dat[1:number,match(ooi, colnames(dat))])
	selection_name2<-gsub(".*\\(","",selection_name2)
	selection<-paste(selection_name2,selection_name1,sep="_")
	copref<-read.delim(copref_path, sep="\t", header=T,comment.char = "#")	
	tab<-as.matrix(tab)
	e<-grep("Annotation", colnames(dat))
	temp_dat<-dat[,4:(e-1)]
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

	orgs<-colnames(dat)[4:(e-1)]
	tabs<-vector("list",length(orgs))
	names(tabs)<-orgs
	nam2<-paste(nam2,orgs,sep="_")
	Annotation<-tab[1:number,e]
	
	for(i in 1:length(orgs)){
		tab2<-tab[1:number,orgs[i]]
		tab2<-gsub("\\(","\\|",tab2)
		tab2<-gsub("\\)","",tab2)
		tab2<-strsplit(tab2,"\\|")
		nu<-which(unlist(lapply(tab2,length))!=0)
		tab3<-matrix(,number,9)
		if(length(nu)>0){
			tab3[nu,]<-do.call(rbind,tab2)
		}
		tab2<-tab3
		tab2<-gsub("GeneID:","",tab2)
		
		tab3<-cbind(1:number,tab[1:number,c(1,2)],tab2, Annotation)
		colnames(tab3)<-c("Rank","CopraRNA_FDR","CopraRNA_p-value","Locus_tag","Gene_name","Energy","IntaRNA_p-value","start_mRNA","end_mRNA","start_sRNA","end_sRNA","GeneID","Annotation")
		tabs[[i]]<-tab3
	}
	
	wd<-getwd()
	co<-readLines(paste(wd,"/CopraRNA_option_file.txt",sep=""))
	noclean<-grep("noclean:", co)
	noclean<-as.numeric(gsub("noclean:","",co[noclean]))
	if(noclean==1){
		enrich1<-"./enriched_heatmap_big.pdf"
		enrich_thumb<-"./enriched_heatmap_big.png"
		aux<-"'./aux_table.csv'"
		mrna_reg1<-"./mRNA_regions_with_histogram.pdf"
		srna_reg1<-"./sRNA_regions_with_histogram.pdf"
		cons1<-"./conservation_heatmap.pdf"
	}

	if(noclean==0){
		enrich1<-"./Enrichment/enriched_heatmap_big.pdf"
		enrich_thumb<-"./Enrichment/enriched_heatmap_big.png"
		aux<-"'./Enrichment/aux_table.csv'"
		mrna_reg1<-"./Regions_plots/mRNA_regions_with_histogram.pdf"
		srna_reg1<-"./Regions_plots/sRNA_regions_with_histogram.pdf"
		cons1<-"./conservation_heatmap.pdf"
	}
	enrich_thumb<-paste0("![](",enrich_thumb,")")
	
	if(file.exists(enrich1)){
		enrich<- paste0("[Functional enrichment](",enrich1,"){target='_blank'}")
		enrich2<- paste0("![Functional enrichment](",enrich1,"#zoom=75){width=100% height=600}")
	} else{
		enrich<-""
		enrich2<-""
	}
	
	mrna_reg<-paste0("[mRNA regions plot](",mrna_reg1,"){target='_blank'}")
	mrna_reg2<-paste0("![mRNA regions plot](",mrna_reg1,"){width=100% height=600}")
	srna_reg<-paste0("[sRNA regions plot](",srna_reg1,"){target='_blank'}")
	srna_reg2<-paste0("![sRNA regions plot](",srna_reg1,"){width=100% height=600}")
	cons<-paste0("[Phylogenetic target conservation](",cons1,"){target='_blank'}")
	cons2<-paste0("![Phylogenetic target conservation](",cons1,"){width=100% height=600}")
	ma<-c(ma,paste("## Overview {.tabset .tabset-fade .tabset-pills}","### Phylogenetic target conservation",cons,"\n",cons2,"\n","### Functional enrichment",enrich,"\n",enrich2,"\n","### mRNA regions plota",mrna_reg,"\n",mrna_reg2,"\n","### sRNA regions plots",srna_reg,"\n",srna_reg2,"\n","### Auxiliary enrichment","\n",sep="\n"))
	
	prefix<-paste("```{r,echo=F}","wd<-getwd()",sep="\n")
	suffix<-paste("kable((tab)) %>% ","kable_styling(c('striped', 'bordered')) %>% ","kable_styling(fixed_thead = T) %>% ","kable_styling(full_width = F)",  "```",sep="\n")
	aux2<-c("tab<-read.csv(",aux,",sep=',')")
	
	ma<-c(ma,prefix,aux2,suffix)
	ints<-paste("./evo_alignments2/",selection,"/interactions.html",sep="")
	ints<- paste0("[interaction](",ints,"){target='_blank'}")
	#d<-dir("./jalview")
	na<-which(is.na(d[match(selection,d)]))
	mrna<-paste("./evo_alignments2/",d[match(selection,d)],"/",d[match(selection,d)],"_mRNA.PNG",sep="")
	mrna<-paste0("[mRNA_alignment](",mrna,"){target='_blank'}")
	srna<-paste("./evo_alignments2/",d[match(selection,d)],"/",d[match(selection,d)],"_sRNA.PNG",sep="")
	srna<-paste0("[sRNA_alignment](",srna,"){target='_blank'}")

	heat<-paste("./evo_alignments2/",selection,"/",selection,"_conservation_heatmap.pdf",sep="")
	#heat<-paste("./heatmaps/",selection,".pdf",sep="")
	heat<- paste0("[conservation](",heat,"){target='_blank'}")
	probcons<-paste("./evo_alignments2/",selection,"/spot_probabilities.png",sep="")
	probcons<- paste0("[conserved_peaks](",probcons,"){target='_blank'}")
	#Links<-paste(heat,mrna,srna,ints, sep=", ")
	Links<-paste(heat,mrna, srna,ints,probcons, sep=", ")
	
	if(length(na)>0){
		Links[na]<-""
	}
	ind_tables<-paste(getwd(),"/evo_alignments2/ind_tables",sep="")
	dir.create(ind_tables)
	
	
	tab<-cbind(tabs[[1]][,1], Links,tabs[[1]][,2:ncol(tabs[[1]])])
	colnames(tab)[1]<-"Rank"
	write.table(tab, file=paste(ind_tables,"/",names(tabs)[1],".txt",sep=""), sep="\t", row.names=FALSE, quote=F)
	tmp<-paste("'./evo_alignments2/ind_tables/",names(tabs)[1],".txt'",sep="")
	tmp<-c("tab<-read.csv(",tmp,",sep='\\t')")
	ma<-c(ma,"\n", paste("## Prediction organisms of interest",paste("### ", nam2[1],sep=""),sep="\n"),prefix,tmp,suffix, "## Predictions other organisms {.tabset .tabset-fade .tabset-pills}")
	
	if(length(tabs)>1){
		for(i in 2:length(tabs)){
			tab<-cbind(tabs[[i]][,1], Links,tabs[[i]][,2:ncol(tabs[[i]])])
			colnames(tab)[1]<-"Rank"
			write.table(tab, file=paste(ind_tables,"/",names(tabs)[i],".txt",sep=""), sep="\t", row.names=FALSE, quote=F)
			tmp<-paste("'./evo_alignments2/ind_tables/",names(tabs)[i],".txt'",sep="")
			tmp<-c("tab<-read.csv(",tmp,",sep='\\t')")
			ma<-c(ma,"\n", paste("### ", nam2[i],sep=""),prefix,tmp,suffix)
		}
	}
	
	writeLines(ma, con = "markdown_final.Rmd", sep = "\n", useBytes = FALSE)
	rmarkdown::render(paste(getwd(),"/markdown_final.Rmd",sep=""),output_file='CopraRNA2_result.html',intermediates_dir=getwd(),knit_root_dir=getwd(),output_dir=getwd(),clean =F)
}

interactions()	
jalview()
html_table()
