

#Call:
# R --slave -f ~/CopraRNA-git/coprarna_aux/CopraRNA2html.r


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

#bibliography
bib<-paste0(path,"bibliography.bib")


# markdown template path
markdown<-paste(path,"markdown_template.Rmd",sep="")

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
	selection_name<-gsub("N/A","NA",selection_name)
	
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
	jobs<-number%/%max_cores
	rest<-number-max_cores*jobs
	jobs<-rep(jobs,max_cores)
	if(rest>0){
		jobs[1:rest]<-jobs[1:rest]+1
	}
	count_vect1<-cumsum(c(1,jobs[1:(length(jobs)-1)]))
	count_vect2<-cumsum(jobs)

	dump1<-foreach(ji=1:max_cores)  %dopar% {
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
	
	opt<-readLines("CopraRNA_option_file.txt")
	ver<-grep("version:", opt)
	ver<-gsub("version:","", opt[ver])
	ver<-paste0("Version: ",ver)
	ma<-append(ma,ver, after=6)	
	tab<-read.csv(inpfile, sep=",")
	selection_name1<-gsub("\\(.*","",dat[1:number,match(ooi, colnames(dat))])
	selection_name2<-gsub("\\|.*","",dat[1:number,match(ooi, colnames(dat))])
	selection_name2<-gsub(".*\\(","",selection_name2)
	selection<-paste(selection_name2,selection_name1,sep="_")
	
	selection<-gsub("N/A","NA",selection)
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
	ooi_nam<-paste(nam[1],orgs[1])
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
		colnames(tab3)<-c("#","FDR","pVal.","Locus_tag","Name","Energy","Int.pVal.","st","en","start_sRNA","end_sRNA","GeneID","Annotation")
		tab3<-tab3[,-c(10,11,12)]
		tabs[[i]]<-tab3
	}
	
	wd<-getwd()
	co<-readLines(paste(wd,"/CopraRNA_option_file.txt",sep=""))
	noclean<-grep("noclean:", co)
	noclean<-as.numeric(gsub("noclean:","",co[noclean]))
	#if(noclean==1){
		enrich1<-"./enriched_heatmap_big.pdf"
		enrich_thumb<-"./enriched_heatmap_big.png"
		aux<-"./aux_table.csv"
		mrna_reg1<-"./mRNA_regions_with_histogram.pdf"
		srna_reg1<-"./sRNA_regions_with_histogram.pdf"
		cons1<-"./conservation_heatmap.pdf"
	#}

	# if(noclean==0){
		# enrich1<-"./Enrichment/enriched_heatmap_big.pdf"
		# enrich_thumb<-"./Enrichment/enriched_heatmap_big.png"
		# aux<-"'./Enrichment/aux_table.csv'"
		# mrna_reg1<-"./Regions_plots/mRNA_regions_with_histogram.pdf"
		# srna_reg1<-"./Regions_plots/sRNA_regions_with_histogram.pdf"
		# cons1<-"./conservation_heatmap.pdf"
	# }
	enrich_thumb<-paste0("![](",enrich_thumb,")")
	
	if(file.exists(enrich1)){
		enrich<- paste0("[Functional enrichment](",enrich1	,"){target='_blank'}")
		enrich2<- paste0("<object data='",enrich1,"#zoom=75' type='application/pdf' width='100%' height='600'><embed src='",enrich1,"#zoom=75' type='application/pdf'> <p>This browser does not support PDFs  </p></embed></object>")
	} else{
		enrich<-""
		enrich2<-""
	}

	if(file.exists(mrna_reg1)){
		mrna_reg<-paste0("[mRNA regions plot](",mrna_reg1,"){target='_blank'}")
		mrna_reg2<-paste0("<object data='",mrna_reg1,"' type='application/pdf' width='100%' height='600'><embed src='",mrna_reg1,"#' type='application/pdf'> <p>This browser does not support PDFs  </p></embed></object>")
	} else{
		mrna_reg<-""
		mrna_reg2<-""
	}
	
	if(file.exists(srna_reg1)){
		srna_reg<-paste0("[sRNA regions plot](",srna_reg1,"){target='_blank'}")
		srna_reg2<-paste0("<object data='",srna_reg1,"' type='application/pdf' width='100%' height='600'><embed src='",srna_reg1,"#' type='application/pdf'> <p>This browser does not support PDFs  </p></embed></object>")
	} else{
		srna_reg<-""
		srna_reg2<-""
	}
	
	if(file.exists(cons1)){
		cons<-paste0("[Phylogenetic target conservation](",cons1,"){target='_blank'}")
		cons2<-paste0("<object data='",cons1,"' type='application/pdf' width='100%' height='600'><embed src='",cons1,"#' type='application/pdf'> <p>This browser does not support PDFs  </p></embed></object>")
	} else{
		cons<-""
		cons2<-""
	}
	
	if(file.exists(aux)){
		aux2<-c("tab<-read.csv(",aux,",sep=',')")
	} else{
		aux2<-"tab<-matrix('no auxiliary enrichment available',1,1)"
	}
	
	ma<-c(ma,paste("## Overview {.tabset .tabset-fade .tabset-pills}","### Phylogenetic target conservation",cons2,"\n","### Functional enrichment",enrich2,"\n","### mRNA regions plots",mrna_reg2,"\n","### sRNA regions plots",srna_reg2,"\n","### Auxiliary enrichment","\n",sep="\n"))
	#ma<-c(ma,paste("## Overview {.tabset .tabset-fade .tabset-pills}","### Phylogenetic target conservation",cons2,"\n","### mRNA regions plots",mrna_reg2,"\n","### sRNA regions plots",srna_reg2,"\n",sep="\n"))
	
	prefix<-paste("```{r,echo=F}","wd<-getwd()",sep="\n")
	#suffix<-paste("kable((tab)) %>% ","kable_styling(c('striped', 'bordered')) %>% ","kable_styling(fixed_thead = T) %>% ","kable_styling(full_width = F)",  "```",sep="\n")
	suffix<-paste("datatable(tab, escape = FALSE,rownames= FALSE, extensions = 'FixedHeader',options = list(fixedHeader = TRUE),class = 'cell-border stripe',width='100%')","  ```"," ",sep="\n")
	
	
	#<a href="http://rstudio.com">RStudio</a>
	ma<-c(ma,prefix,aux2,suffix)
	ints<-paste("./evo_alignments2/",selection,"/interactions.html",sep="")
	#ints<- paste0("[interaction](",ints,"){target='_blank'}")
	ints<- paste0("<a href='",ints,"' target='popup' onclick=\"window.open('",ints,"','popup','width=600,height=600'); return false;\">Interactions</a>")
	d<-dir("./evo_alignments2")
	na<-which(is.na(d[match(selection,d)]))
	mrna<-paste("./evo_alignments2/",d[match(selection,d)],"/",d[match(selection,d)],"_mRNA.PNG",sep="")
	#mrna<-paste0("[mRNA_alignment](",mrna,"){target='_blank'}")
	mrna<-paste0("<a href='",mrna,"' target='popup' onclick=\"window.open('",mrna,"','popup','width=600,height=600'); return false;\">mRNA_alignment</a>")
	srna<-paste("./evo_alignments2/",d[match(selection,d)],"/",d[match(selection,d)],"_sRNA.PNG",sep="")
	#srna<-paste0("[sRNA_alignment](",srna,"){target='_blank'}")
	srna<-paste0("<a href='",srna,"' target='popup' onclick=\"window.open('",srna,"','popup','width=600,height=600'); return false;\">sRNA_alignment</a>")

	#heat<-paste("./evo_alignments2/",selection,"/tree.pdf",sep="")
	heat<-paste("./evo_alignments2/",selection,"/",selection,"_conservation_heatmap.pdf",sep="")
	#heat<-paste("./heatmaps/",selection,".pdf",sep="")
	#heat<- paste0("[conservation](",heat,"){target='_blank'}")
	heat<- paste0("<a href='",heat,"' target='popup' onclick=\"window.open('",heat,"','popup','width=600,height=600'); return false;\">conservation</a>")
	#probcons<-paste("./evo_alignments2/",selection,"/", selection,".html",sep="")
	probcons<-paste("./evo_alignments2/",selection,"/spot_probabilities.png",sep="")
	#probcons<- paste0("[conserved_peaks](",probcons,"){target='_blank'}")
	probcons<- paste0("<a href='",probcons,"' target='popup' onclick=\"window.open('",probcons,"','popup','width=600,height=600'); return false;\">conserved_peaks</a>")
	#Links<-paste(heat,mrna,srna,ints, sep=", ")
	
	tree<-paste("./evo_alignments2/",selection,"/tree.pdf",sep="")
	#probcons<- paste0("[conserved_peaks](",probcons,"){target='_blank'}")
	tree<- paste0("<a href='",tree,"' target='popup' onclick=\"window.open('",tree,"','popup','width=600,height=600'); return false;\">tree</a>")
	#Links<-paste(heat,mrna,srna,ints, sep=", ")
	
	
	
	Links<-paste(heat,mrna, srna,ints,probcons,tree, sep=", ")
	
	if(length(na)>0){
		Links[na]<-""
	}
	ind_tables<-paste(getwd(),"/evo_alignments2/ind_tables",sep="")
	dir.create(ind_tables)
	
	
	tab<-cbind(tabs[[1]][,1], Links,tabs[[1]][,2:ncol(tabs[[1]])])
	colnames(tab)[1]<-" "
	write.table(tab, file=paste(ind_tables,"/",names(tabs)[1],".txt",sep=""), sep="\t", row.names=FALSE, quote=F)
	tmp<-paste("'./evo_alignments2/ind_tables/",names(tabs)[1],".txt'",sep="")
	tmp<-c("tab<-read.csv(",tmp,",sep='\\t')", "\n","tab[,3]<-formatC(as.numeric(tab[,3]), format = 'e', digits = 2)", "\n","tab[,4]<-formatC(as.numeric(tab[,4]), format = 'e', digits = 2)", "\n","tab[,8]<-formatC(as.numeric(tab[,8]), format = 'e', digits = 2)" )
	ma<-c(ma,"\n", paste("## Prediction organisms of interest",paste("### ", ooi_nam,sep=""),sep="\n"),prefix,tmp,suffix, "## Predictions other organisms {.tabset .tabset-fade .tabset-pills}")
	
	load("copra_results_all.Rdata")
	
	if(length(copra_results)>1){
		for(i in 2:length(copra_results)){
				fdr<-copra_results[[i]][1:number,1]
				pval<-copra_results[[i]][1:number,2]
				ranks<-1:number
				annotation<-copra_results[[i]][1:number,"Annotation"]
				tab2<-copra_results[[i]][1:number,names(copra_results)[i]]
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
		
				tab3<-cbind(ranks,fdr,pval,tab2, annotation)
				colnames(tab3)<-c("#","FDR","pVal.","Locus_tag","Name","Energy","Int.pVal.","st","en","start_sRNA","end_sRNA","GeneID","Annotation")
				tab3<-tab3[,-c(10,11,12)]
				tab<-tab3
			
				#tab<-cbind(ranks,fdr,pval,tab)
				#colnames(tab)[1]<-" "
				write.table(tab, file=paste(ind_tables,"/",names(tabs)[i],".txt",sep=""), sep="\t", row.names=FALSE, quote=F)
				tmp<-paste("'./evo_alignments2/ind_tables/",names(tabs)[i],".txt'",sep="")
				tmp<-c("tab<-read.csv(",tmp,",sep='\\t')", "\n","tab[,2]<-formatC(as.numeric(tab[,2]), format = 'e', digits = 2)", "\n","tab[,3]<-formatC(as.numeric(tab[,3]), format = 'e', digits = 2)", "\n","tab[,7]<-formatC(as.numeric(tab[,7]), format = 'e', digits = 2)" )
				ma<-c(ma,"\n", paste("### ", nam2[i],sep=""),prefix,tmp,suffix)
		}
	}
	ma<-c(ma, "## When using CopraRNA please cite:
	
<font size='-1'>
				
Patrick R. Wright, Andreas S. Richter, Kai Papenfort, Martin Mann, Joerg Vogel, Wolfgang R. Hess, Rolf Backofen and Jens Georg  
[Comparative genomics boosts target prediction for bacterial small RNAs](https://www.pnas.org/content/110/37/E3487){target='_blank'}  
Proc Natl Acad Sci USA, 2013, 110 (37), E3487-E3496.

Patrick R. Wright, Jens Georg, Martin Mann, Dragos A. Sorescu, Andreas S. Richter, Steffen Lott, Robert Kleinkauf, Wolfgang R. Hess, and Rolf Backofen  
[CopraRNA and IntaRNA: predicting small RNA targets, networks and interaction domains](https://academic.oup.com/nar/article/42/W1/W119/2435325){target='_blank'}  
Nucleic Acids Research, 2014, 42 (W1), W119-W123.

Martin Raden, Syed M Ali, Omer S Alkhnbashi, Anke Busch, Fabrizio Costa, Jason A Davis, Florian Eggenhofer, Rick Gelhausen, Jens Georg, Steffen Heyne, Michael Hiller, Kousik Kundu, Robert Kleinkauf, Steffen C Lott, Mostafa M Mohamed, Alexander Mattheis, Milad Miladi, Andreas S Richter, Sebastian Will, Joachim Wolff, Patrick R Wright, and Rolf Backofen  
[Freiburg RNA tools: a central online resource for RNA-focused research and teaching](https://academic.oup.com/nar/article/46/W1/W25/5000013){target='_blank'}  
Nucleic Acids Research, 46(W1), W25-W29, 2018.

</font>"
		)
	
	ma<-c(ma,"\n", "## CopraRNA parameter", "\n", paste(opt,"\n",sep=""))
	writeLines(ma, con = "markdown_final.Rmd", sep = "\n", useBytes = FALSE)
	rmarkdown::render("./markdown_final.Rmd",output_file='CopraRNA2_result.html',intermediates_dir=getwd(),knit_root_dir=getwd(),output_dir=getwd(),clean =F)
}

interactions()	
jalview()
html_table()



# create clean html zip folder
# evo_alignments2



#co<-readLines("CopraRNA_option_file.txt")
#noclean<-grep("noclean:", co)
#noclean<-as.numeric(gsub("noclean:","",co[noclean]))
dir.create("./results_html")
# if(noclean==0){
	# dir.create("./results_html/Enrichment")
	# dir.create("./results_html/Regions_plots")
	# file.copy("./Enrichment/enriched_heatmap_big.pdf", "./results_html/Enrichment/enriched_heatmap_big.pdf")
	# file.copy("./Enrichment/enriched_heatmap_big.png", "./results_html/Enrichment/enriched_heatmap_big.png")
	# file.copy("./Enrichment/aux_table.csv", "./results_html/Enrichment/aux_table.csv")
	# file.copy("./Regions_plots/mRNA_regions_with_histogram.pdf", "./results_html/Regions_plots/mRNA_regions_with_histogram.pdf")
	# file.copy("./Regions_plots/sRNA_regions_with_histogram.pdf", "./results_html/Regions_plots/sRNA_regions_with_histogram.pdf")
	# file.copy("./conservation_heatmap.pdf", "./results_html/conservation_heatmap.pdf")
# }

# if(noclean==1){
	file.copy("./enriched_heatmap_big.pdf", "./results_html/enriched_heatmap_big.pdf")
	file.copy("./enriched_heatmap_big.png", "./results_html/enriched_heatmap_big.png")
	file.copy("./aux_table.csv", "./results_html/aux_table.csv")
	file.copy("./mRNA_regions_with_histogram.pdf", "./results_html/mRNA_regions_with_histogram.pdf")
	file.copy("./sRNA_regions_with_histogram.pdf", "./results_html/sRNA_regions_with_histogram.pdf")
	file.copy("./conservation_heatmap.pdf", "./results_html/conservation_heatmap.pdf")
# }
file.copy("./evo_alignments2/","./results_html/", recursive=T)


d<-dir("./results_html/evo_alignments2/")
d<-d[-which(d=="ind_tables")]

for(i in d){
	system(paste("rm ./results_html/evo_alignments2/",i,"/*.txt",sep=""))
	system(paste("rm ./results_html/evo_alignments2/",i,"/*.fasta",sep=""))
}
for(i in d){
	system(paste("rm ./evo_alignments2/",i,"/*.txt",sep=""))
	system(paste("rm ./evo_alignments2/",i,"/*.fasta",sep=""))
	system(paste("rm ./evo_alignments2/",i,"/interactions.html",sep=""))
}

file.copy("./CopraRNA2_result.html", "./results_html/CopraRNA2_result.html")

system("zip -r copra_html.zip ./results_html")
system("rm -r ./results_html")
system("rm -r ./evo_alignments2/ind_tables")
system("rm CopraRNA2_result.html")
system("rm *.md")


