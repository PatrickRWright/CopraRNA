# Parameter:
# ooi as refseq id
# script by Jens Georg

# dependency: CopraRNA_available_organisms.txt


# R --slave -f ../copraRNA2_position_script_for_evo_precalculated_alignments_w_ooi_fast3.r --args NC_000913

#require(data.table)

args <- commandArgs(trailingOnly = TRUE) 
ooi2 <- args[1] 





peakFind<-function(interaction_site_distr, thres=0.4){
	m<-max(interaction_site_distr)
	th<-m*thres
	mea<-which(interaction_site_distr==m)
	
	s<-mea[1]
	while(interaction_site_distr[s]>=th & s >= 1){
	
		s<-s-1
		if(s==0){
		break
		}
	}
	s<-s+1
	
	e<-mea[length(mea)]
	while(interaction_site_distr[e]>=th & e <= length(interaction_site_distr)){
		e<-e+1
		
	}
	e<-e-1
	
	
	peak<-c(s,e)
	}


	
	

peakFind2<-function(interaction_site_distr, thres=0.4, thres2=0.9){
	sorted<-sort(unique(interaction_site_distr), decreasing=T)
	s1<-sorted[1]
	s2<-sorted[2]
	eq<-s2/s1
	
	#m<-max(interaction_site_distr)
	th<-thres*s1
	th2<-thres*s2
	mea<-which(interaction_site_distr==s1)
	mea2<-which(interaction_site_distr==s2)
	s<-mea[1]
	s2<-mea2[1]
	while(interaction_site_distr[s]>=th & s >= 1){
	
		s<-s-1
		if(s==0){
		break
		}
	}
	s<-s+1
	
	e<-mea[length(mea)]
	while(interaction_site_distr[e]>=th & e <= length(interaction_site_distr)){
		e<-e+1
		
	}
	e<-e-1
	
	
	peak<-c(s,e)
	peak2<-c()
	ov<-intersect(mea2,seq(s,e))
	
	if(length(ov)==0){
		if(eq>=thres2){
			while(interaction_site_distr[s2]>=th2 & s2 >= 1){
				s2<-s2-1
				if(s2==0){
					break
				}
			}
			s2<-s2+1
			e2<-mea2[length(mea2)]
			while(interaction_site_distr[e2]>=th & e2 <= length(interaction_site_distr)){
				e2<-e2+1
			}
			e2<-e2-1
			peak2<-c(s2,e2)
		}
		
	}
	
	peak<-list(peak,peak2)
	peak
	}	
	
	
	
	
overlap<-function(peak, start_interaction, end_interaction, thres=0.7){
	out<-rep(FALSE,length(start_interaction))
	if(length(peak)>0){
	for(i in 1:length(start_interaction)){
		temp<-seq(start_interaction[i], end_interaction[i])
		site<-seq(peak[1],peak[2])
		inter<-intersect(temp,site)
		over<-length(inter)/min(length(temp),length(site))
		if(over>=thres){
			out[i]<-TRUE
		}
	}
	}
	out
}	
	




	

build_anno<-function(ooi="NC_000911"){
	
	require(seqinr)
	#fastutr<-read.fasta("utr_seqs.fa")
	#fastnames<-tolower(names(fastutr))
	input<-paste("mafft --maxiterate 1000 --localpair --quiet", " --localpair", " --quiet input_sRNA.fa >", "aligned_sRNA.fa", sep="")
	system(input)
	sRNA_alignment<-read.fasta("aligned_sRNA.fa")
	s_names<-names(sRNA_alignment)
	dir.create("evo_alignments")
	wd<-getwd()
	
	d2<-dir()
	f<-grep("tags.clustered$", d2)
	opt<-read.csv(d2[f], sep=";")
	#dat<-read.csv("CopraRNA2_final_all_evo_table2.csv", sep="\t",header=T)
	#dat<-dat[seq(1,num),]
	#dat<-as.matrix(dat)
	dat<-read.csv("CopraRNA2_prep_anno_addhomologs_padj_amountsamp.csv")
	dat<-as.matrix(dat)
	subopt1<-grep("_subopt.intarna.csv",d2)
	subopt<-c()
	for(i in 1:length(subopt1)){
		temp<-read.csv(d2[subopt1[i]], sep=";")
		subopt<-rbind(subopt,temp)
	}
	subopt[,1]<-tolower(subopt[,1])
	copref<-read.delim("CopraRNA_available_organisms.txt", sep="\t", header=T,comment.char = "#")	
	
	e<-grep("Annotation", colnames(dat))
		#temp<-dat[i,5:(e-1)]
		temp_dat<-dat[,3:(e-1)]
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
	conservation_table<-dat[,3:(e-1)]
	colnames(conservation_table)<-colnames(dat)[3:(e-1)]
	conservation_table[]<-NA
	#p_table<-conservation_table
	#p_table_sub<-p_table
	conservation_table_sub<-conservation_table
	conservation_table_ooi<-conservation_table
	conservation_table_sub_ooi<-conservation_table
	
	for(i in 1:nrow(dat)){
		
		print(i)
		
		
		
		e<-grep("Annotation", colnames(dat))
		#temp<-dat[i,5:(e-1)]
		temp<-dat[i,3:(e-1)]
		temp<-temp[which(temp!="")]
		orgs<-names(temp)
		if(length(temp)>1){
		#tab<-matrix(,length(temp),12)
		
		#c("darkolivegreen1","olivedrab1","olivedrab2","olivedrab3","olivedrab4","orangered1")
		#tabsub<-matrix(,length(temp),12)
		#tabsub[]<-NA
		#colnames(tabsub)<-c("start","end","start_seed", "end_seed","start_sRNA","end_sRNA","start_seed_sRNA", "end_seed_sRNA","name" ,"name2","pvalue","orgs")
		#tab[,12]<-orgs
		#tabsub[,12]<-orgs
		poslist<-c()
		query_list<-c()
		exist<-na.omit(match(ooi,orgs))
		anno_mRNA<-c("p<=0001\tD9F69C", "p<=001\tC2F161","p<01\tBAEF4D", "p<02\tB3EE3A","p<03\tA1D634",  "p>03\tFF4500")
		anno_sRNA<-c("p<=0001\tD9F69C", "p<=001\tC2F161","p<01\tBAEF4D", "p<02\tB3EE3A","p<03\tA1D634",  "p>03\tFF4500")
		query<-gsub("\\(.*","",as.character(temp))
		genename<-gsub(".*\\(","",as.character(temp))
		genename<-gsub("\\|.*", "", genename)
		genename<-gsub("\\/", "", genename)
		
		pos_opt<-(match(query, opt[,1]))
		
		tab<-opt[pos_opt,c("start1","end1","seedStart1","seedEnd1","start2","end2","seedStart2","seedEnd2")]
		tab<-cbind(tab,paste(genename, "_", query,sep=""),query,opt[pos_opt,"p.value"],orgs)
		colnames(tab)<-c("start","end","start_seed", "end_seed","start_sRNA","end_sRNA","start_seed_sRNA", "end_seed_sRNA","name" ,"name2","pvalue","orgs")
		
		pos_sub<-(match(query, (subopt[,1])))
		exist_sub<-which(is.na(pos_sub)==F)
		pos_sub<-na.omit(pos_sub)
		tabsub<-matrix(,1,12)
		if(length(pos_sub)>0){
			tabsub<-subopt[pos_sub,c("start1","end1","seedStart1","seedEnd1","start2","end2","seedStart2","seedEnd2")]
			tabsub<-cbind(tabsub,paste(genename[exist_sub], "_", query[exist_sub],sep=""),query[exist_sub],subopt[pos_sub,"p.value"],orgs[exist_sub])
			colnames(tabsub)<-c("start","end","start_seed", "end_seed","start_sRNA","end_sRNA","start_seed_sRNA", "end_seed_sRNA","name" ,"name2","pvalue","orgs")
		}
	
		tabsub<-na.omit(tabsub)
		#temp2<-na.omit(match(tolower(tab[,10]), fastnames))
		#write.fasta(fastutr[temp2],file="test.fa" ,names=orgs)
		#input<-paste("mafft --maxiterate 1000 --localpair --quiet"," --localpair",  " --quiet test.fa > test2.fa", sep="")
		#system(input)
		pos<-as.numeric(names(sort(table(poslist), decreasing=T))[1])
		#print(pos)
		#alignment<-read.fasta(paste(wd,"/target_alignments/",i, ".aln", sep=""))
		alignment<-read.fasta(paste(wd,"/target_alignments/",i, ".aln", sep=""))
		unlink("test2.fa")
		unlink("test.fa")
		dir.create(paste(wd,"/","evo_alignments/",i,"_",tab[1,9],sep=""))
		m_names<-na.omit( match(orgs, names(dat[i,3:(e-1)])))
		s_names2<-na.omit(match(orgs, s_names))
		sRNA_alignment2<-sRNA_alignment[s_names2]
		write.fasta(sRNA_alignment2, paste(nam2[s_names2],(orgs),tab[,10],sep="_"), file.out=paste(wd,"/evo_alignments/",i,"_" ,tab[1,9],"/", i,"_" ,tab[1,9], "_sRNA.fasta", sep=""))
		write.fasta(alignment, paste(nam2[m_names],(orgs),tab[,10],sep="_"), file.out=paste(wd,"/evo_alignments/",i,"_" ,tab[1,9],"/", i,"_" ,tab[1,9], "_mRNA.fasta", sep=""))
		#file.copy(paste(wd,"/alignments/", pos, "_linsi", sep=""),paste(wd,"/evo_alignments/",i,"_" ,tab[1,9],"/", i,"_" ,tab[1,9], "_mRNA.fasta", sep=""))
		tab_aligned<-tab
		tabsub_aligned<-tabsub
		 if(nrow(tabsub)>0){
		 for(j in 1:nrow(tabsub)){
		  #print(paste("align",j))
			align_num_sRNA<-which(names(sRNA_alignment2)==tabsub[j,"orgs"])
			align_num<-which(names(alignment)==tabsub[j,"name2"])
			 nogaps<-which(alignment[[align_num]]!="-")
			 nogaps2<-which(sRNA_alignment2[[align_num_sRNA]]!="-")
			 tabsub_aligned[j,1]<-nogaps[as.numeric(tabsub[j,1])]
			 tabsub_aligned[j,2]<-nogaps[as.numeric(tabsub[j,2])]
			 #tabsub_aligned[j,3]<-nogaps[as.numeric(tabsub[j,3])]
			# tabsub_aligned[j,4]<-nogaps[as.numeric(tabsub[j,4])]
			 tabsub_aligned[j,5]<-nogaps2[as.numeric(tabsub[j,5])]
			 tabsub_aligned[j,6]<-nogaps2[as.numeric(tabsub[j,6])]
			 #tabsub_aligned[j,7]<-nogaps2[as.numeric(tabsub[j,7])]
			 #tabsub_aligned[j,8]<-nogaps2[as.numeric(tabsub[j,8])]
			 
			 }
			
			a1<-which(as.numeric(tabsub[,11])<=0.001)
			if(length(a1)>0){
				m<-paste(tabsub[a1,11],"\t", paste(nam2[m_names[a1]],(orgs)[a1],tabsub[a1,10],sep="_"), "\t-1\t", tabsub[a1,1], "\t", tabsub[a1,2] ,"\tp<=0001", sep="")
				s<-paste(tabsub[a1,11],"\t", paste(nam2[m_names[a1]],(orgs)[a1],tabsub[a1,10],sep="_"), "\t-1\t", tabsub[a1,5], "\t", tabsub[a1,6] ,"\tp<=0001", sep="")
				anno_mRNA<-c(anno_mRNA, m)
				anno_sRNA<-c(anno_sRNA, s)
			}
			
			a1<-which(as.numeric(tabsub[,11])>0.001 & as.numeric(tabsub[,11])<=0.01)
			if(length(a1)>0){
				m<-paste(tabsub[a1,11],"\t", paste(nam2[m_names[a1]],(orgs)[a1],tabsub[a1,10],sep="_"), "\t-1\t", tabsub[a1,1], "\t", tabsub[a1,2] ,"\tp<=001", sep="")
				s<-paste(tabsub[a1,11],"\t", paste(nam2[m_names[a1]],(orgs)[a1],tabsub[a1,10],sep="_"), "\t-1\t", tabsub[a1,5], "\t", tabsub[a1,6] ,"\tp<=001", sep="")
				anno_mRNA<-c(anno_mRNA, m)
				anno_sRNA<-c(anno_sRNA, s)
			}
			
			a1<-which(as.numeric(tabsub[,11])>0.01 & as.numeric(tabsub[,11])<=0.1)
			if(length(a1)>0){
				m<-paste(tabsub[a1,11],"\t", paste(nam2[m_names[a1]],(orgs)[a1],tabsub[a1,10],sep="_"), "\t-1\t", tabsub[a1,1], "\t", tabsub[a1,2] ,"\tp<=01", sep="")
				s<-paste(tabsub[a1,11],"\t", paste(nam2[m_names[a1]],(orgs)[a1],tabsub[a1,10],sep="_"), "\t-1\t", tabsub[a1,5], "\t", tabsub[a1,6] ,"\tp<=01", sep="")
				anno_mRNA<-c(anno_mRNA, m)
				anno_sRNA<-c(anno_sRNA, s)
			}
			
			a1<-which(as.numeric(tabsub[,11])>0.1 & as.numeric(tabsub[,11])<=0.2)
			if(length(a1)>0){
				m<-paste(tabsub[a1,11],"\t", paste(nam2[m_names[a1]],(orgs)[a1],tabsub[a1,10],sep="_"), "\t-1\t", tabsub[a1,1], "\t", tabsub[a1,2] ,"\tp<=02", sep="")
				s<-paste(tabsub[a1,11],"\t", paste(nam2[m_names[a1]],(orgs)[a1],tabsub[a1,10],sep="_"), "\t-1\t", tabsub[a1,5], "\t", tabsub[a1,6] ,"\tp<=02", sep="")
				anno_mRNA<-c(anno_mRNA, m)
				anno_sRNA<-c(anno_sRNA, s)
			}
			a1<-which(as.numeric(tabsub[,11])>0.2 & as.numeric(tabsub[,11])<=0.3)
			if(length(a1)>0){
				m<-paste(tabsub[a1,11],"\t", paste(nam2[m_names[a1]],(orgs)[a1],tabsub[a1,10],sep="_"), "\t-1\t", tabsub[a1,1], "\t", tabsub[a1,2] ,"\tp<=03", sep="")
				s<-paste(tabsub[a1,11],"\t", paste(nam2[m_names[a1]],(orgs)[a1],tabsub[a1,10],sep="_"), "\t-1\t", tabsub[a1,5], "\t", tabsub[a1,6] ,"\tp<=03", sep="")
				anno_mRNA<-c(anno_mRNA, m)
				anno_sRNA<-c(anno_sRNA, s)
			}
			a1<-which(as.numeric(tabsub[,11])>0.3 )
			if(length(a1)>0){
				m<-paste(tabsub[a1,11],"\t", paste(nam2[m_names[a1]],(orgs)[a1],tabsub[a1,10],sep="_"), "\t-1\t", tabsub[a1,1], "\t", tabsub[a1,2] ,"\tp>03", sep="")
				s<-paste(tabsub[a1,11],"\t", paste(nam2[m_names[a1]],(orgs)[a1],tabsub[a1,10],sep="_"), "\t-1\t", tabsub[a1,5], "\t", tabsub[a1,6] ,"\tp>03", sep="")
				anno_mRNA<-c(anno_mRNA, m)
				anno_sRNA<-c(anno_sRNA, s)
			}
			
		 
		 }
		 
		 for(j in 1:nrow(tab)){
		#print(paste("align2",j))
			align_num_sRNA<-which(names(sRNA_alignment2)==tab[j,"orgs"])
			align_num<-which(names(alignment)==tab[j,"name2"])
			 nogaps<-which(alignment[[align_num]]!="-")
			 nogaps2<-which(sRNA_alignment2[[align_num_sRNA]]!="-")
			 tab_aligned[j,1]<-nogaps[as.numeric(tab[j,1])]
			 tab_aligned[j,2]<-nogaps[as.numeric(tab[j,2])]
			 #tab_aligned[j,3]<-nogaps[as.numeric(tab[j,3])]
			#tab_aligned[j,4]<-nogaps[as.numeric(tab[j,4])]
			 tab_aligned[j,5]<-nogaps2[as.numeric(tab[j,5])]
			 tab_aligned[j,6]<-nogaps2[as.numeric(tab[j,6])]
			 #tab_aligned[j,7]<-nogaps2[as.numeric(tab[j,7])]
			 #tab_aligned[j,8]<-nogaps2[as.numeric(tab[j,8])]
			 
			}
			
			a1<-which(as.numeric(tab[,11])<=0.001)
			if(length(a1)>0){
				m<-paste(tab[a1,11],"\t", paste(nam2[m_names[a1]],(orgs)[a1],tab[a1,10],sep="_"), "\t-1\t", tab[a1,1], "\t", tab[a1,2] ,"\tp<=0001", sep="")
				s<-paste(tab[a1,11],"\t", paste(nam2[m_names[a1]],(orgs)[a1],tab[a1,10],sep="_"), "\t-1\t", tab[a1,5], "\t", tab[a1,6] ,"\tp<=0001", sep="")
				anno_mRNA<-c(anno_mRNA, m)
				anno_sRNA<-c(anno_sRNA, s)
			}
			
			a1<-which(as.numeric(tab[,11])>0.001 & as.numeric(tab[,11])<=0.01)
			if(length(a1)>0){
				m<-paste(tab[a1,11],"\t", paste(nam2[m_names[a1]],(orgs)[a1],tab[a1,10],sep="_"), "\t-1\t", tab[a1,1], "\t", tab[a1,2] ,"\tp<=001", sep="")
				s<-paste(tab[a1,11],"\t", paste(nam2[m_names[a1]],(orgs)[a1],tab[a1,10],sep="_"), "\t-1\t", tab[a1,5], "\t", tab[a1,6] ,"\tp<=001", sep="")
				anno_mRNA<-c(anno_mRNA, m)
				anno_sRNA<-c(anno_sRNA, s)
			}
			
			a1<-which(as.numeric(tab[,11])>0.01 & as.numeric(tab[,11])<=0.1)
			if(length(a1)>0){
				m<-paste(tab[a1,11],"\t", paste(nam2[m_names[a1]],(orgs)[a1],tab[a1,10],sep="_"), "\t-1\t", tab[a1,1], "\t", tab[a1,2] ,"\tp<=01", sep="")
				s<-paste(tab[a1,11],"\t", paste(nam2[m_names[a1]],(orgs)[a1],tab[a1,10],sep="_"), "\t-1\t", tab[a1,5], "\t", tab[a1,6] ,"\tp<=01", sep="")
				anno_mRNA<-c(anno_mRNA, m)
				anno_sRNA<-c(anno_sRNA, s)
			}
			
			a1<-which(as.numeric(tab[,11])>0.1 & as.numeric(tab[,11])<=0.2)
			if(length(a1)>0){
				m<-paste(tab[a1,11],"\t", paste(nam2[m_names[a1]],(orgs)[a1],tab[a1,10],sep="_"), "\t-1\t", tab[a1,1], "\t", tab[a1,2] ,"\tp<=02", sep="")
				s<-paste(tab[a1,11],"\t", paste(nam2[m_names[a1]],(orgs)[a1],tab[a1,10],sep="_"), "\t-1\t", tab[a1,5], "\t", tab[a1,6] ,"\tp<=02", sep="")
				anno_mRNA<-c(anno_mRNA, m)
				anno_sRNA<-c(anno_sRNA, s)
			}
			a1<-which(as.numeric(tab[,11])>0.2 & as.numeric(tab[,11])<=0.3)
			if(length(a1)>0){
				m<-paste(tab[a1,11],"\t", paste(nam2[m_names[a1]],(orgs)[a1],tab[a1,10],sep="_"), "\t-1\t", tab[a1,1], "\t", tab[a1,2] ,"\tp<=03", sep="")
				s<-paste(tab[a1,11],"\t", paste(nam2[m_names[a1]],(orgs)[a1],tab[a1,10],sep="_"), "\t-1\t", tab[a1,5], "\t", tab[a1,6] ,"\tp<=03", sep="")
				anno_mRNA<-c(anno_mRNA, m)
				anno_sRNA<-c(anno_sRNA, s)
			}
			a1<-which(as.numeric(tab[,11])>0.3 )
			if(length(a1)>0){
				m<-paste(tab[a1,11],"\t", paste(nam2[m_names[a1]],(orgs)[a1],tab[a1,10],sep="_"), "\t-1\t", tab[a1,1], "\t", tab[a1,2] ,"\tp>03", sep="")
				s<-paste(tab[a1,11],"\t", paste(nam2[m_names[a1]],(orgs)[a1],tab[a1,10],sep="_"), "\t-1\t", tab[a1,5], "\t", tab[a1,6] ,"\tp>03", sep="")
				anno_mRNA<-c(anno_mRNA, m)
				anno_sRNA<-c(anno_sRNA, s)
			}
			
			
			# seed<-seq(as.numeric(tab[j,3]),as.numeric(tab[j,4]))
			# site<-seq(as.numeric(tab[j,1]),as.numeric(tab[j,2]))
			
			# ma<-match(seed, site)[c(1,length(seed))]
			# ms<-paste("Binding site\t", paste((orgs)[j],tab[j,10],sep="_"), "\t-1\t", tab[j,3], "\t", tab[j,4] ,"\tseed", sep="")
			
			# m1<-paste("Binding site\t", paste((orgs)[j],tab[j,10],sep="_"), "\t-1\t", tab[j,1], "\t", site[ma[1]]-1 ,"\topt", sep="")
			# m2<-paste("Binding site\t", paste((orgs)[j],tab[j,10],sep="_"), "\t-1\t", site[ma[2]]+1, "\t", tab[j,2] ,"\topt", sep="")
			
			
			# seed<-seq(as.numeric(tab[j,7]),as.numeric(tab[j,8]))
			# site<-seq(as.numeric(tab[j,5]),as.numeric(tab[j,6]))
			# ma<-match(seed, site)[c(1,length(seed))]
			
			# ms_s<-paste("Binding site\t", paste((orgs)[j],tab[j,10],sep="_"), "\t-1\t", tab[j,7], "\t", tab[j,8] ,"\tseed", sep="")
			# m_s1<-paste("Binding site\t", paste((orgs)[j],tab[j,10],sep="_"), "\t-1\t", tab[j,5], "\t", site[ma[1]]-1 ,"\topt", sep="")
			# m_s2<-paste("Binding site\t", paste((orgs)[j],tab[j,10],sep="_"), "\t-1\t", site[ma[2]]+1, "\t", tab[j,6] ,"\topt", sep="")
			
			# anno_mRNA<-c(anno_mRNA, m1,m2, ms)
			# anno_sRNA<-c(anno_sRNA, m_s1,m_s2, ms_s)
			
		 #}
		 
		 tab_aligned[,11]<-1-as.numeric(tab_aligned[,11])
		tab<-as.matrix(tab)
		 
		 
		 align_table_mRNA<-matrix(,length(alignment[[1]]),(nrow(tab)+nrow(tabsub)))
		 align_table_sRNA<-matrix(,length(sRNA_alignment[[1]]),(nrow(tab)+nrow(tabsub)))
		 align_table_mRNA[]<-0
		 align_table_sRNA[]<-0
		 for(j in 1:nrow(tab_aligned)){
		# print(paste("distri",j))
			#m1<-seq(as.numeric(tab_aligned[j,1]),as.numeric(tab_aligned[j,2]))
			align_table_mRNA[as.numeric(tab_aligned[j,1]):as.numeric(tab_aligned[j,2]),j]<-as.numeric(tab_aligned[j,11])
			#m1<-seq(as.numeric(tab_aligned[j,5]),as.numeric(tab_aligned[j,6]))
			align_table_sRNA[as.numeric(tab_aligned[j,5]):as.numeric(tab_aligned[j,6]),j]<-as.numeric(tab_aligned[j,11])
			
			
		 }
		 if(nrow(tabsub_aligned)>0){
		 tabsub_aligned[,11]<-1-as.numeric(tabsub_aligned[,11])
		tabsub<-as.matrix(tabsub)
		 for(j in 1:nrow(tabsub_aligned)){
		 #print(j)
			
			count<-nrow(tab_aligned)
			#m1<-seq(as.numeric(tabsub_aligned[j,1]),as.numeric(tabsub_aligned[j,2]))
			align_table_mRNA[as.numeric(tabsub_aligned[j,1]):as.numeric(tabsub_aligned[j,2]),j+count]<-as.numeric(tabsub_aligned[j,11])
			#m1<-seq(as.numeric(tabsub_aligned[j,5]),as.numeric(tabsub_aligned[j,6]))
			align_table_sRNA[as.numeric(tabsub_aligned[j,5]):as.numeric(tabsub_aligned[j,6]),j+count]<-as.numeric(tabsub_aligned[j,11])
			
			
		 }
		 }
		 
		 
		 align_table_sRNA<-rowSums(align_table_sRNA)
		 align_table_mRNA<-rowSums(align_table_mRNA)
		 peak_mRNA<-peakFind(align_table_mRNA, thres=0.40)
		 peak_sRNA<-peakFind2(align_table_sRNA, thres=0.40)
		 
		cons_sRNA<-overlap(peak_sRNA[[1]], as.numeric(tab_aligned[,5]),as.numeric(tab_aligned[,6]),thres=0.6)
		cons_sRNA2<-overlap(peak_sRNA[[2]], as.numeric(tab_aligned[,5]),as.numeric(tab_aligned[,6]),thres=0.6)
		cons_mRNA<-overlap(peak_mRNA, as.numeric(tab_aligned[,1]),as.numeric(tab_aligned[,2]),thres=0.6)
		res<-cbind(tab,cons_mRNA,cons_sRNA,cons_sRNA2)
		p_sub<-rep(NA,nrow(tab))
		cons_sRNA_sub<-rep(NA,nrow(tab))
		cons_sRNA_sub2<-rep(NA,nrow(tab))
		cons_mRNA_sub<-rep(NA,nrow(tab))
		if(nrow(tabsub_aligned)>0){
			
			cons_sRNAsu<-overlap(peak_sRNA[[1]], as.numeric(tabsub_aligned[,5]),as.numeric(tabsub_aligned[,6]),thres=0.6)
			cons_sRNAsu2<-overlap(peak_sRNA[[2]], as.numeric(tabsub_aligned[,5]),as.numeric(tabsub_aligned[,6]),thres=0.6)
			cons_mRNAsu<-overlap(peak_mRNA, as.numeric(tabsub_aligned[,1]),as.numeric(tabsub_aligned[,2]),thres=0.6)
			tabsub<-cbind(tabsub,cons_mRNAsu,cons_sRNAsu)
			
			ov<-na.omit(match(as.character(tabsub[,"name"]),as.character(tab[,"name"])))
			
			p_sub[ov]<-tabsub[,11]
			cons_sRNA_sub[ov]<-cons_sRNAsu
			cons_sRNA_sub2[ov]<-cons_sRNAsu2
			cons_mRNA_sub[ov]<-cons_mRNAsu
			res<-cbind(res,p_sub,cons_mRNA_sub,cons_sRNA_sub,cons_sRNA_sub2)
			
			
		}
		 
		 if(nrow(tabsub_aligned)==0){
			res<-cbind(res,p_sub,cons_mRNA_sub,cons_sRNA_sub)
		 }
		 
		 
		
		write.table(anno_mRNA, file=paste(wd,"/evo_alignments/",i,"_" ,tab[1,9],"/", i,"_" ,tab[1,9], "_mRNA_anno.txt", sep=""), quote=F, row.names=F, col.names=F)
		write.table(anno_sRNA, file=paste(wd,"/evo_alignments/",i,"_" ,tab[1,9],"/", i,"_" ,tab[1,9], "_sRNA_anno.txt", sep=""), quote=F, row.names=F, col.names=F)
		#########################
		 align_anno_mRNA<-("JALVIEW_ANNOTATION")
		 bar_mRNA<-paste(align_table_mRNA, collapse="|")
		 bar_mRNA<-paste("BAR_GRAPH\tweighted interaction region\t", bar_mRNA)
		 cons_site<-rep("E",peak_mRNA[2]-peak_mRNA[1])
		 cons_site<-paste(cons_site, collapse="|")
		 cons_site<-c(rep("|",peak_mRNA[1]-1),"E,consensus interaction site",cons_site,rep("|",length(alignment[[1]])-peak_mRNA[2]))
		 cons_site<-paste(cons_site, collapse="")
		 cons_site<-paste("NO_GRAPH\t \t", cons_site)
		 align_anno_mRNA<-c(align_anno_mRNA, bar_mRNA, cons_site)
		 write.table(align_anno_mRNA, file=paste(wd,"/evo_alignments/",i,"_" ,tab[1,9],"/", i,"_" ,tab[1,9], "_mRNA_annotation.txt", sep=""), quote=F, row.names=F, col.names=F)
		 
		 
		 align_anno_sRNA<-("JALVIEW_ANNOTATION")
		 bar_sRNA<-paste(align_table_sRNA, collapse="|")
		 bar_sRNA<-paste("BAR_GRAPH\tweighted interaction region\t", bar_sRNA)
		 cons_site<-rep("E",peak_sRNA[[1]][2]-peak_sRNA[[1]][1])
		 cons_site<-paste(cons_site, collapse="|")
		 cons_site<-c(rep("|",peak_sRNA[[1]][1]-1),"E,int_site_1",cons_site,rep("|",length(alignment[[1]])-peak_sRNA[[1]][2]))
		 cons_site<-paste(cons_site, collapse="")
		 cons_site<-paste("NO_GRAPH\t \t", cons_site)
		 align_anno_sRNA<-c(align_anno_sRNA, bar_sRNA, cons_site)
		 if(length(peak_sRNA[[2]]>0)){
		 cons_site<-rep("E",peak_sRNA[[2]][2]-peak_sRNA[[2]][1])
		 cons_site<-paste(cons_site, collapse="|")
		 cons_site<-c(rep("|",peak_sRNA[[2]][1]-1),"E,int_site_2",cons_site,rep("|",length(alignment[[1]])-peak_sRNA[[1]][2]))
		 cons_site<-paste(cons_site, collapse="")
		 cons_site<-paste("NO_GRAPH\t \t", cons_site)
		 align_anno_sRNA<-c(align_anno_sRNA, cons_site)
		 }
		 write.table(align_anno_sRNA, file=paste(wd,"/evo_alignments/",i,"_" ,tab[1,9],"/", i,"_" ,tab[1,9], "_sRNA_annotation.txt", sep=""), quote=F, row.names=F, col.names=F)
		
		###############
	
	# #cons_res<-rep(NA,nrow(res))
	# #p_res<-rep(NA,nrow(res))
	
	
	# #cons_res_sub<-rep(NA,nrow(res))
	# #p_res_sub<-rep(NA,nrow(res))
	# for(j in 1:nrow(res)){
	# #print(paste("table",j))
		# #print(j)
		
		# m_opt<-"FALSE"
		# s1_opt<-"FALSE"
		# s2_opt<-"FALSE"
		
		# m_opt_sub<-"FALSE"
		# s1_opt_sub<-"FALSE"
		# s2_opt_sub<-"FALSE"
		# p1<-as.numeric(res[j,"pvalue"])
		
		
		
		# #if(as.numeric(res[j,"pvalue"])<=0.3){
		# #	tem<-tem+1
		# #}
		
		# if((res[j,"cons_mRNA"])=="TRUE"){
		
			
			# m_opt<-"TRUE"
		# }
		 # if((res[j,"cons_sRNA"])=="TRUE"){
			 # s1_opt<-"TRUE"
		 # }
		 # if((res[j,"cons_sRNA2"])=="TRUE"){
			 # s2_opt<-"TRUE"
		 # }
		 
		 
		# if(is.na(res[j,"p_sub"])==F){
		# #if(as.numeric(res[j,"p_sub"])<=0.3){
		# #	tem2<-tem2+1
		
		# p2<-as.numeric(res[j,"p_sub"])
		# #p_res_sub[j]<-p2
		# if((res[j,"cons_mRNA_sub"])=="TRUE"){
			# m_opt_sub<-"TRUE"
		# }
		# if((res[j,"cons_sRNA_sub"])=="TRUE"){
			 # s1_opt_sub<-"TRUE"
		 # }
		 # if((res[j,"cons_sRNA_sub2"])=="TRUE"){
			 # s2_opt_sub<-"TRUE"
		 # }
		 
		 # sub_temp<-paste(p2,m_opt_sub ,s1_opt_sub,s2_opt_sub, sep="|")
		 # cons_res_sub[j]<-sub_temp
		# }
		
		# #p_res[j]<-p1
		
		# opt_temp<-paste(p1,m_opt ,s1_opt,s2_opt, sep="|")
		
		
		# cons_res[j]<-opt_temp
		
		
	# }
	cons_res<-paste(res[,"pvalue"],res[,"cons_mRNA"],res[,"cons_sRNA"],res[,"cons_sRNA2"], sep="|")
	
	res2<-na.omit(res)
	if(nrow(res2)>0){
		cons_res_sub<-paste(res2[,"p_sub"],res2[,"cons_mRNA_sub"],res2[,"cons_sRNA_sub"],res2[,"cons_sRNA_sub2"], sep="|")
		conservation_table_sub[i,na.omit(match(res2[,"orgs"],colnames(conservation_table)))]<-cons_res_sub
	}
	
	write.table(res, file=paste(wd,"/evo_alignments/",i,"_" ,tab[1,9],"/", i,"_" ,tab[1,9], "_mapping_result.txt", sep=""), sep="\t")
	conservation_table[i,na.omit(match(res[,"orgs"],colnames(conservation_table)))]<-cons_res
	
	#p_table[i,na.omit(match(res[,"orgs"],colnames(p_table)))]<-p_res
	#p_table_sub[i,na.omit(match(res[,"orgs"],colnames(p_table)))]<-p_res_sub
	
	########################################
	if(length(exist)>0){
		peak_mRNA<-as.numeric(c(tab_aligned[exist,1], tab_aligned[exist,2]))
		peak_sRNA<-list(as.numeric(c(tab_aligned[exist,5], tab_aligned[exist,6])),NULL)
		 
		cons_sRNA<-overlap(peak_sRNA[[1]], as.numeric(tab_aligned[,5]),as.numeric(tab_aligned[,6]),thres=0.6)
		cons_sRNA2<-overlap(peak_sRNA[[2]], as.numeric(tab_aligned[,5]),as.numeric(tab_aligned[,6]),thres=0.6)
		cons_mRNA<-overlap(peak_mRNA, as.numeric(tab_aligned[,1]),as.numeric(tab_aligned[,2]),thres=0.6)
		res<-cbind(tab,cons_mRNA,cons_sRNA,cons_sRNA2)
		p_sub<-rep(NA,nrow(tab))
		cons_sRNA_sub<-rep(NA,nrow(tab))
		cons_sRNA_sub2<-rep(NA,nrow(tab))
		cons_mRNA_sub<-rep(NA,nrow(tab))
		if(nrow(tabsub_aligned)>0){
			
			cons_sRNAsu<-overlap(peak_sRNA[[1]], as.numeric(tabsub_aligned[,5]),as.numeric(tabsub_aligned[,6]),thres=0.6)
			cons_sRNAsu2<-overlap(peak_sRNA[[2]], as.numeric(tabsub_aligned[,5]),as.numeric(tabsub_aligned[,6]),thres=0.6)
			cons_mRNAsu<-overlap(peak_mRNA, as.numeric(tabsub_aligned[,1]),as.numeric(tabsub_aligned[,2]),thres=0.6)
			tabsub<-cbind(tabsub,cons_mRNAsu,cons_sRNAsu)
			
			ov<-na.omit(match(as.character(tabsub[,"name"]),as.character(tab[,"name"])))
			
			p_sub[ov]<-tabsub[,11]
			cons_sRNA_sub[ov]<-cons_sRNAsu
			cons_sRNA_sub2[ov]<-cons_sRNAsu2
			cons_mRNA_sub[ov]<-cons_mRNAsu
			res<-cbind(res,p_sub,cons_mRNA_sub,cons_sRNA_sub,cons_sRNA_sub2)
			
			
		}
		 
		 if(nrow(tabsub_aligned)==0){
			res<-cbind(res,p_sub,cons_mRNA_sub,cons_sRNA_sub)
		 }
		 
		 
		
		write.table(anno_mRNA, file=paste(wd,"/evo_alignments/",i,"_" ,tab[1,9],"/", i,"_" ,tab[1,9], "_mRNA_anno.txt", sep=""), quote=F, row.names=F, col.names=F)
		write.table(anno_sRNA, file=paste(wd,"/evo_alignments/",i,"_" ,tab[1,9],"/", i,"_" ,tab[1,9], "_sRNA_anno.txt", sep=""), quote=F, row.names=F, col.names=F)
	
	
	# cons_res_ooi<-rep(NA,nrow(res))
	# #p_res<-rep(NA,nrow(res))
	
	
	# cons_res_sub_ooi<-rep(NA,nrow(res))
	# #p_res_sub<-rep(NA,nrow(res))
	# for(j in 1:nrow(res)){
		# #print(j)
		
		# m_opt<-"FALSE"
		# s1_opt<-"FALSE"
		# s2_opt<-"FALSE"
		
		# m_opt_sub<-"FALSE"
		# s1_opt_sub<-"FALSE"
		# s2_opt_sub<-"FALSE"
		# p1<-as.numeric(res[j,"pvalue"])
		
		
		
		# #if(as.numeric(res[j,"pvalue"])<=0.3){
		# #	tem<-tem+1
		# #}
		
		# if((res[j,"cons_mRNA"])=="TRUE"){
		
			
			# m_opt<-"TRUE"
		# }
		 # if((res[j,"cons_sRNA"])=="TRUE"){
			 # s1_opt<-"TRUE"
		 # }
		 # if((res[j,"cons_sRNA2"])=="TRUE"){
			 # s2_opt<-"TRUE"
		 # }
		 
		 
		# if(is.na(res[j,"p_sub"])==F){
		# #if(as.numeric(res[j,"p_sub"])<=0.3){
		# #	tem2<-tem2+1
		
		# p2<-as.numeric(res[j,"p_sub"])
		# #p_res_sub[j]<-p2
		# if((res[j,"cons_mRNA_sub"])=="TRUE"){
			# m_opt_sub<-"TRUE"
		# }
		# if((res[j,"cons_sRNA_sub"])=="TRUE"){
			 # s1_opt_sub<-"TRUE"
		 # }
		 # if((res[j,"cons_sRNA_sub2"])=="TRUE"){
			 # s2_opt_sub<-"TRUE"
		 # }
		 
		 # sub_temp<-paste(p2,m_opt_sub ,s1_opt_sub,s2_opt_sub, sep="|")
		 # cons_res_sub_ooi[j]<-sub_temp
		# }
		
		# #p_res[j]<-p1
		
		# opt_temp<-paste(p1,m_opt ,s1_opt,s2_opt, sep="|")
		
		
		# cons_res_ooi[j]<-opt_temp
		
		
	# }
	
	
	cons_res<-paste(res[,"pvalue"],res[,"cons_mRNA"],res[,"cons_sRNA"],res[,"cons_sRNA2"], sep="|")
	
	res2<-na.omit(res)
	if(nrow(res2)>0){
		cons_res_sub<-paste(res2[,"p_sub"],res2[,"cons_mRNA_sub"],res2[,"cons_sRNA_sub"],res2[,"cons_sRNA_sub2"], sep="|")
		conservation_table_sub_ooi[i,na.omit(match(res2[,"orgs"],colnames(conservation_table)))]<-cons_res_sub
	}
	
	write.table(res, file=paste(wd,"/evo_alignments/",i,"_" ,tab[1,9],"/", i,"_" ,tab[1,9], "_mapping_result.txt", sep=""), sep="\t")
	conservation_table_ooi[i,na.omit(match(res[,"orgs"],colnames(conservation_table)))]<-cons_res
	
	#write.table(res, file=paste(wd,"/evo_alignments/",i,"_" ,tab[1,9],"/", i,"_" ,tab[1,9], "_mapping_result.txt", sep=""), sep="\t")
	# conservation_table_ooi[i,na.omit(match(res[,"orgs"],colnames(conservation_table)))]<-cons_res_ooi
	# conservation_table_sub_ooi[i,na.omit(match(res[,"orgs"],colnames(conservation_table)))]<-cons_res_sub_ooi
		
	}
	###########################
	}
	}
	out<-list(conservation_table,conservation_table_sub,conservation_table_ooi,conservation_table_sub_ooi )
	out
}

conservation_table<-build_anno(ooi=ooi2)
save(conservation_table, file="conservation_table.Rdata")
#################################
























