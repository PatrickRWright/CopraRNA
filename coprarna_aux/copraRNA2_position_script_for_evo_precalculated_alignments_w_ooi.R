# Parameter:
# ooi as refseq id
# script by Jens Georg

# dependency: CopraRNA_available_organisms.txt

# R --slave -f ../copraRNA2_position_script_for_evo_precalculated_alignments_w_ooi.R --args NC_000913

args <- commandArgs(trailingOnly = TRUE) 

ooi2 <- args[1] 



seqle <- function(x,incr=1) { 
  if(!is.numeric(x)) x <- as.numeric(x) 
  n <- length(x)  
  y <- x[-1L] != x[-n] + incr 
  i <- c(which(y|is.na(y)),n) 
  list(lengths = diff(c(0L,i)),
       values = x[head(c(0L,i)+1L,-1L)]) 
} 


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
	
	
	
	th<-thres*s1
	
	mea<-which(interaction_site_distr==s1)
	
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
	peak2<-c()
	sorted2<-sort(unique(interaction_site_distr[-seq(s,e)]), decreasing=T)
	if(length(sorted2)>10){
	s2<-sorted2[1]
	th2<-thres*s2
	eq<-s2/s1
	mea2<-which(interaction_site_distr==s2)
	s2<-mea2[1]
	
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
	
	input<-paste("mafft --maxiterate 1000 --localpair --quiet", " --localpair", " --quiet input_sRNA.fa >", "aligned_sRNA.fa", sep="")
	system(input)
	sRNA_alignment<-read.fasta("aligned_sRNA.fa")
	s_names<-names(sRNA_alignment)
	dir.create("evo_alignments")
	wd<-getwd()
	
	d2<-dir()
	f<-grep("tags.clustered$", d2)
	opt<-read.csv(d2[f], sep=";")
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
	conservation_table_sub<-conservation_table
	conservation_table_ooi<-conservation_table
	conservation_table_sub_ooi<-conservation_table
	conservation_position<-conservation_table
	conservation_position_sub<-conservation_table
	conservation_position_sRNA<-conservation_table
	conservation_position_sRNA_sub<-conservation_table
	consensus_mRNA<-vector("list", nrow(conservation_table))
	consensus_sRNA<-vector("list", nrow(conservation_table))
	
	for(i in 1:nrow(dat)){
		
		#print(i)
		
		e<-grep("Annotation", colnames(dat))
		temp<-dat[i,3:(e-1)]
		temp<-temp[which(temp!="")]
		orgs<-names(temp)
		if(length(temp)>1){
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
		tab<-cbind(tab,paste(genename, "_", query,sep=""),query,opt[pos_opt,"p.value"],orgs,opt[pos_opt,"hybridDP"])
		colnames(tab)<-c("start","end","start_seed", "end_seed","start_sRNA","end_sRNA","start_seed_sRNA", "end_seed_sRNA","name" ,"name2","pvalue","orgs","hybridDP")
		
		pos_sub<-(match(query, (subopt[,1])))
		exist_sub<-which(is.na(pos_sub)==F)
		pos_sub<-na.omit(pos_sub)
		tabsub<-matrix(,1,12)
		if(length(pos_sub)>0){
			tabsub<-subopt[pos_sub,c("start1","end1","seedStart1","seedEnd1","start2","end2","seedStart2","seedEnd2")]
			tabsub<-cbind(tabsub,paste(genename[exist_sub], "_", query[exist_sub],sep=""),query[exist_sub],subopt[pos_sub,"p.value"],orgs[exist_sub],subopt[pos_sub,"hybridDP"])
			colnames(tabsub)<-c("start","end","start_seed", "end_seed","start_sRNA","end_sRNA","start_seed_sRNA", "end_seed_sRNA","name" ,"name2","pvalue","orgs","hybridDP")
		}
	
		tabsub<-na.omit(tabsub)
		pos<-as.numeric(names(sort(table(poslist), decreasing=T))[1])
		alignment<-read.fasta(paste(wd,"/target_alignments/",i, ".aln", sep=""))
		unlink("test2.fa")
		unlink("test.fa")
		dir.create(paste(wd,"/","evo_alignments/",i,"_",tab[1,9],sep=""))
		m_names<-na.omit( match(orgs, names(dat[i,3:(e-1)])))
		s_names2<-na.omit(match(orgs, s_names))
		sRNA_alignment2<-sRNA_alignment[s_names2]
		write.fasta(sRNA_alignment2, paste(nam2[s_names2],(orgs),tab[,10],sep="_"), file.out=paste(wd,"/evo_alignments/",i,"_" ,tab[1,9],"/", i,"_" ,tab[1,9], "_sRNA.fasta", sep=""))
		write.fasta(alignment, paste(nam2[m_names],(orgs),tab[,10],sep="_"), file.out=paste(wd,"/evo_alignments/",i,"_" ,tab[1,9],"/", i,"_" ,tab[1,9], "_mRNA.fasta", sep=""))
		tab_aligned<-tab
		tabsub_aligned<-tabsub
		
		
		 if(nrow(tabsub)>0){
		 ############
		 for(j in 1:nrow(tabsub)){
		 
			align_num_sRNA<-which(names(sRNA_alignment2)==tabsub[j,"orgs"])
			align_num<-which(names(alignment)==tabsub[j,"name2"])
			 nogaps<-which(alignment[[align_num]]!="-")
			 nogaps2<-which(sRNA_alignment2[[align_num_sRNA]]!="-")
			 tabsub_aligned[j,1]<-nogaps[as.numeric(tabsub[j,1])]
			 tabsub_aligned[j,2]<-nogaps[as.numeric(tabsub[j,2])]
			 tabsub_aligned[j,5]<-nogaps2[as.numeric(tabsub[j,5])]
			 tabsub_aligned[j,6]<-nogaps2[as.numeric(tabsub[j,6])]
			 
			 int_postion<-seq(tabsub[j, "start"],tabsub[j, "end"])
			in_interaction<-gregexpr("\\(", tabsub[j,"hybridDP"])[[1]]
			in_interaction<-int_postion[in_interaction]
			in2<-seqle(in_interaction)
			int_list<-c()

			for(jj in 1:length(in2$lengths)){
				temp_int<-c(in2$values[jj],in2$values[jj]+in2$lengths[jj]-1)
				int_list[[jj]]<-temp_int
			}
			
			int_postion_sRNA<-seq(tabsub[j, "start_sRNA"],tabsub[j, "end_sRNA"])
			sRNA_string<-gsub(".*&","",tabsub[j,"hybridDP"])
			sRNA_string<-paste(rev(strsplit(sRNA_string,"")[[1]]), collapse="")
			in_interaction_sRNA<-gregexpr("\\)", sRNA_string)[[1]]
			in_interaction_sRNA<-int_postion_sRNA[in_interaction_sRNA]
			in2_sRNA<-seqle(in_interaction_sRNA)
			int_list_sRNA<-c()

			for(jj in 1:length(in2_sRNA$lengths)){
				temp_int<-c(in2_sRNA$values[jj],in2_sRNA$values[jj]+in2_sRNA$lengths[jj]-1)
				int_list_sRNA[[jj]]<-temp_int
			}
			
			a1<-j
			if(as.numeric(tabsub[j,11])<=0.001){
				for(jj in 1:length(int_list)){
					m<-paste(tabsub[a1,11],"\t", paste(nam2[m_names[a1]],(orgs)[a1],tabsub[a1,10],sep="_"), "\t-1\t", int_list[[jj]][1], "\t", int_list[[jj]][2] ,"\tp<=0001", sep="")
					anno_mRNA<-c(anno_mRNA, m)
					}
					for(jj in 1:length(int_list_sRNA)){
					s<-paste(tabsub[a1,11],"\t", paste(nam2[m_names[a1]],(orgs)[a1],tabsub[a1,10],sep="_"), "\t-1\t", int_list_sRNA[[jj]][1], "\t", int_list_sRNA[[jj]][2] ,"\tp<=0001", sep="")
					anno_sRNA<-c(anno_sRNA, s)
					}
					
					
				
			}
			
			
			if(as.numeric(tabsub[j,11])>0.001 & as.numeric(tabsub[j,11])<=0.01){
				for(jj in 1:length(int_list)){
				m<-paste(tabsub[a1,11],"\t", paste(nam2[m_names[a1]],(orgs)[a1],tabsub[a1,10],sep="_"), "\t-1\t", int_list[[jj]][1], "\t", int_list[[jj]][2] ,"\tp<=001", sep="")
				anno_mRNA<-c(anno_mRNA, m)
				}
				for(jj in 1:length(int_list_sRNA)){
				s<-paste(tabsub[a1,11],"\t", paste(nam2[m_names[a1]],(orgs)[a1],tabsub[a1,10],sep="_"), "\t-1\t", int_list_sRNA[[jj]][1], "\t", int_list_sRNA[[jj]][2] ,"\tp<=001", sep="")
				anno_sRNA<-c(anno_sRNA, s)
				}
				
				
				
			}
			
			
			if(as.numeric(tabsub[j,11])>0.01 & as.numeric(tabsub[j,11])<=0.1){
				for(jj in 1:length(int_list)){
				m<-paste(tabsub[a1,11],"\t", paste(nam2[m_names[a1]],(orgs)[a1],tabsub[a1,10],sep="_"), "\t-1\t", int_list[[jj]][1], "\t", int_list[[jj]][2] ,"\tp<=01", sep="")
				anno_mRNA<-c(anno_mRNA, m)
				}
				for(jj in 1:length(int_list_sRNA)){
				s<-paste(tabsub[a1,11],"\t", paste(nam2[m_names[a1]],(orgs)[a1],tabsub[a1,10],sep="_"), "\t-1\t", int_list_sRNA[[jj]][1], "\t", int_list_sRNA[[jj]][2] ,"\tp<=01", sep="")
				anno_sRNA<-c(anno_sRNA, s)
				}
				
				
				
			}
			
			
			if(as.numeric(tabsub[j,11])>0.1 & as.numeric(tabsub[j,11])<=0.2){
				for(jj in 1:length(int_list)){
				m<-paste(tabsub[a1,11],"\t", paste(nam2[m_names[a1]],(orgs)[a1],tabsub[a1,10],sep="_"), "\t-1\t", int_list[[jj]][1], "\t", int_list[[jj]][2] ,"\tp<=02", sep="")
				anno_mRNA<-c(anno_mRNA, m)
				}
				for(jj in 1:length(int_list_sRNA)){
				s<-paste(tabsub[a1,11],"\t", paste(nam2[m_names[a1]],(orgs)[a1],tabsub[a1,10],sep="_"), "\t-1\t", int_list_sRNA[[jj]][1], "\t", int_list_sRNA[[jj]][2] ,"\tp<=02", sep="")
				anno_sRNA<-c(anno_sRNA, s)
				}
				
				
				
			}
			
			if(as.numeric(tabsub[j,11])>0.2 & as.numeric(tabsub[j,11])<=0.3){
				for(jj in 1:length(int_list)){
					m<-paste(tabsub[a1,11],"\t", paste(nam2[m_names[a1]],(orgs)[a1],tabsub[a1,10],sep="_"), "\t-1\t", int_list[[jj]][1], "\t", int_list[[jj]][2] ,"\tp<=03", sep="")
					anno_mRNA<-c(anno_mRNA, m)
				}
				for(jj in 1:length(int_list_sRNA)){
					s<-paste(tabsub[a1,11],"\t", paste(nam2[m_names[a1]],(orgs)[a1],tabsub[a1,10],sep="_"), "\t-1\t", int_list_sRNA[[jj]][1], "\t", int_list_sRNA[[jj]][2] ,"\tp<=03", sep="")
					anno_sRNA<-c(anno_sRNA, s)
				}
				
				
				
			}
			
			if(as.numeric(tabsub[j,11])>0.3){
				for(jj in 1:length(int_list)){
				m<-paste(tabsub[a1,11],"\t", paste(nam2[m_names[a1]],(orgs)[a1],tabsub[a1,10],sep="_"), "\t-1\t", int_list[[jj]][1], "\t", int_list[[jj]][2] ,"\tp>03", sep="")
				anno_mRNA<-c(anno_mRNA, m)
				}
				for(jj in 1:length(int_list_sRNA)){
				s<-paste(tabsub[a1,11],"\t", paste(nam2[m_names[a1]],(orgs)[a1],tabsub[a1,10],sep="_"), "\t-1\t", int_list_sRNA[[jj]][1], "\t", int_list_sRNA[[jj]][2] ,"\tp>03", sep="")
				anno_sRNA<-c(anno_sRNA, s)
				}
				
				
				
			}
			
			
			}
		 ##################
		 }
		 
		 for(j in 1:nrow(tab)){
			align_num_sRNA<-which(names(sRNA_alignment2)==tab[j,"orgs"])
			align_num<-which(names(alignment)==tab[j,"name2"])
			 nogaps<-which(alignment[[align_num]]!="-")
			 nogaps2<-which(sRNA_alignment2[[align_num_sRNA]]!="-")
			 tab_aligned[j,1]<-nogaps[as.numeric(tab[j,1])]
			 tab_aligned[j,2]<-nogaps[as.numeric(tab[j,2])]
			 tab_aligned[j,5]<-nogaps2[as.numeric(tab[j,5])]
			 tab_aligned[j,6]<-nogaps2[as.numeric(tab[j,6])]
			 
			int_postion<-seq(tab[j, "start"],tab[j, "end"])
			in_interaction<-gregexpr("\\(", tab[j,"hybridDP"])[[1]]
			in_interaction<-int_postion[in_interaction]
			in2<-seqle(in_interaction)
			int_list<-c()

			for(jj in 1:length(in2$lengths)){
				temp_int<-c(in2$values[jj],in2$values[jj]+in2$lengths[jj]-1)
				int_list[[jj]]<-temp_int
			}
			
			int_postion_sRNA<-seq(tab[j, "start_sRNA"],tab[j, "end_sRNA"])
			sRNA_string<-gsub(".*&","",tab[j,"hybridDP"])
			sRNA_string<-paste(rev(strsplit(sRNA_string,"")[[1]]), collapse="")
			in_interaction_sRNA<-gregexpr("\\)", sRNA_string)[[1]]
			in_interaction_sRNA<-int_postion_sRNA[in_interaction_sRNA]
			in2_sRNA<-seqle(in_interaction_sRNA)
			int_list_sRNA<-c()

			for(jj in 1:length(in2_sRNA$lengths)){
				temp_int<-c(in2_sRNA$values[jj],in2_sRNA$values[jj]+in2_sRNA$lengths[jj]-1)
				int_list_sRNA[[jj]]<-temp_int
			}
			
			a1<-j
			if(as.numeric(tab[j,11])<=0.001){
				for(jj in 1:length(int_list)){
					m<-paste(tab[a1,11],"\t", paste(nam2[m_names[a1]],(orgs)[a1],tab[a1,10],sep="_"), "\t-1\t", int_list[[jj]][1], "\t", int_list[[jj]][2] ,"\tp<=0001", sep="")
					anno_mRNA<-c(anno_mRNA, m)
					}
					for(jj in 1:length(int_list_sRNA)){
					s<-paste(tab[a1,11],"\t", paste(nam2[m_names[a1]],(orgs)[a1],tab[a1,10],sep="_"), "\t-1\t", int_list_sRNA[[jj]][1], "\t", int_list_sRNA[[jj]][2] ,"\tp<=0001", sep="")
					anno_sRNA<-c(anno_sRNA, s)
					}
					
					
				
			}
			
			
			if(as.numeric(tab[j,11])>0.001 & as.numeric(tab[j,11])<=0.01){
				for(jj in 1:length(int_list)){
				m<-paste(tab[a1,11],"\t", paste(nam2[m_names[a1]],(orgs)[a1],tab[a1,10],sep="_"), "\t-1\t", int_list[[jj]][1], "\t", int_list[[jj]][2] ,"\tp<=001", sep="")
				anno_mRNA<-c(anno_mRNA, m)
				}
				for(jj in 1:length(int_list_sRNA)){
				s<-paste(tab[a1,11],"\t", paste(nam2[m_names[a1]],(orgs)[a1],tab[a1,10],sep="_"), "\t-1\t", int_list_sRNA[[jj]][1], "\t", int_list_sRNA[[jj]][2] ,"\tp<=001", sep="")
				anno_sRNA<-c(anno_sRNA, s)
				}
				
				
				
			}
			
			
			if(as.numeric(tab[j,11])>0.01 & as.numeric(tab[j,11])<=0.1){
				for(jj in 1:length(int_list)){
				m<-paste(tab[a1,11],"\t", paste(nam2[m_names[a1]],(orgs)[a1],tab[a1,10],sep="_"), "\t-1\t", int_list[[jj]][1], "\t", int_list[[jj]][2] ,"\tp<=01", sep="")
				anno_mRNA<-c(anno_mRNA, m)
				}
				for(jj in 1:length(int_list_sRNA)){
				s<-paste(tab[a1,11],"\t", paste(nam2[m_names[a1]],(orgs)[a1],tab[a1,10],sep="_"), "\t-1\t", int_list_sRNA[[jj]][1], "\t", int_list_sRNA[[jj]][2] ,"\tp<=01", sep="")
				anno_sRNA<-c(anno_sRNA, s)
				}
				
				
				
			}
			
			
			if(as.numeric(tab[j,11])>0.1 & as.numeric(tab[j,11])<=0.2){
				for(jj in 1:length(int_list)){
				m<-paste(tab[a1,11],"\t", paste(nam2[m_names[a1]],(orgs)[a1],tab[a1,10],sep="_"), "\t-1\t", int_list[[jj]][1], "\t", int_list[[jj]][2] ,"\tp<=02", sep="")
				anno_mRNA<-c(anno_mRNA, m)
				}
				for(jj in 1:length(int_list_sRNA)){
				s<-paste(tab[a1,11],"\t", paste(nam2[m_names[a1]],(orgs)[a1],tab[a1,10],sep="_"), "\t-1\t", int_list_sRNA[[jj]][1], "\t", int_list_sRNA[[jj]][2] ,"\tp<=02", sep="")
				anno_sRNA<-c(anno_sRNA, s)
				}
				
				
				
			}
			
			if(as.numeric(tab[j,11])>0.2 & as.numeric(tab[j,11])<=0.3){
				for(jj in 1:length(int_list)){
					m<-paste(tab[a1,11],"\t", paste(nam2[m_names[a1]],(orgs)[a1],tab[a1,10],sep="_"), "\t-1\t", int_list[[jj]][1], "\t", int_list[[jj]][2] ,"\tp<=03", sep="")
					anno_mRNA<-c(anno_mRNA, m)
				}
				for(jj in 1:length(int_list_sRNA)){
					s<-paste(tab[a1,11],"\t", paste(nam2[m_names[a1]],(orgs)[a1],tab[a1,10],sep="_"), "\t-1\t", int_list_sRNA[[jj]][1], "\t", int_list_sRNA[[jj]][2] ,"\tp<=03", sep="")
					anno_sRNA<-c(anno_sRNA, s)
				}
				
				
				
			}
			
			if(as.numeric(tab[j,11])>0.3){
				for(jj in 1:length(int_list)){
				m<-paste(tab[a1,11],"\t", paste(nam2[m_names[a1]],(orgs)[a1],tab[a1,10],sep="_"), "\t-1\t", int_list[[jj]][1], "\t", int_list[[jj]][2] ,"\tp>03", sep="")
				anno_mRNA<-c(anno_mRNA, m)
				}
				for(jj in 1:length(int_list_sRNA)){
				s<-paste(tab[a1,11],"\t", paste(nam2[m_names[a1]],(orgs)[a1],tab[a1,10],sep="_"), "\t-1\t", int_list_sRNA[[jj]][1], "\t", int_list_sRNA[[jj]][2] ,"\tp>03", sep="")
				anno_sRNA<-c(anno_sRNA, s)
				
				}
				
				
				
			}
			
			
			}
			
		write.table(anno_mRNA, file=paste(wd,"/evo_alignments/",i,"_" ,tab[1,9],"/", i,"_" ,tab[1,9], "_mRNA_features.txt", sep=""), quote=F, row.names=F, col.names=F)
		write.table(anno_sRNA, file=paste(wd,"/evo_alignments/",i,"_" ,tab[1,9],"/", i,"_" ,tab[1,9], "_sRNA_features.txt", sep=""), quote=F, row.names=F, col.names=F)
	
			
#MARTIN		write.table(anno_mRNA, file=paste(wd,"/evo_alignments/",i,"_" ,tab[1,9],"/", i,"_" ,tab[1,9], "_mRNA_anno.txt", sep=""), quote=F, row.names=F, col.names=F)
#MARTIN write.table(anno_sRNA, file=paste(wd,"/evo_alignments/",i,"_" ,tab[1,9],"/", i,"_" ,tab[1,9], "_sRNA_anno.txt", sep=""), quote=F, row.names=F, col.names=F)	
		tab_aligned[,11]<-1-as.numeric(tab_aligned[,11])
		tab<-as.matrix(tab)
		 
		 
		 align_table_mRNA<-matrix(,length(alignment[[1]]),(nrow(tab)+nrow(tabsub)))
		 align_table_sRNA<-matrix(,length(sRNA_alignment[[1]]),(nrow(tab)+nrow(tabsub)))
		 align_table_mRNA[]<-0
		 align_table_sRNA[]<-0
		 for(j in 1:nrow(tab_aligned)){
			align_table_mRNA[as.numeric(tab_aligned[j,1]):as.numeric(tab_aligned[j,2]),j]<-as.numeric(tab_aligned[j,11])
			align_table_sRNA[as.numeric(tab_aligned[j,5]):as.numeric(tab_aligned[j,6]),j]<-as.numeric(tab_aligned[j,11])
			
		 }
		 if(nrow(tabsub_aligned)>0){
		 tabsub_aligned[,11]<-1-as.numeric(tabsub_aligned[,11])
		tabsub<-as.matrix(tabsub)
		 for(j in 1:nrow(tabsub_aligned)){
			
			count<-nrow(tab_aligned)
			align_table_mRNA[as.numeric(tabsub_aligned[j,1]):as.numeric(tabsub_aligned[j,2]),j+count]<-as.numeric(tabsub_aligned[j,11])
			align_table_sRNA[as.numeric(tabsub_aligned[j,5]):as.numeric(tabsub_aligned[j,6]),j+count]<-as.numeric(tabsub_aligned[j,11])
			
		 }
		 }
		 
		 
		 align_table_sRNA<-rowSums(align_table_sRNA)
		 align_table_mRNA<-rowSums(align_table_mRNA)
		 peak_mRNA<-peakFind(align_table_mRNA, thres=0.40)
		 peak_sRNA<-peakFind2(align_table_sRNA, thres=0.40)
		 
		consensus_mRNA[[i]]<-peak_mRNA	
		consensus_sRNA[[i]]<-peak_sRNA
			
		names(consensus_mRNA)[i]<-paste(i,"_" ,tab[1,9],  sep="")
		names(consensus_sRNA)[i]<-paste(i,"_" ,tab[1,9],  sep="")
			
		cons_sRNA<-overlap(peak_sRNA[[1]], as.numeric(tab_aligned[,5]),as.numeric(tab_aligned[,6]),thres=0.6)
		cons_sRNA2<-overlap(peak_sRNA[[2]], as.numeric(tab_aligned[,5]),as.numeric(tab_aligned[,6]),thres=0.6)
		cons_mRNA<-overlap(peak_mRNA, as.numeric(tab_aligned[,1]),as.numeric(tab_aligned[,2]),thres=0.6)
		min_mRNA_opt<-min(tab_aligned[which(cons_mRNA==TRUE),1])
		min_sRNA_opt<-min(c(tab_aligned[which(cons_sRNA==TRUE),5],tab_aligned[which(cons_sRNA2==TRUE),5]))
		max_mRNA_opt<-max(tab_aligned[which(cons_mRNA==TRUE),2])
		max_sRNA_opt<-max(c(tab_aligned[which(cons_sRNA==TRUE),6],tab_aligned[which(cons_sRNA2==TRUE),6]))
		
		min_mRNA_sub<-c()
		min_sRNA_sub<-c()
		max_mRNA_sub<-c()
		max_sRNA_sub<-c()
		
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
			
			min_mRNA_sub<-min(tabsub_aligned[which(cons_mRNAsu==TRUE),1])
			min_sRNA_sub<-min(c(tabsub_aligned[which(cons_sRNAsu==TRUE),5],tabsub_aligned[which(cons_sRNAsu2==TRUE),5]))
			max_mRNA_sub<-max(tabsub_aligned[which(cons_mRNAsu==TRUE),2])
			max_sRNA_sub<-max(c(tabsub_aligned[which(cons_sRNAsu==TRUE),6],tabsub_aligned[which(cons_sRNAsu2==TRUE),6]))
			
			
			ov<-na.omit(match(as.character(tabsub[,"name"]),as.character(tab[,"name"])))
			
			p_sub[ov]<-tabsub[,11]
			cons_sRNA_sub[ov]<-cons_sRNAsu
			cons_sRNA_sub2[ov]<-cons_sRNAsu2
			cons_mRNA_sub[ov]<-cons_mRNAsu
			res<-cbind(res,p_sub,cons_mRNA_sub,cons_sRNA_sub,cons_sRNA_sub2)
			
		}
		

		min_cons_int_start_mRNA<-max(0,min(peak_mRNA[1]-5,min_mRNA_opt,min_mRNA_sub))
		max_cons_int_start_mRNA<-min(length(alignment[[1]]),max(peak_mRNA[2]+5,max_mRNA_opt,max_mRNA_sub))
		min_cons_int_start_sRNA<-max(0,min(min(peak_sRNA[[1]][1],peak_sRNA[[2]][1],na.rm=T)-5,min_sRNA_opt,min_sRNA_sub))
		max_cons_int_start_sRNA<-min(length(sRNA_alignment2[[1]]),max(max(peak_sRNA[[1]][2],peak_sRNA[[2]][2],na.rm=T)+5,max_sRNA_opt,max_sRNA_sub))
		
		m1<-vector("list",length(alignment))
		s1<-vector("list",length(sRNA_alignment2))
		
		for(jj in 1:length(alignment)){
			m1[[jj]]<-alignment[[jj]][min_cons_int_start_mRNA:max_cons_int_start_mRNA]
		}
	
		for(jj in 1:length(sRNA_alignment2)){
			s1[[jj]]<-sRNA_alignment2[[jj]][min_cons_int_start_sRNA:max_cons_int_start_sRNA]
		}
		
		s1_rev<-lapply(s1,rev)
		
		conc<-vector("list",length(alignment))
		
		for(jj in 1:length(conc)){
			conc[[jj]]<-c("5\'...",m1[[jj]],"...3\'-&-3\'...",s1_rev[[jj]],"...5\'")
		}
		
		write.fasta(conc, paste(nam2[m_names],(orgs),tab[,10],sep="_"), file.out=paste(wd,"/evo_alignments/",i,"_" ,tab[1,9],"/", i,"_" ,tab[1,9], "_mRNA+sRNA.fasta", sep=""))
		
		anno_mRNA<-c("p<=0001\tD9F69C", "p<=001\tC2F161","p<01\tBAEF4D", "p<02\tB3EE3A","p<03\tA1D634",  "p>03\tFF4500")
		for(j in 1:nrow(tab)){
			
			 
			int_postion2<-seq(tab_aligned[j, "start"],tab_aligned[j, "end"])
			int_postion<-seq(tab[j, "start"],tab[j, "end"])
			if(length(intersect(int_postion2, seq(min_cons_int_start_mRNA,max_cons_int_start_mRNA)))>0){
			
				in_interaction<-gregexpr("\\(", tab[j,"hybridDP"])[[1]]
				in_interaction<-int_postion[in_interaction]
				in2<-seqle(in_interaction)
				int_list<-c()

				for(jj in 1:length(in2$lengths)){
					temp_int<-c(in2$values[jj],in2$values[jj]+in2$lengths[jj]-1)
					int_list[[jj]]<-temp_int
				}
				
				int_postion_sRNA<-seq(tab[j, "start_sRNA"],tab[j, "end_sRNA"])
				sRNA_string<-gsub(".*&","",tab[j,"hybridDP"])
				sRNA_string<-paste((strsplit(sRNA_string,"")[[1]]), collapse="")
				in_interaction_sRNA<-gregexpr("\\)", sRNA_string)[[1]]
				in_interaction_sRNA<-int_postion_sRNA[in_interaction_sRNA]
				in2_sRNA<-seqle(in_interaction_sRNA)
				int_list_sRNA<-c()

				for(jj in 1:length(in2_sRNA$lengths)){
					temp_int<-c(in2_sRNA$values[jj],in2_sRNA$values[jj]+in2_sRNA$lengths[jj]-1)
					int_list_sRNA[[jj]]<-temp_int
				}
				
				le_sRNA<-length(which(sRNA_alignment2[[j]]!="-"))
				le_mRNA<-length(which(alignment[[j]][min_cons_int_start_mRNA:max_cons_int_start_mRNA]!="-"))
				del_mRNA<-length(which(alignment[[j]][1:max((min_cons_int_start_mRNA-1),1)]!="-"))-2
				sRNA_rev<-rev(sRNA_alignment2[[j]])
				del_sRNA<-le_sRNA-length(which(sRNA_alignment2[[j]][1:(max_cons_int_start_sRNA)+1]!="-"))-7-le_mRNA
				
				a1<-j
				if(as.numeric(tab[j,11])<=0.001){
					for(jj in 1:length(int_list)){
						
						m<-paste(tab[a1,11],"\t", paste(nam2[m_names[a1]],(orgs)[a1],tab[a1,10],sep="_"), "\t-1\t", int_list[[jj]][1]-del_mRNA, "\t", min(le_mRNA+4,int_list[[jj]][2]-del_mRNA) ,"\tp<=0001", sep="")
						anno_mRNA<-c(anno_mRNA, m)
						}
						for(jj in 1:length(int_list_sRNA)){
						s<-paste(tab[a1,11],"\t", paste(nam2[m_names[a1]],(orgs)[a1],tab[a1,10],sep="_"), "\t-1\t", max(le_mRNA+6,le_sRNA-int_list_sRNA[[jj]][2]-del_sRNA), "\t", le_sRNA-int_list_sRNA[[jj]][1]-del_sRNA ,"\tp<=0001", sep="")
						anno_mRNA<-c(anno_mRNA, s)
						}
						
						
					
				}
				
				
				if(as.numeric(tab[j,11])>0.001 & as.numeric(tab[j,11])<=0.01){
					for(jj in 1:length(int_list)){
					m<-paste(tab[a1,11],"\t", paste(nam2[m_names[a1]],(orgs)[a1],tab[a1,10],sep="_"), "\t-1\t", int_list[[jj]][1]-del_mRNA, "\t", min(le_mRNA+4,int_list[[jj]][2]-del_mRNA) ,"\tp<=001", sep="")
					anno_mRNA<-c(anno_mRNA, m)
					}
					for(jj in 1:length(int_list_sRNA)){
					s<-paste(tab[a1,11],"\t", paste(nam2[m_names[a1]],(orgs)[a1],tab[a1,10],sep="_"), "\t-1\t", max(le_mRNA+6,le_sRNA-int_list_sRNA[[jj]][2]-del_sRNA), "\t", le_sRNA-int_list_sRNA[[jj]][1]-del_sRNA ,"\tp<=001", sep="")
					anno_mRNA<-c(anno_mRNA, s)
					}
					
					
					
				}
				
				
				if(as.numeric(tab[j,11])>0.01 & as.numeric(tab[j,11])<=0.1){
					for(jj in 1:length(int_list)){
					m<-paste(tab[a1,11],"\t", paste(nam2[m_names[a1]],(orgs)[a1],tab[a1,10],sep="_"), "\t-1\t", int_list[[jj]][1]-del_mRNA, "\t", min(le_mRNA+4,int_list[[jj]][2]-del_mRNA) ,"\tp<=01", sep="")
					anno_mRNA<-c(anno_mRNA, m)
					}
					for(jj in 1:length(int_list_sRNA)){
					s<-paste(tab[a1,11],"\t", paste(nam2[m_names[a1]],(orgs)[a1],tab[a1,10],sep="_"), "\t-1\t", max(le_mRNA+6,le_sRNA-int_list_sRNA[[jj]][2]-del_sRNA), "\t", le_sRNA-int_list_sRNA[[jj]][1]-del_sRNA,"\tp<=01", sep="")
					anno_mRNA<-c(anno_mRNA, s)
					}
					
					
					
				}
				
				
				if(as.numeric(tab[j,11])>0.1 & as.numeric(tab[j,11])<=0.2){
					for(jj in 1:length(int_list)){
					m<-paste(tab[a1,11],"\t", paste(nam2[m_names[a1]],(orgs)[a1],tab[a1,10],sep="_"), "\t-1\t", int_list[[jj]][1]-del_mRNA, "\t", min(le_mRNA+4,int_list[[jj]][2]-del_mRNA) ,"\tp<=02", sep="")
					anno_mRNA<-c(anno_mRNA, m)
					}
					for(jj in 1:length(int_list_sRNA)){
					s<-paste(tab[a1,11],"\t", paste(nam2[m_names[a1]],(orgs)[a1],tab[a1,10],sep="_"), "\t-1\t", max(le_mRNA+6,le_sRNA-int_list_sRNA[[jj]][2]-del_sRNA), "\t", le_sRNA-int_list_sRNA[[jj]][1]-del_sRNA ,"\tp<=02", sep="")
					anno_mRNA<-c(anno_mRNA, s)
					}
					
					
					
				}
				
				if(as.numeric(tab[j,11])>0.2 & as.numeric(tab[j,11])<=0.3){
					for(jj in 1:length(int_list)){
						m<-paste(tab[a1,11],"\t", paste(nam2[m_names[a1]],(orgs)[a1],tab[a1,10],sep="_"), "\t-1\t", int_list[[jj]][1]-del_mRNA, "\t", min(le_mRNA+4,int_list[[jj]][2]-del_mRNA) ,"\tp<=03", sep="")
						anno_mRNA<-c(anno_mRNA, m)
					}
					for(jj in 1:length(int_list_sRNA)){
						s<-paste(tab[a1,11],"\t", paste(nam2[m_names[a1]],(orgs)[a1],tab[a1,10],sep="_"), "\t-1\t", max(le_mRNA+6,le_sRNA-int_list_sRNA[[jj]][2]-del_sRNA), "\t", le_sRNA-int_list_sRNA[[jj]][1]-del_sRNA ,"\tp<=03", sep="")
						anno_mRNA<-c(anno_mRNA, s)
					}
					
					
					
				}
				
				if(as.numeric(tab[j,11])>0.3){
					for(jj in 1:length(int_list)){
					m<-paste(tab[a1,11],"\t", paste(nam2[m_names[a1]],(orgs)[a1],tab[a1,10],sep="_"), "\t-1\t", int_list[[jj]][1]-del_mRNA, "\t", min(le_mRNA+4,int_list[[jj]][2]-del_mRNA) ,"\tp>03", sep="")
					anno_mRNA<-c(anno_mRNA, m)
					}
					for(jj in 1:length(int_list_sRNA)){
					s<-paste(tab[a1,11],"\t", paste(nam2[m_names[a1]],(orgs)[a1],tab[a1,10],sep="_"), "\t-1\t", max(le_mRNA+6,le_sRNA-int_list_sRNA[[jj]][2]-del_sRNA), "\t", le_sRNA-int_list_sRNA[[jj]][1]-del_sRNA ,"\tp>03", sep="")
					anno_mRNA<-c(anno_mRNA, s)
					
					}
					
					
					
				}
				
			}
			}
			
###########################################################			
			if(nrow(tabsub)>0){
			for(j in 1:nrow(tabsub)){	
			int_postion2<-seq(tabsub_aligned[j, "start"],tabsub_aligned[j, "end"])
			int_postion<-seq(tabsub[j, "start"],tabsub[j, "end"])
			if(length(intersect(int_postion2, seq(min_cons_int_start_mRNA,max_cons_int_start_mRNA)))>0){
			
				in_interaction<-gregexpr("\\(", tabsub[j,"hybridDP"])[[1]]
				in_interaction<-int_postion[in_interaction]
				in2<-seqle(in_interaction)
				int_list<-c()

				for(jj in 1:length(in2$lengths)){
					temp_int<-c(in2$values[jj],in2$values[jj]+in2$lengths[jj]-1)
					int_list[[jj]]<-temp_int
				}
				
				int_postion_sRNA<-seq(tabsub[j, "start_sRNA"],tabsub[j, "end_sRNA"])
				sRNA_string<-gsub(".*&","",tabsub[j,"hybridDP"])
				sRNA_string<-paste((strsplit(sRNA_string,"")[[1]]), collapse="")
				in_interaction_sRNA<-gregexpr("\\)", sRNA_string)[[1]]
				in_interaction_sRNA<-int_postion_sRNA[in_interaction_sRNA]
				in2_sRNA<-seqle(in_interaction_sRNA)
				int_list_sRNA<-c()

				for(jj in 1:length(in2_sRNA$lengths)){
					temp_int<-c(in2_sRNA$values[jj],in2_sRNA$values[jj]+in2_sRNA$lengths[jj]-1)
					int_list_sRNA[[jj]]<-temp_int
				}
				
				le_sRNA<-length(which(sRNA_alignment2[[j]]!="-"))
				le_mRNA<-length(which(alignment[[j]][min_cons_int_start_mRNA:max_cons_int_start_mRNA]!="-"))
				del_mRNA<-length(which(alignment[[j]][1:max(min_cons_int_start_mRNA-1,1)]!="-"))-2
				sRNA_rev<-rev(sRNA_alignment2[[j]])
				del_sRNA<-le_sRNA-length(which(sRNA_alignment2[[j]][1:(max_cons_int_start_sRNA)+1]!="-"))-7-le_mRNA
				
				a1<-j
				if(as.numeric(tabsub[j,11])<=0.001){
					for(jj in 1:length(int_list)){
						
						m<-paste(tabsub[a1,11],"\t", paste(nam2[m_names[a1]],(orgs)[a1],tabsub[a1,10],sep="_"), "\t-1\t", int_list[[jj]][1]-del_mRNA, "\t", min(le_mRNA+4,int_list[[jj]][2]-del_mRNA ),"\tp<=0001", sep="")
						anno_mRNA<-c(anno_mRNA, m)
						}
						for(jj in 1:length(int_list_sRNA)){
						s<-paste(tabsub[a1,11],"\t", paste(nam2[m_names[a1]],(orgs)[a1],tabsub[a1,10],sep="_"), "\t-1\t", max(le_mRNA+6,le_sRNA-int_list_sRNA[[jj]][2]-del_sRNA), "\t", le_sRNA-int_list_sRNA[[jj]][1]-del_sRNA ,"\tp<=0001", sep="")
						anno_mRNA<-c(anno_mRNA, s)
						}
						
						
					
				}
				
				
				if(as.numeric(tabsub[j,11])>0.001 & as.numeric(tabsub[j,11])<=0.01){
					for(jj in 1:length(int_list)){
					m<-paste(tabsub[a1,11],"\t", paste(nam2[m_names[a1]],(orgs)[a1],tabsub[a1,10],sep="_"), "\t-1\t", int_list[[jj]][1]-del_mRNA, "\t", min(le_mRNA+4,int_list[[jj]][2]-del_mRNA) ,"\tp<=001", sep="")
					anno_mRNA<-c(anno_mRNA, m)
					}
					for(jj in 1:length(int_list_sRNA)){
					s<-paste(tabsub[a1,11],"\t", paste(nam2[m_names[a1]],(orgs)[a1],tabsub[a1,10],sep="_"), "\t-1\t", max(le_mRNA+6,le_sRNA-int_list_sRNA[[jj]][2]-del_sRNA), "\t", le_sRNA-int_list_sRNA[[jj]][1]-del_sRNA ,"\tp<=001", sep="")
					anno_mRNA<-c(anno_mRNA, s)
					}
					
					
					
				}
				
				
				if(as.numeric(tabsub[j,11])>0.01 & as.numeric(tabsub[j,11])<=0.1){
					for(jj in 1:length(int_list)){
					m<-paste(tabsub[a1,11],"\t", paste(nam2[m_names[a1]],(orgs)[a1],tabsub[a1,10],sep="_"), "\t-1\t", int_list[[jj]][1]-del_mRNA, "\t", min(le_mRNA+4,int_list[[jj]][2]-del_mRNA) ,"\tp<=01", sep="")
					anno_mRNA<-c(anno_mRNA, m)
					}
					for(jj in 1:length(int_list_sRNA)){
					s<-paste(tabsub[a1,11],"\t", paste(nam2[m_names[a1]],(orgs)[a1],tabsub[a1,10],sep="_"), "\t-1\t", max(le_mRNA+6,le_sRNA-int_list_sRNA[[jj]][2]-del_sRNA), "\t", le_sRNA-int_list_sRNA[[jj]][1]-del_sRNA,"\tp<=01", sep="")
					anno_mRNA<-c(anno_mRNA, s)
					}
					
					
					
				}
				
				
				if(as.numeric(tabsub[j,11])>0.1 & as.numeric(tabsub[j,11])<=0.2){
					for(jj in 1:length(int_list)){
					m<-paste(tabsub[a1,11],"\t", paste(nam2[m_names[a1]],(orgs)[a1],tabsub[a1,10],sep="_"), "\t-1\t", int_list[[jj]][1]-del_mRNA, "\t", min(le_mRNA+4,int_list[[jj]][2]-del_mRNA) ,"\tp<=02", sep="")
					anno_mRNA<-c(anno_mRNA, m)
					}
					for(jj in 1:length(int_list_sRNA)){
					s<-paste(tabsub[a1,11],"\t", paste(nam2[m_names[a1]],(orgs)[a1],tabsub[a1,10],sep="_"), "\t-1\t", max(le_mRNA+6,le_sRNA-int_list_sRNA[[jj]][2]-del_sRNA), "\t", le_sRNA-int_list_sRNA[[jj]][1]-del_sRNA ,"\tp<=02", sep="")
					anno_mRNA<-c(anno_mRNA, s)
					}
					
					
					
				}
				
				if(as.numeric(tabsub[j,11])>0.2 & as.numeric(tabsub[j,11])<=0.3){
					for(jj in 1:length(int_list)){
						m<-paste(tabsub[a1,11],"\t", paste(nam2[m_names[a1]],(orgs)[a1],tabsub[a1,10],sep="_"), "\t-1\t", int_list[[jj]][1]-del_mRNA, "\t", min(le_mRNA+4,int_list[[jj]][2]-del_mRNA) ,"\tp<=03", sep="")
						anno_mRNA<-c(anno_mRNA, m)
					}
					for(jj in 1:length(int_list_sRNA)){
						s<-paste(tabsub[a1,11],"\t", paste(nam2[m_names[a1]],(orgs)[a1],tabsub[a1,10],sep="_"), "\t-1\t", max(le_mRNA+6,le_sRNA-int_list_sRNA[[jj]][2]-del_sRNA), "\t", le_sRNA-int_list_sRNA[[jj]][1]-del_sRNA ,"\tp<=03", sep="")
						anno_mRNA<-c(anno_mRNA, s)
					}
					
					
					
				}
				
				if(as.numeric(tabsub[j,11])>0.3){
					for(jj in 1:length(int_list)){
					m<-paste(tabsub[a1,11],"\t", paste(nam2[m_names[a1]],(orgs)[a1],tabsub[a1,10],sep="_"), "\t-1\t", int_list[[jj]][1]-del_mRNA, "\t", min(le_mRNA+4,int_list[[jj]][2]-del_mRNA) ,"\tp>03", sep="")
					anno_mRNA<-c(anno_mRNA, m)
					}
					for(jj in 1:length(int_list_sRNA)){
					s<-paste(tabsub[a1,11],"\t", paste(nam2[m_names[a1]],(orgs)[a1],tabsub[a1,10],sep="_"), "\t-1\t", max(le_mRNA+6,le_sRNA-int_list_sRNA[[jj]][2]-del_sRNA), "\t", le_sRNA-int_list_sRNA[[jj]][1]-del_sRNA ,"\tp>03", sep="")
					anno_mRNA<-c(anno_mRNA, s)
					
					}
					
					
					
				}
				
			}
			}
			}
		write.table(anno_mRNA, file=paste(wd,"/evo_alignments/",i,"_" ,tab[1,9],"/", i,"_" ,tab[1,9], "_mRNA+sRNA_anno.txt", sep=""), quote=F, row.names=F, col.names=F)
			
		
		 if(nrow(tabsub_aligned)==0){
			res<-cbind(res,p_sub,cons_mRNA_sub,cons_sRNA_sub)
		 }
		 
		 
		
		
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
		 cons_site<-c(rep("|",peak_sRNA[[1]][1]-1),"E,int_site_1",cons_site,rep("|",length(sRNA_alignment[[1]])-peak_sRNA[[1]][2]))
		 cons_site<-paste(cons_site, collapse="")
		 cons_site<-paste("NO_GRAPH\t \t", cons_site)
		 align_anno_sRNA<-c(align_anno_sRNA, bar_sRNA, cons_site)
		 if(length(peak_sRNA[[2]]>0)){
		 cons_site<-rep("E",peak_sRNA[[2]][2]-peak_sRNA[[2]][1])
		 cons_site<-paste(cons_site, collapse="|")
		 cons_site<-c(rep("|",peak_sRNA[[2]][1]-1),"E,int_site_2",cons_site,rep("|",length(sRNA_alignment[[1]])-peak_sRNA[[1]][2]))
		 cons_site<-paste(cons_site, collapse="")
		 cons_site<-paste("NO_GRAPH\t \t", cons_site)
		 align_anno_sRNA<-c(align_anno_sRNA, cons_site)
		 }
		 write.table(align_anno_sRNA, file=paste(wd,"/evo_alignments/",i,"_" ,tab[1,9],"/", i,"_" ,tab[1,9], "_sRNA_annotation.txt", sep=""), quote=F, row.names=F, col.names=F)
		
		###############
	cons_res<-paste(res[,"pvalue"],res[,"cons_mRNA"],res[,"cons_sRNA"],res[,"cons_sRNA2"], sep="|")
	
	res2<-na.omit(res)
	if(nrow(res2)>0){
		cons_res_sub<-paste(res2[,"p_sub"],res2[,"cons_mRNA_sub"],res2[,"cons_sRNA_sub"],res2[,"cons_sRNA_sub2"], sep="|")
		conservation_table_sub[i,na.omit(match(res2[,"orgs"],colnames(conservation_table)))]<-cons_res_sub
		res2_position<-paste(as.numeric(tabsub_aligned[,1]),as.numeric(tabsub_aligned[,2]), sep="|")
		res2_position_sRNA<-paste(as.numeric(tabsub_aligned[,5]),as.numeric(tabsub_aligned[,6]), sep="|")
		conservation_position_sub[i,na.omit(match(res2[,"orgs"],colnames(conservation_table)))]<-res2_position
		conservation_position_sRNA_sub[i,na.omit(match(res2[,"orgs"],colnames(conservation_table)))]<-res2_position_sRNA
	}
	
	write.table(res, file=paste(wd,"/evo_alignments/",i,"_" ,tab[1,9],"/", i,"_" ,tab[1,9], "_mapping_result.txt", sep=""), sep="\t")
	conservation_table[i,na.omit(match(res[,"orgs"],colnames(conservation_table)))]<-cons_res
	res_position<-paste(as.numeric(tab_aligned[,1]),as.numeric(tab_aligned[,2]), sep="|")
	res_position_sRNA<-paste(as.numeric(tab_aligned[,5]),as.numeric(tab_aligned[,6]), sep="|")
	conservation_position[i,na.omit(match(res[,"orgs"],colnames(conservation_table)))]<-res_position
	conservation_position_sRNA[i,na.omit(match(res[,"orgs"],colnames(conservation_table)))]<-res_position_sRNA
			
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
		 
		 
		
		#write.table(anno_mRNA, file=paste(wd,"/evo_alignments/",i,"_" ,tab[1,9],"/", i,"_" ,tab[1,9], "_mRNA
   # xt", sep=""), quote=F, row.names=F, col.names=F)
		#write.table(anno_sRNA, file=paste(wd,"/evo_alignments/",i,"_" ,tab[1,9],"/", i,"_" ,tab[1,9], "_sRNA_anno.txt", sep=""), quote=F, row.names=F, col.names=F)
	
	cons_res<-paste(res[,"pvalue"],res[,"cons_mRNA"],res[,"cons_sRNA"],res[,"cons_sRNA2"], sep="|")
	
	res2<-na.omit(res)
	if(nrow(res2)>0){
		cons_res_sub<-paste(res2[,"p_sub"],res2[,"cons_mRNA_sub"],res2[,"cons_sRNA_sub"],res2[,"cons_sRNA_sub2"], sep="|")
		conservation_table_sub_ooi[i,na.omit(match(res2[,"orgs"],colnames(conservation_table)))]<-cons_res_sub
		
	}
	
	write.table(res, file=paste(wd,"/evo_alignments/",i,"_" ,tab[1,9],"/", i,"_" ,tab[1,9], "_mapping_result.txt", sep=""), sep="\t")
	conservation_table_ooi[i,na.omit(match(res[,"orgs"],colnames(conservation_table)))]<-cons_res
	
	
	}
	###########################
	}
	}
	consensus_both<-list(consensus_mRNA,consensus_sRNA)
	interaction_positions<-list(conservation_position,conservation_position_sub,conservation_position_sRNA,conservation_position_sRNA_sub)
	save(consensus_both, file="consensus_positions.Rdata")
	save(interaction_positions, file="interaction_positions.Rdata")
	out<-list(conservation_table,conservation_table_sub,conservation_table_ooi,conservation_table_sub_ooi )
	out
}

conservation_table<-build_anno(ooi=ooi2)
save(conservation_table, file="conservation_table.Rdata")
#################################
