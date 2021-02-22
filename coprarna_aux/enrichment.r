require(genbankr)
require(topGO)
require(tidyr)
require(GOSemSim)
require(AnnotationForge)
require(genbankr)
#load("enrichment.Rdata")
#load("16S_rooted_tree.Rdata")
#weight<-read.csv("weights.txt", sep="\t")

#require(ggplot2)
#require(ggtree)

#R --slave -f ~/media/jens@margarita/CopraRNA-git/coprarna_aux/enrichment.r --args num=200 genomes_path=~/media/jens@margarita/For_CopraRNA2.0/Genomes/  

# get absolute path  
#path="/home/jens/CopraRNA-git/coprarna_aux/"
#path="~/media/jens@margarita/CopraRNA-git/coprarna_aux/"
#genomes_path="~/media/jens@margarita/For_CopraRNA2.0/Genomes/"
genomes_path=getwd()
genomes_path<-paste0(genomes_path,"/")
initial.options <- commandArgs(trailingOnly = FALSE)
path<-initial.options [4]
path<-sub("functional_enrichment.r","",path)
print(path)

load("copra_results_all.Rdata")
all_orgs<-names(copra_results)

num=200   # number of top CopraRNA predictions used for the enrichment


# transforming arguments into valid variables 
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

num<-as.numeric(num)


# preset path to required files, path can also be specified as argument
copref_path<-paste(path,"coprarana-aux/CopraRNA_available_organisms.txt",sep="")  ########### fix path


enrich_list<-vector("list",length(all_orgs))
names(enrich_list)<-all_orgs

for(i in 1:1){
print(i)
cop_all<-as.matrix(copra_results[[i]])
present<-which(cop_all[,all_orgs[i]]!="")
cop_all<-cop_all[present,]
e<-grep("Annotation", colnames(cop_all))

cop<-cop_all[,4:(e-1)]

# int<-cop[,all_orgs[i]]
# int<-strsplit(int, "\\|")
# int<-as.numeric(do.call(rbind,int)[,3])
# ab<-which(int>=int_thres)
# if(length(ab)>0){
	# cop<-cop[-ab,]
# }
gbk<-paste(genomes_path,all_orgs[i],".gb.gz",sep="")
gbk1<-readLines(gbk)
org<-grep("ORGANISM",gbk1)
org<-gbk1[org]
org<-gsub(".*ORGANISM. ","",org)
org<-gsub(" .*","",org)
gb<-readGenBank(gbk)

cd<-cds(gb)
cd<-GenomicRanges::mcols(cd)

url <- "https://www.uniprot.org/uploadlists/"

cop_ids<-unname(cop[,all_orgs[i]])
cop_ids<-gsub("\\(.*","",cop_ids)


id_pos<-na.omit(match(tolower(cop_ids), tolower(cd$locus_tag)))

ids<-unlist(cd$locus_tag)[id_pos]
ids_keep<-ids


id_map<-cbind(ids,ids)
if(all_orgs[i]=="NC_002505"){
	ids<-gsub("VC","VC_",ids)
	id_map[,1]<-ids
}
query_string=paste(na.omit(ids), collapse=" ")
#get go_terms, uniprot id, go_ids
params = list(
  from = "GENENAME",
  to = "ACC",
  columns = "go,go-id,database(KEGG),id,organism,organism-id,database(InterPro),comment(PATHWAY),genes(PREFERRED),genes(ALTERNATIVE),genes,protein names",
  format = "tab",
  query = query_string
)

r <- httr::POST(url, body = params, encode = "form")
anno<-httr::content(r,encoding ="UTF-8")
if(is.null(anno)){
	ids<-unlist(cd$old_locus_tag)[id_pos]	
	id_map<-cbind(ids,unlist(cd$locus_tag)[id_pos])
	query_string=paste(na.omit(ids), collapse=" ")
	#get go_terms, uniprot id, go_ids
	params = list(
	  from = "GENENAME",
	  to = "ACC",
	 
	  columns = "go,go-id,database(KEGG),id,organism,organism-id,database(InterPro),comment(PATHWAY),genes(PREFERRED),genes(ALTERNATIVE),genes,protein names",
	  format = "tab",
	  query = query_string
	)

	r <- httr::POST(url, body = params, encode = "form")
	anno<-httr::content(r,encoding ="UTF-8")
}


if(is.null(anno)){
	j<-1
	r2<-"\n"
	while(r2 == "\n" & j <= 10){
		r <- httr::GET(paste("http://rest.kegg.jp/find/genes/",na.omit(unlist(cd$locus_tag))[j] ,sep=""))
		r2<-httr::content(r,encoding ="UTF-8")
		j<-j+1
		if(is.null(r2)){
			r2<-"\n"
		}
	}
	
	if(r2 == "\n"){
		
		j<-1
		while(r2 == "\n" & j <= 10){
		#print(r2)
			r <- httr::GET(paste("http://rest.kegg.jp/find/genes/",na.omit(unlist(cd$old_locus_tag))[j] ,sep=""))
			r2<-httr::content(r,encoding ="UTF-8")
			j<-j+1
			if(is.null(r2)){
				r2<-"\n"
			}
		}	
		if(r2 != "\n"){
			r2<-gsub(":.*","",r2)
			ids<-unlist(cd$old_locus_tag)[id_pos]	
			ids<-paste(r2,":",ids,sep="")
			id_map<-cbind(ids,unlist(cd$locus_tag)[id_pos])
			query_string=paste(na.omit(ids), collapse=" ")
			#get go_terms, uniprot id, go_ids
			params = list(
			  from = "KEGG_ID",
			  to = "ACC",
			   columns = "go,go-id,database(KEGG),id,organism,organism-id,database(InterPro),comment(PATHWAY),genes(PREFERRED),genes(ALTERNATIVE),genes,protein names",
			  format = "tab",
			  query = query_string
			)

			r <- httr::POST(url, body = params, encode = "form")
			anno<-httr::content(r,encoding ="UTF-8")
		} else {
			next
		}
	} else {
		r2<-gsub(":.*","",r2)
		ids<-unlist(cd$locus_tag)[id_pos]	
		ids<-paste(r2,":",ids,sep="")
		id_map<-cbind(ids,unlist(cd$locus_tag)[id_pos])
		query_string=paste(na.omit(ids), collapse=" ")
		#get go_terms, uniprot id, go_ids
		params = list(
		  from = "KEGG_ID",
		  to = "ACC",
		   columns = "go,go-id,database(KEGG),id,organism,organism-id,database(InterPro),comment(PATHWAY),genes(PREFERRED),genes(ALTERNATIVE),genes,protein names",
		  format = "tab",
		  query = query_string
		)

		r <- httr::POST(url, body = params, encode = "form")
		anno<-httr::content(r,encoding ="UTF-8")
	}
}

if(is.null(anno)==F ){
write.table(anno,file="tmp_uniprot.txt",row.names=F, col.names=F, quote=F)
anno<-read.csv("tmp_uniprot.txt", header=F, skip=1, sep="\t")
if(is.null(anno)==F & nrow(anno) > 100){
path<-gsub(",","",anno[,7])
path<-gsub(" ","_",path)
path<-gsub("PATHWAY:_","",path)

term<-anno[,2]
term<-gsub(" ","",term)
term<-gsub(";;","",term)
term<-gsub(";",", ",term)


geneid<-anno[,6]

mat<-match(as.character(anno[,ncol(anno)]), id_map[,1])
an<-id_map[mat,2]

coprarank<-match(toupper(an),gsub("\\(.*","",toupper(unlist(cop[,all_orgs[i]]))))
term<-cbind(toupper(an),term)
geneid<-cbind(toupper(an),geneid)

term2<-cbind(anno,id_map[mat,],coprarank)
write.table(term, file="enrich_anno.txt", sep="\t", quote=F, row.names=F,col.names=F)
geneID2GO <- readMappings("enrich_anno.txt")
geneNames <- names(geneID2GO)
#myInterestingGenes<-cop_ids[1:num]
#geneList<-as.numeric(cop_all[,2])
geneList<-1:nrow(cop_all)
names(geneList)<-toupper(cop_ids)


geneList<-na.omit(geneList)
#geneList<-na.omit(geneList)

#geneNames <- names(geneID2GO)
#myInterestingGenes <- geneNames[1:200]
#geneList <- factor(as.integer(geneNames %in% myInterestingGenes))
#names(geneList) <- geneNames


topDiffGenes <- function(allScore) { 
	#return(allScore < 0.01) 
	return(allScore < (num+1))
}
genesOfInterest<-topDiffGenes(geneList)

# topDiffGenes2 <- function(allScore) { 
	# return(allScore < 0.2) 
	# #return(allScore < 201)
# }
# genesOfInterest2<-topDiffGenes2(geneList)


#geneList <- factor(as.integer(geneNames %in% myInterestingGenes))
#names(geneList) <- geneNames
go_cat<-c("BP","MF","CC")
temp_enrich<-list()
temp_go_dat<-list()
nodes<-c()
for(jj in go_cat){
GOdata <- new("topGOdata",ontology = jj,allGenes = geneList,geneSel = topDiffGenes,annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = 2)
temp_go_dat[[jj]]<-GOdata
#result01FS <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")
#result01ks <- runTest(GOdata, algorithm = "weight01", statistic = "ks")
classicFS <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
#classicKS <- runTest(GOdata, algorithm = "classic", statistic = "ks")
#nodes<-c(nodes,result01ks@geneData[4])
#nodes<-c(nodes,result01FS@geneData[4])
nodes<-c(nodes,classicFS@geneData[4])
#nodes<-c(nodes,classicKS@geneData[4])
#GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList,annot = annFUN.gene2GO, gene2GO = geneID2GO)


allRes <- GenTable(GOdata,classicFS=classicFS ,ranksOf = "classicFS", topNodes = min(nodes)) # classicKS = classicKS , weight01FS=result01FS, weight01ks=result01ks


allRes$genes <- sapply(allRes$GO.ID, function(x)
{
  genes<-genesInTerm(GOdata, x) 
  paste(genes[[1]][genes[[1]] %in% names(which(genesOfInterest))],collapse=",")
})

allRes$ranks <- sapply(allRes$genes, function(x)
{
	if(length(x)>0){
		x<-strsplit(x,",")[[1]]
		posi<-paste(match(x, toupper(gsub("\\(.*","",cop[,all_orgs[i]]))), collapse=",")
		posi
  }
})


#allRes <- GenTable(GOdata,   ks=resultKs,  topNodes = 1000)
#write.table(allRes, file="enrichment.txt", sep="\t")
temp_enrich[[jj]]<-allRes
}

enrich_list[[i]][[1]]<-temp_enrich
enrich_list[[i]][[2]]<-temp_go_dat
enrich_list[[i]][[3]]<-term2
enrich_list[[i]][[4]]<-anno
enrich_list[[i]][[5]]<-cd
}
}
}
save (enrich_list, file="enrichment.Rdata")





# cd<-cds(gb)
# cd<-GenomicRanges::mcols(cd)





# parameters for term reduction
thres<-0.05		# threshold for enrichment p-Value
max_num=500		# maximal number of members for a term. Removes very broad general terms
min_num=2		# minimal number of members for a term. Renmoves very narrow terms
min_num2=2 		# minimal number of significant genes in a term.
cutoff<-0.6		# overlap threshold for term condensation based on geometric mean of best enrichment p-Value, information content of term and number of significant genes in term
thres_ov<-0.6	# overlap threshold for the number of shared genes between two terms. Proportion of the term with less members.



cd<-cds(gb)
cd<-GenomicRanges::mcols(cd)
ooi<-readLines("ncrna.fa")
ooi<-gsub(">ncRNA_","",ooi[1])
anno<-enrich_list[[ooi]][[4]]


pack<-require(paste0("org.",gsub("_","",ooi),".eg.db"),character.only=TRUE)

if(pack==F){

	fSym<-cbind(seq(1:nrow(cd)),seq(1:nrow(cd)),toupper(as.character(cd[,"locus_tag"])),as.character(cd[,"gene"]))
	colnames(fSym) <- c("GID","ENTREZID","SYMBOL","GENENAME")
	
	dup<-which(duplicated(unlist(fSym[,3])))
	if(length(dup)>0){
		fSym<-fSym[-dup,]
	}
	na<-which(is.na(unlist(fSym[,1])))
	if(length(na)>0){
		fSym<-fSym[-na,]
	}
	na<-which(is.na(fSym[,3]))
	 if(length(na)>0){
		fSym<-fSym[-na,]
	 }
	na<-which(is.na(fSym[,4]))
	 if(length(na)>0){
		fSym[na,4]<-"DD"
	 }
	 fSym<-data.frame(GID=as.integer(fSym[,1]),ENTREZID=as.integer(fSym[,1]),SYMBOL=as.character(fSym[,3]),GENENAME=as.character(fSym[,4]),stringsAsFactors =F)
	 
	 
	 fChr <- cbind(fSym[,1], rep(1,nrow(fSym)))
	 fChr<-data.frame(GID=as.integer(fChr[,1]),CHROMOSOME=fChr[,2])
	 colnames(fChr) <- c("GID","CHROMOSOME")



fGO<-c()
for(i in 1:nrow(fSym)){
	pos<-na.omit(match(fSym[i,"locus_tag"],anno[,13]))
	if(length(pos)>0){
		go<-anno[pos,2]
		go<-strsplit(go, split="; ")[[1]]
		tmp<-cbind(rep(i,length(go)),go)
		fGO<-rbind(fGO,tmp)
	}
}

if(is.null(fGO)){
	
	fSym <- cd[,c(4,2)]
	fSym<-cbind(seq(1:nrow(cd)),seq(1:nrow(cd)),toupper(as.character(cd[,"old_locus_tag"])),as.character(cd[,"gene"]))
	colnames(fSym) <- c("GID","ENTREZID","SYMBOL","GENENAME")
	
	dup<-which(duplicated(unlist(fSym[,3])))
	if(length(dup)>0){
		fSym<-fSym[-dup,]
	}
	na<-which(is.na(unlist(fSym[,1])))
	if(length(na)>0){
		fSym<-fSym[-na,]
	}
	na<-which(is.na(fSym[,3]))
	 if(length(na)>0){
		fSym<-fSym[-na,]
	 }
	na<-which(is.na(fSym[,4]))
	 if(length(na)>0){
		fSym[na,4]<-"DD"
	 }
	 fSym<-data.frame(GID=as.integer(fSym[,1]),ENTREZID=as.integer(fSym[,1]),SYMBOL=as.character(fSym[,3]),GENENAME=as.character(fSym[,4]),stringsAsFactors =F)
	 
	 
	 fChr <- cbind(fSym[,1], rep(1,nrow(fSym)))
	 fChr<-data.frame(GID=as.integer(fChr[,1]),CHROMOSOME=fChr[,2])
	 colnames(fChr) <- c("GID","CHROMOSOME")
	
	fGO<-c()
	for(i in 1:nrow(fSym)){
		pos<-na.omit(match(unlist(fSym[i,3]),toupper(anno[,13])))
		if(length(pos)>0){
			go<-anno[pos,2]
			go<-strsplit(go, split="; ")[[1]]
			tmp<-cbind(rep(i,length(go)),go)
			fGO<-rbind(fGO,tmp)
		}
	}

}
fGO<-data.frame(GID=as.integer(fGO[,1]),GO=fGO[,2], EVIDENCE=rep("NO",nrow(fGO)))

out<-list(fGO,fSym, fChr)
#save(out, file="out.Rdata")


## Then call the function
makeOrgPackage(gene_info=fSym, chromosome=fChr, go=fGO,
               version="0.1",
               maintainer="Some One <so@someplace.org>",
               author="Some One <so@someplace.org>",
               outputDir = ".",
               tax_id="blaba",
               genus="",
               species=gsub("_","",ooi),
               goTable="go")
install.packages(paste0("./org.",gsub("_","",ooi),".eg.db"), repos=NULL)
require(paste0("org.",gsub("_","",ooi),".eg.db"),character.only=TRUE)


}

term_list<-list()
go<-c("BP","MF","CC")
out_en<-list()
#node<-setdiff(tree[[1]][,1],tree[[1]][,2])  # root node
for(ii in 1:3){


ecoGO <- godata(paste0("org.",gsub("_","",ooi),".eg.db"), ont=go[ii])


IC<-c(ecoGO@IC)

	best<-c()
	pval<-c()
	term2<-c()
	# inner<-unique(tree[[1]][,1])
	# sel<-which(inner==node)
	# leaf<-leafs(inner[sel], tree)
	# leaf<-tree[[3]][leaf]
	# selection<-leaf
	# selction<-ord
	term<-c()
	
	#for(i in 1:min(length(selection),length(enrich_list))){
		tmp3<-enrich_list[[1]][[1]]
		#tmp<-enrich_list[[i]][[1]][[1]]
		#for(j in 1:3){
		j<-ii
		tmp2<-tmp3[[j]]
		#th<-which(tmp[,"classicFS"]<=thres )
		th2<-which(as.numeric(tmp2[,"classicFS"])<=thres & tmp2[,"Annotated"] > min_num & tmp2[,"Significant"] > min_num2 &  tmp2[,"Annotated"] < max_num)
		term2<-c(term2,tmp2[th2,1])
		pval<-c(pval,tmp2[th2,"classicFS"])
		tmp<-tmp2[th2,2]
		#tmp<-tmp[1:10]
		best<-c(best,tmp2[1:1,2])
		term<-c(term,tmp2[th2,2])
		#}
	#}
	k<-which(pval=="< 1e-30")
	if(length(k)>0){
		pval[k]<-"1e-30"
	}
	pval<-as.numeric(pval)
	term<-term[order(pval)]
	term2<-term2[order(pval)]
	pval<-sort(pval)
	count<-table(term)
	term_dup<-which(duplicated(term2))
	numbers<-table(term2)
	if(length(term_dup>0)){
		term<-term[-term_dup]
		term2<-term2[-term_dup]
		pval<-pval[-term_dup]
	}

	best<-unique(best)





term2<-unique(term2)
term_list[[ii]]<-term2
ma<-match(term2,names(numbers))

m2<-match(term2,names(IC))
IC_temp=IC[m2]

na<-which(is.na(IC_temp))

names(IC_temp)[na]<-"GO:1905039"
IC_temp[na]<-0
sim <- mgoSim(term2, term2,
                  semData=ecoGO,
                  measure="Wang",#Jiang
                  combine=NULL)
res<-data.frame(ID=term2,term,pVAL=pval,count=(numbers)[ma], IC=IC_temp)
res[which(is.finite(res$IC)==F),"IC"]<-max(is.finite(res$IC))

full_anno<-c()
for(i in 1:length(enrich_list)){
	tmp3<-enrich_list[[i]][[1]][[ii]]
	if(is.null(tmp3)==F){
		full_anno<-rbind(full_anno,as.matrix(tmp3))
	}
}


gene_cont<-vector("list", length(term2))
names(gene_cont)<-term2
for(i in 1:length(term2)){
	t1<-full_anno[grep(term2[i],full_anno[,1]),"genes"]
	t1<-unlist(strsplit(t1,","))
	gene_cont[[i]]<-t1
}

term_over<-matrix(,length(term2),length(term2))
for(i in 1:(length(term2))){
#print(i)
	for(j in (i):length(term2)){
		t1<-gene_cont[[i]]
		t2<-gene_cont[[j]]
		tmp<-length(intersect(t1,t2))/min(length(t1),length(t2))
		term_over[i,j]<-term_over[j,i]<-tmp
	}
}

number_sig_genes<-rep(NA, length(term2))
names(number_sig_genes)<-term2
for(i in term2){
	tmp<-which(full_anno[,1]==i & full_anno[,"classicFS"] <= thres)
	tmp<-full_anno[tmp,"genes"]
	tmp<-length(unlist(strsplit(tmp, ",")))
	number_sig_genes[i]<-tmp
}


com<-(((-log(as.numeric(res$pVAL))/max(-log(as.numeric(res$pVAL))))**1)*((res$count.Freq/max(res$count.Freq)))*res$IC/max(res$IC[is.finite(res$IC)]))**(1/3)


res$both<-com
res$number_sig_genes<-number_sig_genes

colnames(term_over)<-rownames(term_over)<-term2
ID<-term2	
go1 <- go2 <- overlap <- NULL

ov.df <- as.data.frame(term_over)
ov.df$go1 <- row.names(ov.df)
ov.df <- gather(ov.df, go2, overlap, -go1)
ov.df <- ov.df[!is.na(ov.df$overlap),]				  

	  
				  
by="both"				  


				  
go1 <- go2 <- similarity <- NULL

sim.df <- as.data.frame(sim)
sim.df$go1 <- row.names(sim.df)
sim.df <- gather(sim.df, go2, similarity, -go1)
sim.df <- sim.df[!is.na(sim.df$similarity),]

pos<-na.omit(match(paste(sim.df[,1],sim.df[,2]),paste(ov.df[,1],ov.df[,2])))

sim.df$overlap<-ov.df$overlap[pos]

## feature 'by' is attached to 'go1'
sim.df <- merge(sim.df, res[, c("ID", by)], by.x="go1", by.y="ID")



sim.df$go2 <- as.character(sim.df$go2)


out_en[[ii]]<-list(list(res,sim.df),best)

}

for(i in 1:3){
	tmp<-out_en[[i]][[1]][[1]]
	write.table(tmp, file=paste(go[i],"_enrichment_initial_set.txt",sep=""),sep="\t", quote=F, row.names=FALSE)
}

out_res<-list()

for(j in 1:3){

sim.df<-out_en[[j]][[1]][[2]]
res<-out_en[[j]][[1]][[1]]
best<-out_en[[j]][[2]]
ID<-res$ID
select_fun<-max	


sim.df<-sim.df[order(sim.df$both, decreasing=T),]

GO_to_remove <- character()


    for (i in seq_along(ID)) {
        ii <- which(sim.df$go2 == ID[i] & sim.df$similarity > cutoff & sim.df$overlap > thres_ov )

        if (length(ii) < 2)
            next

        sim_subset <- sim.df[ii,]

        jj <- which(sim_subset[, by] == select_fun(sim_subset[, by]))

        ## sim.df <- sim.df[-ii[-jj]]
        GO_to_remove <- c(GO_to_remove, sim_subset$go1[-jj]) %>% unique
    }				  
				  
out2<-res[!res$ID %in% GO_to_remove, ]		
best<-setdiff(best,out2[,2])
if(length(best)>0){
	best<-na.omit(match(best,res[,2]))
	if(length(best)>0){
		out2<-rbind(out2,res[best,])
	}
}
	out2<-out2[order(out2[,"both"], decreasing=T),]	  
	out_res[[go[j]]]<-out2	  
}

names(out_res)<-go


enri<-c()
# 	for(i in 1:length(selection)){
		tmp<-enrich_list[[1]][[1]][[1]]
		tmp2<-enrich_list[[1]][[1]][[2]]
		tmp3<-enrich_list[[1]][[1]][[3]]
		enri<-rbind(enri,tmp,tmp2,tmp3)


terml<-c()
pvall<-c()
fc<-c()
out<-c()
term_list<-list()
type<-c()
for(i in 1:3){
	temp<-out_res[[i]]
	
	if(nrow(temp)>0){
	type<-c(type,rep(i,nrow(temp)))
		for(j in 1:nrow(temp)){
			terml<-c(terml,paste0(temp[j,"term"],"~",temp[j,"ID"]))
			pvall<-c(pvall, temp[j,"pVAL"])
			#pos<-match(temp[j,"ID"], enri[,"GO.ID"])
			term<-temp[j,"ID"]	

			anno<-cbind(anno,an)
			tmp<-grep(term,enri[,1])
			tmp_fc<-enri[tmp,"Significant"]/enri[tmp,"Expected"]
			fc<-c(fc,tmp_fc)
			if(length(tmp)>0){
				genes<-strsplit(enri[tmp,"genes"],split=",")
				le<-unlist(lapply(genes,length))
				genes<-unlist(genes)
				pos_cop<-unlist(lapply(genes, grep,cop[,1], ignore.case=T))
				pval_cop<-copra_results[[1]][pos_cop,"p-value"]
				pval<-rep(enri[tmp,"classicFS"],le)
				ranks<-unlist(strsplit(enri[tmp,"ranks"],split=","))
				pos<-match(toupper(genes),toupper(anno[,ncol(anno)]))
				gene<-as.character(anno[pos,9])
				des<-as.character(anno[pos,12])
				org<-as.character(anno[pos,5])
				tmp1<-cbind(genes,gene,pval,org,ranks,des)
				tmp1<-cbind(rep(948746,le), pval_cop, genes, paste0(genes," - " , gene),paste0(genes," - " , gene),rep(NA,le),rep(NA,le),rep(NA,le),rep(NA,le),rep(NA,le),rep(NA,le))
				out<-rbind(out,tmp1)
				term_list[[length(term_list)+1]]<-genes
				names(term_list)[length(term_list)]<-paste0(temp[j,"term"],"~",temp[j,"ID"])
			}
		}
	}
}
na<-c("Entrez GeneID","CopraRNA_pValue","Locustag","GeneName","","IntaRNA energy","IntaRNA p-Value","start mRNA","end mRNA","start ncRNA","end ncRNA")


out<-unique(out)
for(i in 1:length(term_list)){
	tmp<-rep(0,nrow(out))
	tmp[match(term_list[[i]], out[,3])]<-1
	out<-cbind(out,tmp)
	na<-c(na, paste0(round(fc[i],digits=2)," ", terml[i]))
}

score<--log10(as.numeric(pvall))
score<-c(score,rep("",11))

clas<-c(nrow(out), length(terml), 1:length(terml), rep("",9))

out<-rbind(score,clas,na,out)

write.table(out, file="enrichment.txt", sep="\t",col.names = FALSE,row.names = FALSE, quote=FALSE)

