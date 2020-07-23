#!/usr/bin/env Rscript

#call
# R --slave -f ./copraRNA2_conservation_heatmaps.r --args top=50

inputfile="CopraRNA_result_all.csv"


suppressPackageStartupMessages(require(ComplexHeatmap))
suppressPackageStartupMessages(require(RColorBrewer))
suppressPackageStartupMessages(require(phangorn))
suppressPackageStartupMessages(require(seqinr))
suppressPackageStartupMessages(require(parallel))


# get absolute path
initial.options <- commandArgs(trailingOnly = FALSE)
path<-initial.options [4]
path<-sub("copraRNA2_conservation_heatmaps.r","",path)

# preset path to required files, path can also be specified as argument
coprarna_reference_file<-paste(path,"CopraRNA_available_organisms.txt",sep="")

# register cores for parallel processing
co<-readLines("CopraRNA_option_file.txt") 
max_cores<-as.numeric(gsub("core count:","",co[grep("core count:", co)]))

# number of top predictions which should be investigated
top<-as.numeric(gsub("top count:","",co[grep("top count:", co)]))

# if True the heatmaps are not done on the top predictions but on the locus tags given in the genelist.txt file
select=FALSE

# case sensitve list of locus tags (requires that the calculations from the "CopraRNA2_find_conserved_sites.r" are done on the selected genes
genelist="genelist.txt"

# Defines if the tree for the heatmap is calcualted based on the ribosomal RNAs or on the sRNA sequences
clustering="ribosomal" # alternative: "sRNA"

# empiric threshold for not being a target in the single organisms target prediction
# if the IntaRNA p-Vale is smaller than "int_p_thres" a black box is drawn in the heatmap
int_p_thres=0.35

# read the organism of interest (ooi) from the ncRNA fasta file. The sRNA of the ooi is considered to be the first sequence.
ooi<-gsub("ncRNA_","",names(read.fasta("ncrna.fa"))[1])
oois<-ooi


# pre-calculated objects from "CopraRNA2_find_conserved_sites.r"
load("int_sites.Rdata")
load("peak_list.Rdata")

# table of protein homolog clusters
clus_tab<-read.csv("cluster.tab", sep="\t", row.names=1)

# reference for row positions
evo_analysis<-read.csv("CopraRNA2_prep_anno_addhomologs_padj_amountsamp.csv",sep=",", header=T) 
evo_analysis<-as.matrix(evo_analysis)

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

int_p_thres<-as.numeric(int_p_thres)
top<-as.numeric(top)
num<-top

# function from ape package to transform a phylo to a hclust object
as.hclust.phylo2 <- function(x, ...)
{
    if (!is.ultrametric(x)) stop("the tree is not ultrametric")
    if (!is.binary.phylo(x)) stop("the tree is not binary")
    if (!is.rooted(x)) stop("the tree is not rooted")
    n <- length(x$tip.label)
    x$node.label <- NULL # by Jinlong Zhang (2010-12-15)
    bt <- branching.times(x)
    N <- n - 1L

    x <- reorder(x, "postorder")
    m <- matrix(x$edge[, 2], N, 2, byrow = TRUE)
    anc <- x$edge[c(TRUE, FALSE), 1]
    bt <- bt[as.character(anc)] # 1st, reorder
    ## 2nd, sort keeping the root branching time in last (in case of
    ## rounding error if there zero-lengthed branches nead the root)
    bt <- c(sort(bt[-N]), bt[N])
    o <- match(names(bt), anc)
    m <- m[o, ]

    ## first renumber the tips:
    TIPS <- m <= n
    m[TIPS] <- -m[TIPS]

    ## then renumber the nodes:
    oldnodes <- as.numeric(names(bt))[-N]
    m[match(oldnodes, m)] <- 1:(N - 1)

    names(bt) <- NULL
    obj <- list(merge = m, height = 2*bt, order = 1:n, labels = x$tip.label,
                call = match.call(), method = "unknown")
    class(obj) <- "hclust"
    obj
}

# wrapper to call mafft
mafft<-function(filename="ncrna.fa", outname="ncrna_aligned.fa", mode="accurate"){
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
	#fas<-read.fasta(outname)
	#fas
} 

# select the homologs set for drawing the heatmap
find_sRNA<-function(x, y){
	e<-grep("Annotation", colnames(y))-1
	temp<-evo_analysis[,3:e]
	genes<-apply(temp, 1, paste, collapse="")
	x<-as.matrix(x)
	out<-unlist(apply((x), 1, find1, genes))
	out
}
find1<-function(x, genes){
	temp<-grep(x, genes, ignore.case=F)
	if(length(temp)>0){
		temp<-paste(temp,x,sep=",")
	}
	temp
}

if(select==T){
	sel=genelist # case sensitve list of locus tags
	genelist<-read.csv(sel, sep="\t", header=F)[,1]
	selection2<-unique(find_sRNA(genelist, evo_analysis))
	selection<-as.numeric(gsub(",.*","",selection2))
	genelist<-gsub(".*,","",selection2)
	
}
if(select==F){
	res<-read.csv(inputfile, sep=",")
	num<-min(num,nrow(res))
	selection<-as.character(res[1:num,4])
	selection<-as.numeric(as.character(res[1:num,"initial_sorting"]))
	genelist<-gsub("\\(.*","",res[1:num,4])

}

ee<-grep("Annotation", colnames(evo_analysis))-1
homologs<-evo_analysis[selection,3:ee]
org_names<-colnames(evo_analysis)[3:ee]
all_orgs<-colnames(evo_analysis)[3:ee]


# transform list of sites to matrix
to_table1<-function(x, pos=1, norgs){
	x<-gsub("\\|.*","",x[pos,])
	if(length(x)==0){
		x<-rep(NA,norgs)
	}
	x
}

to_table2<-function(x, pos=1, norgs){
	y<-lapply(x,to_table1, pos=pos, norgs=norgs)
	out<-c()
	for(i in 1:length(y)){
		out<-rbind(out,y[[i]])
	}
	out
}

int_opt<-to_table2(int_sites[selection], pos=1,norgs=length(org_names))
int_sub<-to_table2(int_sites[selection], pos=2,norgs=length(org_names))
colnames(int_opt)<-colnames(evo_analysis)[3:ee]
colnames(int_sub)<-colnames(evo_analysis)[3:ee]
peaks<-peak_list[selection]



# calculate phylogenetic tree for heatmap clustering
if(clustering=="sRNA"){
	temp_align<-tempfile()
	temp_align2<-tempfile()
	mafft(filename="input_sRNA.fa", outname=temp_align)
	tempf<-read.fasta(temp_align)
	write.fasta(tempf, file.out=temp_align2, names=names(tempf), nbchar=100000)
	dat<-read.phyDat(temp_align2, format="fasta", type="DNA")
	dm <- dist.ml(dat, model="F81")
	treeNJ <- NJ(dm)
	fitJC = pml(treeNJ, data=dat)
	fit2<-midpoint(fitJC$tree)
	fit2<-chronos(fit2)
	clus<-as.hclust.phylo2((fit2))
	ord<-clus$label
	ord2<-match(ord, gsub("\\..*","",colnames(int_opt)))
}
if(clustering=="ribosomal"){
	
	load("16S_tree.Rdata")
	fit2<-midpoint(fit2)
	if(length(fit2$tip.label)==(length(unique(fit2[[1]][,1]))+1)){
		fit2<-chronos(fit2)
	} else {
		dat<-read.phyDat("16S_aligned.fa", format="fasta", type="DNA")
		dm <- dist.ml(dat, model="F81")
		treeNJ <- NJ(dm)
		fitJC = pml(treeNJ, data=dat)
		fit2<-fitJC$tree
		fit2<-midpoint(fit2)
		fit2<-chronos(fit2)
	}
	clus<-as.hclust.phylo2((fit2))
	ord<-clus$label
	ord2<-match(ord, gsub("\\..*","",colnames(int_opt)))
}

# calculate p_values based on IntaRNA energies
intarna_link<-function(all_orgs){
		tmp<-grep(paste(all_orgs,"_.*_opt.intarna.csv", sep=""), dir())
		tmp_sub<-grep(paste(all_orgs,"_.*_subopt.intarna.csv", sep=""), dir())
		tmp<-as.matrix(read.csv(dir()[tmp],sep=";"))
		tmp_sub<-as.matrix(read.csv(dir()[tmp_sub],sep=";"))
		tmp_out<-matrix(,nrow(tmp),2)
		tmp_out[]<-0
		row.names(tmp_out)<-tolower(tmp[,1])
		tmp_out[,1]<-as.numeric(tmp[,"E"])
		tmp_match<-match(tmp[,1],tmp_sub[,1])
		tmp_ee<-c()
			for(j in 1:length(tmp_match)){
				if(is.na(tmp_match[j])==F){
					tmp_out[j,2]<-as.numeric(tmp_sub[tmp_match[j],"E"])
				}
			}
			
	tmp_out
}


e<-grep("Annotation", colnames(evo_analysis))
extreme<-mclapply(all_orgs, intarna_link,mc.cores=max_cores)
names(extreme)<-all_orgs



pvalue<-function(int_opt,int_sub, evo_analysis, selection, thres=0.35){
	int_optt<-int_opt
	int_subt<-int_sub
	for(i in 1:nrow(int_opt)){
		gene_names<-gsub("\\(.*","",evo_analysis[selection[i],3:ee])
		for(j in 1:length(gene_names)){
			if(gene_names[j]!=""){
				temp<-extreme_ind[[names(gene_names)[j]]][gene_names[j],]
				if(temp[1]<=thres){
					int_optt[i,j]<-paste(int_opt[i,j],"_sig1",sep="")
				}
				if(temp[2]<=thres){
					int_subt[i,j]<-paste(int_sub[i,j],"_sig1",sep="")
				}
			}
		}
	}
	out<-list(int_optt,int_subt)
	out
}

pval<-function(energy_vector, color=1){
	empty<-which(energy_vector>=0)
	if(length(empty)>0){
		energy_vector<-energy_vector[-empty]
	}
	energy_vector<-energy_vector*(-1)
	gevenergies <- gev(energy_vector) # fit
	xi <- gevenergies$par.ests[1] # read chi from fitting
	sigma <- gevenergies$par.ests[2] # read sigma
	mu <- gevenergies$par.ests[3] # read mu
	out<-list(xi=xi,sigma=sigma,mu=mu)
	out
}
 
gev <- function (data, block = NA, ...)
{
    n.all <- NA
    if (!is.na(block)) {
        n.all <- length(data)
        if (is.character(block)) {
            times <- as.POSIXlt(attributes(data)$times)
            if (block %in% c("semester", "quarter")) {
                sem <- quart <- times$mon
                sem[sem %in% 0:5] <- quart[quart %in% 0:2] <- 0
                sem[sem %in% 6:11] <- quart[quart %in% 3:5] <- 1
                quart[quart %in% 6:8] <- 2
                quart[quart %in% 9:11] <- 3
            }
            grouping <- switch(block, semester = paste(times$year,
                sem), quarter = paste(times$year, quart), month = paste(times$year,
                times$mon), year = times$year, stop("unknown time period"))
            data <- tapply(data, grouping, max)
        }
        else {
            data <- as.numeric(data)
            nblocks <- (length(data)%/%block) + 1
            grouping <- rep(1:nblocks, rep(block, nblocks))[1:length(data)]
            data <- tapply(data, grouping, max)
        }
    }
    data <- as.numeric(data)
    n <- length(data)
    sigma0 <- sqrt(6 * var(data))/pi
    mu0 <- mean(data) - 0.57722 * sigma0
    xi0 <- 0.1
    theta <- c(xi0, sigma0, mu0)
    negloglik <- function(theta, tmp) {
        y <- 1 + (theta[1] * (tmp - theta[3]))/theta[2]
        if ((theta[2] < 0) || (min(y) < 0))
            out <- 1e+06
        else {
            term1 <- length(tmp) * logb(theta[2])
            term2 <- sum((1 + 1/theta[1]) * logb(y))
            term3 <- sum(y^(-1/theta[1]))
            out <- term1 + term2 + term3
        }
        out
    }
    fit <- optim(theta, negloglik, hessian = TRUE, ..., tmp = data)
    if (fit$convergence)
        warning("optimization may not have succeeded")
    par.ests <- fit$par
    varcov <- solve(fit$hessian)
    par.ses <- sqrt(diag(varcov))
    out <- list(n.all = n.all, n = n, data = data, block = block,
        par.ests = par.ests, par.ses = par.ses, varcov = varcov,
        converged = fit$convergence, nllh.final = fit$value)
    names(out$par.ests) <- c("xi", "sigma", "mu")
    names(out$par.ses) <- c("xi", "sigma", "mu")
    class(out) <- "gev"
    out
}

pgev <- function (q, xi = 1, mu = 0, sigma = 1)
{
    exp(-(1 + (xi * (q - mu))/sigma)^(-1/xi))
}

p_vals<-list()
	
for(i in 1:length(extreme)){
	opti<-pval(c(extreme[[i]][,1]))
	if(min(extreme[[i]][,2]<0)){
		subopti<-pval(c(extreme[[i]][,2]))
	} else {
		subopti<-rep(NA,length(extreme[[i]][,2]))
	}
	li<-list(opt=opti,subopti=subopti)
	p_vals[[i]]<-li
}
names(p_vals)<-names(extreme)

extreme_ind<-extreme

for(i in 1:length(extreme)){
	extreme_ind[[i]][,1]<-1-pgev(extreme[[i]][,1]*(-1),p_vals[[i]][[1]]$xi,p_vals[[i]][[1]]$mu, p_vals[[i]][[1]]$sigma )
	extreme_ind[[i]][,2]<-1-pgev(extreme[[i]][,2]*(-1),p_vals[[i]][[1]]$xi,p_vals[[i]][[1]]$mu, p_vals[[i]][[1]]$sigma )
}

# read availbale organism file for annotation purposes
copref<-read.delim(coprarna_reference_file, sep="\t", header=T,comment.char = "#")	
nam<-c()
for(i in 1:ncol(int_opt)){
	tnam<-grep(gsub("\\..*","",colnames(int_opt)[i]),copref[,1])
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
nam2<-paste(nam2, colnames(int_opt), sep="_")
selection<-na.omit(selection)


gene_anno<-c()
gene_anno2<-c()
evo_analysis<-as.matrix(evo_analysis)
for(i in 1:length(selection)){
	co<-grep(genelist[i],evo_analysis[selection[i],])[1]
	temp<-evo_analysis[selection[i],co]
	if(length(temp)>0){
		locus_tag<-gsub("\\(.*","",as.character(temp))
		genename<-gsub(".*\\(","",as.character(temp))
		genename<-gsub("\\|.*", "", genename)
		genename<-gsub("\\/", "", genename)
		genename2<-genename
		if(length(genename)==1 && genename[1]=="NA" ){
			genename<-evo_analysis[selection[i], "Annotation"]
			genename<-substr(genename,1,min(nchar(genename),20))
		}
		na<-paste(genename, locus_tag,sep="_")
		na2<-paste(genename2, locus_tag,sep="_")
	} else{
		na<-selection[i]
	}
	gene_anno<-c(gene_anno, na)
	gene_anno2<-c(gene_anno2, na2)
}

l<-2
len<-length(which(duplicated(gene_anno)))
while(len>0){
	gene_anno[which((duplicated(gene_anno)))]<-paste(gene_anno[which((duplicated(gene_anno)))],l,sep="_")
	l<-l+1
	len<-length(which(duplicated(gene_anno)))
}
row.names(int_opt)<-gene_anno



temp<-pvalue(int_opt,int_sub,evo_analysis, selection, thres=int_p_thres)
int_opt2<-temp[[1]]
int_sub2<-temp[[2]]
	

int_opt<-int_opt[,ord2]
int_sub<-int_sub[,ord2]
int_opt2<-int_opt2[,ord2]
int_sub2<-int_sub2[,ord2]
colnames(int_opt)<-nam2[ord2]

# tranform the conserved peaks into integers and assign colors for heatmap drawing 
to_number<-function(int_opt, int_sub,homologs, ooi, peaks){
	colo2<-c("#caff70","#7aa9df","#cf97bb" ,"#fede7e","#fc9187","#d3648f","#4c9998","#988377","#576e81","#9fe1e5")
	colo<-c("#799943","#5680b0","#9c6488","#e4c054","#e26a5f","#b43768","#006362","#614736","#294258","#77b1b5")
	ooi_pos<-match(ooi,colnames(int_sub))
	int<-matrix(,nrow(int_opt),ncol(int_opt))
	rownames(int)<-rownames(int_opt)
	colnames(int)<-colnames(int_opt)
	int[]<--1
	int_o<-int
	int_s<-int
	color_vect<-c(-1,0)
	names(color_vect)<-c("#FFFFFF", "#b0b0b0")
	for(i in 1:nrow(int_opt)){
		peaks1<-peaks[[i]]
		exist<-which(homologs[i,]!="")
		if(length(exist)>0){
			int_o[i,exist]<-0
		}
		temp1<-na.omit(unique(c(int_opt[i,])))
		temp2<-na.omit(unique(c(int_sub[i,])))
		temp<-unique(sort(c(temp1,temp2)))
		
		ooi_peak<-int_opt[i,ooi_pos]
		
		if(is.na(ooi_peak)==F){
			ooi_p2<-grep(ooi_peak,peaks1)
			peaks1<-c(peaks1[ooi_p2],peaks1[-ooi_p2])
			
		}
		le<-length(temp)
		if(length(peaks1)>0){
			colo3<-c()
			for(j in 1:length(peaks1)){
				fc <- colorRampPalette(c(colo[j], colo2[j]))
				tmp<-fc(length(peaks1[[j]]))
				names(tmp)<-peaks1[[j]]
				colo3<-c(colo3,tmp)
			}
			no_ex<-which(is.element(colo3, names(color_vect))==F)
			if(length(no_ex)>0){
				tmp<-(max(color_vect)+1):(max(color_vect)+length(no_ex))
				names(tmp)<-colo3[no_ex]
				color_vect<-c(color_vect,tmp)
			}
			for(j in 1:le){
				tmp<-match(temp[j],names(colo3))
				tmp<-match(colo3[tmp],names(color_vect))
				tmp<-color_vect[tmp]
				pos<-which(int_opt[i,]==temp[j])
				
				int_o[i,pos]<-tmp
				pos<-which(int_sub[i,]==temp[j])
				int_s[i,pos]<-tmp
			}
		}
	}
	out<-list(int_o, int_s, color_vect)
	out
}

int_opt2<-gsub(".*_","",int_opt2)
int_sub2<-gsub(".*_","",int_sub2)
int_opt3<-int_opt
int_sub3<-int_sub
int_opt3[]<-"-"
int_sub3[]<-"-"
int_opt_col<-int_opt
int_opt_col[]<-"orangered2"
int_sub_col<-int_sub
int_sub_col[]<-"orangered2"
s1<-which(int_opt2=="sig1",arr.ind=T)
int_opt3[s1]<-"++"
int_opt_col[s1]<-"white"
s2<-which(int_opt2=="sig2",arr.ind=T)
int_opt3[s2]<-"+"
int_opt_col[s2]<-"white"
s1<-which(int_sub2=="sig1",arr.ind=T)
int_sub3[s1]<-"++"
int_sub_col[s1]<-"white"
s2<-which(int_sub2=="sig2",arr.ind=T)
int_sub3[s2]<-"+"
int_sub_col[s2]<-"white"

int<-to_number(int_opt,int_sub,homologs,ooi, peaks)
int_opt<-int[[1]]
int_sub<-int[[2]]
min_numb<-min(as.vector(int_opt), na.rm=T)
min_numb<-which(int[[3]]==min_numb)
my_palette<-names(int[[3]])[min_numb:length(int[[3]])]
my_palette2<-my_palette

mp<-my_palette2

# heat map legend
leg<-c("no homolog","not conserved")			
legen<-c("#FFFFFF", "#b0b0b0")
lgd = Legend(labels = leg, title = 'site\nconservation',title_position = "topcenter", legend_gp = gpar(fill = legen),border = "black",title_gp = gpar(font=2))
lgd_list = list(lgd)

# draw combined heatmap
pdf("conservation_heatmap.pdf", width = 2.3+ncol(int_opt)*0.2, height = 3.4+nrow(int_opt)*0.2,useDingbats=F)
col_col<-rep("black", ncol(int_opt))
col_col[clus[[3]][match(oois, clus[[4]])]]<-"orangered"
c_lty<-rep(1, ncol(int_opt))
c_lty[clus[[3]][match(oois, clus[[4]])]]<-2
a<-Heatmap(	int_opt,
			col=my_palette,
			cluster_rows = F, 
			cluster_columns =clus, 
			show_heatmap_legend = F,
			column_names_gp = gpar( rot = 30, cex=0.7, col = col_col, font=c_lty),
			row_names_gp = gpar( rot = 0, cex=0.75,font=3),
			cell_fun=function(j, i, x, y, width=width, height=width, fill){
				s = min(unit.c(convertWidth(width, "cm"), convertHeight(height, "cm")))
				
				x1<-as.numeric(x)
				y1<-as.numeric(y)
				w1<-as.numeric(width)
				h1<-as.numeric(height)
				
				rect_gp = gpar(col = "lightgrey", lty = 1, lwd = 0.15)
				if(int_sub[i,j]>-1){
					grid.polygon(x = c(x1-0.45*w1,x1+0.45*w1,x1-0.45*w1), y = c(y1-0.45*h1,y1+0.45*h1,y1+0.45*h1),  
						gp = gpar(fill=my_palette2[int_sub[i,j]+2], col=my_palette2[int_sub[i,j]+2]))
				}
				if((int_sub[i,j]>=1 &  int_sub3[i,j]=="++") ){
					grid.polygon(x = c(x1-0.42*w1,x1+0.42*w1,x1+0.42*w1,x1-0.42*w1), y = c(y1-0.42*h1,y1-0.42*h1,y1+0.42*h1,y1+0.42*h1),  
						gp = gpar( col=1, lwd=0.9))
				}
				if( int_opt[i,j]>=0 & int_opt3[i,j]=="++"  ){
					grid.polygon(x = c(x1-0.42*w1,x1+0.42*w1,x1+0.42*w1,x1-0.42*w1), y = c(y1-0.42*h1,y1-0.42*h1,y1+0.42*h1,y1+0.42*h1),  
						gp = gpar( col=1, lwd=0.9))
				}
				},
			rect_gp=gpar(col=c("grey77"),lwd=c(1),lty=c(1))
			
			)
draw(a, annotation_legend_list = lgd_list)
dev.off()


# draw cluster specific heatmap
for(jj in 1:length(selection)){	
	int_opt2<-t(as.matrix(int_opt[jj,]))
	int_sub2<-t(as.matrix(int_sub[jj,]))
	int_opt4<-t(as.matrix(int_opt3[jj,]))
	int_sub4<-t(as.matrix(int_sub3[jj,]))
	col_col<-rep("black", ncol(int_opt2))
	col_col[clus[[3]][match(oois, clus[[4]])]]<-"orangered"
	c_lty<-rep(1, ncol(int_opt2))
	c_lty[clus[[3]][match(oois, clus[[4]])]]<-2
	se<-seq(-1,100)
	s<-match(min(int_opt2),se)
	e<-match(max(int_opt2),se)
	my_palette3<-my_palette[s:e]
	s<-match(min(int_sub2),se)
	e<-match(max(int_sub2),se)
	my_palette4<-my_palette2[s:e]
	na3<-paste("./evo_alignments2/",gene_anno2[jj],"/",gene_anno2[jj],"_conservation_heatmap.pdf",sep="")
	if(file.exists(paste("./evo_alignments2/",gene_anno2[jj],"/",sep=""))){
		pdf(na3, width = 2.3+ncol(int_opt2)*0.2, height = 3.4+nrow(int_opt2)*0.2,useDingbats=F)
		a<-Heatmap(	int_opt2,
					col=my_palette3,
					cluster_rows = F, 
					cluster_columns =clus, 
					column_title=rownames(int_opt)[jj],
					show_heatmap_legend = F,
					column_names_gp = gpar( rot = 30, cex=0.7, col = col_col, font=c_lty),
					row_names_gp = gpar( rot = 0, cex=0.75),
					cell_fun=function(j, i, x, y, width=width, height=width, fill){
						s = min(unit.c(convertWidth(width, "cm"), convertHeight(height, "cm")))
						x1<-as.numeric(x)
						y1<-as.numeric(y)
						w1<-as.numeric(width)
						h1<-as.numeric(height)
						rect_gp = gpar(col = "lightgrey", lty = 1, lwd = 0.15)
						if(int_sub2[i,j]>-1){
							grid.polygon(x = c(x1-0.45*w1,x1+0.45*w1,x1-0.45*w1), y = c(y1-0.45*h1,y1+0.45*h1,y1+0.45*h1),  
							gp = gpar(fill=my_palette4[int_sub2[i,j]+2], col=my_palette4[int_sub2[i,j]+2]))
						}
						if((int_sub2[i,j]>=1 &  int_sub4[i,j]=="++") ){
							grid.polygon(x = c(x1-0.42*w1,x1+0.42*w1,x1+0.42*w1,x1-0.42*w1), y = c(y1-0.42*h1,y1-0.42*h1,y1+0.42*h1,y1+0.42*h1),  
							gp = gpar( col=1, lwd=0.9))
						}
						if( int_opt2[i,j]>=0 & int_opt4[i,j]=="++"  ){
							grid.polygon(x = c(x1-0.42*w1,x1+0.42*w1,x1+0.42*w1,x1-0.42*w1), y = c(y1-0.42*h1,y1-0.42*h1,y1+0.42*h1,y1+0.42*h1),  
							gp = gpar( col=1, lwd=0.9))
						}
					},
					rect_gp=gpar(col=c("grey77"),lwd=c(1),lty=c(1))
					)
		draw(a, annotation_legend_list = lgd_list)
		dev.off()
	}
}
