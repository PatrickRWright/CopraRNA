
# GLASSgo postprocessing script

# dependencies: 
suppressPackageStartupMessages(require(phangorn))
suppressPackageStartupMessages(require(seqinr))

#call:
# R --slave -f /home/jens/CopraRNA-git/coprarna_aux/coprarna_selection_fast.r --args mindis=0 maxorgs=15 maxdis=0.5


# parameters from function call:

mindis<--0.1	# required minimal distance of organisms
maxorgs<-15		# number of final organisms for prediction
maxdis<-0.5		# maximal distance between organisms in the prediction

# get absolute path
initial.options <- commandArgs(trailingOnly = FALSE)
path<-initial.options [4]
path<-sub("coprarna_selection_fast.r","",path)
print(path)

# preset path to required files, path can also be specified as argument
cop_path<-paste(path,"CopraRNA_available_organisms.txt",sep="")
print(cop_path)


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

mindis<-as.numeric(mindis)
maxorgs<-as.numeric(maxorgs)
maxdis<-as.numeric(maxdis)

outorg="internal"
exact_tree=TRUE


divide<-function(number,candlist,p4,max_step=5,ooi){
	ooi_p<-grep(ooi, colnames(p4))
	#p5<-names(sort(p4[,ooi_p))
	more<-0
	stepsize<-floor(length(candlist)/number)
	stepsize<-min(max_step,stepsize)
	if(number>1){
		sel<-seq(1,length(candlist)*2,by=stepsize)
		sel<-sel[1:number]
		sel<-na.omit(sel)
		more<-number-length(sel)
		if(length(candlist)>1){
			cands<-names(sort(p4[candlist,ooi_p]))
		}
		if(length(candlist)==1){
			cands<-rownames(p4)[candlist]
		}
	}
	if(number==1){
		sel<-1
		#cands<-rownames(p4)[candlist]
		if(length(candlist)>1){
			cands<-names(sort(p4[candlist,ooi_p]))
		}
		if(length(candlist)==1){
			cands<-rownames(p4)[candlist]
		}
	}
	
	sel<-cands[sel]
	if(more>0){
		r<-setdiff(cands,sel)
		sel<-c(sel,r[1:more])
		}	
	#print(c(i,sel))
	sel	
}

# call mafft for MSA
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

# function to exclude very similar organism based on a phylogentic tree to reduce complexity
exclude_similars<-function(dis, thres=0.01,ooi){
i<-1

o<-grep(ooi, rownames(dis))
nam<-colnames(dis)[o]
	temp<-which(dis[,o]<=thres)
	nam<-na.omit(match(nam, rownames(dis)))
	nam<-na.omit(match(nam, temp))
	if(length(nam)>0){
		temp<-temp[-nam]
	}
	
	if(length(temp)>0){
		dis<-dis[-temp,-temp]
		
	}
	



while(nrow(dis)>i){
	
	nam<-colnames(dis)[i]
	temp<-which(dis[,i]<=thres)
	nam<-na.omit(match(nam, rownames(dis)))
	nam<-na.omit(match(nam, temp))
	if(length(nam)>0){
		temp<-temp[-nam]
	}
	
	if(length(temp)>0){
		dis<-dis[-temp,-temp]
		
	}
	i<-i+1
}
dis
}





mafft(filename="16s_sequences.fa",outname="16s_aligned.fasta")

# create a ML tree based on the alignment
tempf<-read.fasta("16s_aligned.fasta")
write.fasta(tempf, file.out="16s_aligned.fasta", names=names(tempf), nbchar=100000)
dat<-read.phyDat("16s_aligned.fasta", format="fasta", type="DNA")
dm <- dist.ml(dat, model="F81")
dm2<-as.matrix(dm)
	

treeNJ <- NJ(dm)
if(exact_tree==T){
	fitStart = pml(treeNJ, dat, k=4)
	fitJC = optim.pml(fitStart, model="GTR", optGamma=T, rearrangement="stochastic",ratchet.par = list(iter = 5L, maxit = 20L, prop = 1/3),control = pml.control(epsilon = 1e-08, maxit = 10,
	trace = 1L))
}
if(exact_tree==F){
	fitJC = pml(treeNJ, data=dat)
}

fit2<-fitJC
fitJC<-fit2
lab<-fitJC$tree$tip.label
lab2<-lab

copref<-read.delim(cop_path, sep="\t", header=T,comment.char = "#")	
nam<-c()
for(i in 1:length(tempf)){
	tnam<-grep(gsub("\\..*","",names(tempf)[i]),copref[,1])
	nam<-c(nam,as.character(copref[tnam,2]))
	
	
}
nam2<-c()
for(i in 1:length(nam)){
	temp1<-substr(nam[i],1,3)
	temp2<-strsplit(nam[i],"_")[[1]]
	temp1<-paste(temp1,"_",temp2[2], sep="")
	if(length(temp2)>2){
		temp1<-paste(temp1, temp2[length(temp2)], sep="_")
	}
	nam2<-c(nam2,temp1)
}
nam2<-paste(nam2,names(tempf), sep="_")

orgs_selection<-vector("list",length(tempf))

pdf("coprarna_selection_tree.pdf")
for(jj in 1:length(tempf)){	
	ooi<-names(tempf)[jj]

p<-cophenetic(fitJC$tree)
if(outorg=="internal"){
	ooip2<-grep(ooi,colnames(p))
	outorg<-names(sort(p[,ooip2])[nrow(p)])
}
p3<-exclude_similars(p, mindis,ooi=ooi)

if(nrow(p)>=maxorgs){
	while(nrow(p3)<maxorgs){
		mindis<-mindis-mindis/100
		p3<-exclude_similars(p, mindis,ooi=ooi)
	}
}

ooip<-grep(ooi,colnames(p3))
sel<-names(sort(p3[,ooip])[which(sort(p3[,ooip])<=maxdis)][1:maxorgs])
tree<-fitJC$tree


sel<-unique(c(ooi,sel))
color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]

pos<-match( tree$tip.label, names(tempf))
tree$tip.label<-nam2[pos]
selected1<-match(unique(c(sel)),lab2)
selected3<-grep(ooi,lab2)[1]
colo<-rep("1",length(lab))

colo[selected1]<-"dodgerblue"
colo[selected3]<-"olivedrab3"
 

plot(tree, tip.color=colo, cex=0.5)
add.scale.bar()


orgs_selection[[jj]]<-na.omit(sel)
}
dev.off()

names(orgs_selection)<-names(tempf)

save(orgs_selection, file="orgs_selection.Rdata")


        