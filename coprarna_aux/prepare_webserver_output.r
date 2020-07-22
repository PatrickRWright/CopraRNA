# script by Jens Georg

# dependencies:

#call:
#R --slave -f /home/jens/CopraRNA-git/coprarna_aux/prepare_webserver_output.r 

suppressPackageStartupMessages(require(seqinr))

# get absolute path
initial.options <- commandArgs(trailingOnly = FALSE)
path<-initial.options [4]
path<-sub("copraRNA2_find_conserved_sites.r","",path)
#print(path)

# preset path to required files, path can also be specified as argument
copref_path<-paste(path,"CopraRNA_available_organisms.txt",sep="")

# read the organism of interest (ooi) from the ncRNA fasta file. The sRNA of the ooi is considered to be the first sequence.
ooi<-gsub("ncRNA_","",names(read.fasta("ncrna.fa"))[1])

# number of top predictions which should be investigated
co<-readLines("CopraRNA_option_file.txt") 
top<-as.numeric(gsub("top count:","",co[grep("top count:", co)]))

# path to CopraRNA result file (in order to investigate the top predictions)
copra_result<-"CopraRNA_result_all.csv"

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


copra<-read.csv(copra_result, sep=",")
top_tags<-gsub("\\(.*","",copra[,ooi])[1:top]
gene<-gsub(".*\\(","",copra[,ooi])[1:top]
gene<-gsub("\\|.*","",gene)

d<-dir()
int_ooi<-grep(paste(ooi, ".*.fa.intarna.csv",sep=""),d)
int_ooi<-read.csv(d[int_ooi],sep=";")
top_pos<-match(toupper(top_tags),toupper(int_ooi[,1]))
web_out<-int_ooi[top_pos,]
web_out<-cbind(web_out,gene)
web_out<-cbind(web_out,copra[1:top,"fdr"])
web_out<-cbind(web_out,copra[1:top,"p.value"])
web_out<-cbind(web_out,copra[1:top,"Annotation"])
colnames(web_out)[(ncol(web_out)-2):ncol(web_out)]<-c("fdr","p.value","Annotation")

write.table(web_out, file="coprarna_websrv_table.csv", quote=F, row.names=F, sep=",")