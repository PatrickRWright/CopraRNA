args <- commandArgs(trailingOnly = TRUE) ## edit 2.0.5.1
inputFile <- args[1] ## edit 2.0.5.1

pvals <- read.table(inputFile, header=TRUE)
pvals <- unlist(pvals)
padj <- p.adjust(pvals, method="BH") ## edit 1.2.9 Benjamini Hochberg instead of qvalue
write.table(padj,file="padj.csv", row.names=F, col.names='fdr', quote=F) ## edit 1.2.9 changed col.names to fdr

