args <- commandArgs(trailingOnly = TRUE)
inputFile <- args[1]

pvals <- read.table(inputFile, header=TRUE)
pvals <- unlist(pvals)
padj <- p.adjust(pvals, method="BH")
write.table(padj,file="padj.csv", row.names=F, col.names='fdr', quote=F)

