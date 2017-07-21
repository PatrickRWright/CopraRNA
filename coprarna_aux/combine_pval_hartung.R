## edit 2.0.4.1 // changed combination to the hartung method without root weights
weights <- weights[order(weights$V1),]
##weights$V2 <- (weights$V2)^(1/2.5)
weights <- weights$V2

curr_clust_indices <- which(data$clusternumber==cl_num)
curr_clust_data <- data[curr_clust_indices,]
curr_clust_ltags <- curr_clust_data$d1 ## edit 2.0.5.1 changed to d1

curr_clust_data <- curr_clust_data[order(curr_clust_data$id2),] ## edit 2.0.5.1 chaned to id2
pvals <- curr_clust_data$p.value
# "A Note on Combning Dependent Tests of Significance" Joachim Hartung
combined_pvalue <- Hartung(pvals, lambda=weights)$p.value ## edit 2.0.4.2

final_line <- c(combined_pvalue,as.vector(curr_clust_ltags),"") # the "" is to add a semicolon to the end of the line
final_line <- toString(gsub(", ",";",toString(final_line)))

write.table(final_line, file="CopraRNA_result_no_pvsample.csv", append=T, quote=F, col.names=F, row.names=F)

