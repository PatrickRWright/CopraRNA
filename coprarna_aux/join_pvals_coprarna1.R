## edit 2.0.4 // changed script name to join_pvals_coprarna1.R

args <- commandArgs(trailingOnly = TRUE) ## edit 2.0.5.1
inputFile <- args[1] ## edit 2.0.5.1
prep_for_cop2 <- args[2] ## edit 2.0.6

####################### retrieve cov matrix from max clusters

options <- read.table("CopraRNA_option_file.txt", sep=":")
root <- as.numeric(as.character(options$V2[14]))

fasta <- read.table("ncrna.fa")
orgcnt <- length(grep("^>", fasta$V1))

therna <- "ncrna_"

# read weights
weightdata <- read.table("zscore.weight", sep=";", header=FALSE)
weight <- ((as.numeric(weightdata$V2))**(1/(root))) ## edit 2.0.5.1 // parametrized root function
weightdata <- cbind(weightdata,weight)
weightdata <- weightdata[order(weightdata$V1),]
rownames(weightdata)<-NULL
refids <- as.character(weightdata$V1)

# read clusters
data <- read.table(inputFile, sep=";", header=TRUE) # get data from hash clusters output // ## edit 2.0.5.1 // generic input file

# find full clusters
int <- data$clusternumber # get the clusternumbers as vect
maxnums <- c()
for (i in 1:int[length(int)]) { # fuer jede clusternumber
    if(length(int[which(int==i)]) == orgcnt) {
        maxnums <- c(maxnums, i)
    }
}

# extract the lines which participate in maxclusters
maxclustdata <- data[data$clusternumber %in% maxnums, ]

#get the organism names (refseqs)
names <- unique(maxclustdata$id2) ## edit 2.0.5.1 // changed to id2

#re-sort names vector to be the same as refids vector
namesresorted <- c()
for(i in 1:length(refids)) {
    check <- grep(refids[i], names, value=FALSE)
    namesresorted[i] <- as.character(names[check])
}

# do it for first as specific to initialize the dataframe
trueindices <- which(maxclustdata$id2==namesresorted[1]) ## edit 2.0.5.1 // changed to id2
truelines <- maxclustdata[trueindices, ]
thename <- sub(therna, "", namesresorted[1])
# transform the p-values to a normal distribution
newframe <- data.frame(qnorm(truelines$p.value,lower.tail=FALSE))
names(newframe)[1]<-paste(thename)

for (i in 2:orgcnt) {
    # extract line indices from dataframe for org i
    trueindices <- which(maxclustdata$id2==namesresorted[i]) ## edit 2.0.5.1 // changed to id2
    # put the actual lines in a vector
    truelines <- maxclustdata[trueindices, ]
    thename <- sub(therna, "", namesresorted[i])
    # transform the p-values to a normal distribution and add to dataframe
    newframe <- data.frame(newframe, qnorm(truelines$p.value, lower.tail=FALSE))
    names(newframe)[i]<-paste(thename)
}

# remove the inf lines
newframe <- newframe[is.finite(rowSums(newframe)),]

covmatrix <- cov(newframe)
means <- c()

means <- colMeans(newframe)

################################# end cov matrix retrieval


# initialize condNormal function

condNormal <- function(x.given, mu, sigma, given.ind, req.ind){
    # Returns conditional mean and variance of x[req.ind]
    # Given x[given.ind] = x.given
    # where X is multivariate Normal with
    # mean = mu and covariance = sigma

    B <- sigma[req.ind, req.ind]
    C <- sigma[req.ind, given.ind, drop=FALSE]
    D <- sigma[given.ind, given.ind]
    CDinv <- C %*% solve(D)
    cMu <- c(mu[req.ind] + CDinv %*% (x.given - mu[given.ind]))
    cVar <- B - CDinv %*% t(C)
    list(condMean=cMu, condVar=cVar)
}


# intervals seq(0,1,by=0.1), seq(prelim.rho-0.1,prelim.rho+0.1,by=0.01)               ### first call                           ### second call
compute_rho <- function(interval, data, weightdata, means, covmatrix) {  # left:0 right:1 data:tags.clusterd // left:prelim.rho-0.1 right:prelim.rho+0.1 data:tags.clusterd
    error <- 1000000000 # initialize
    optimal_rho <- interval[1]
    final.out_lines <- c()
    # for every rho
    for (rho in interval) {
        all_pvalues <- c()
        final.out_lines_temp <- c()
        # for every cluster
        for (i in unique(data$clusternumber)) { ## edit 2.0.0

            curr_cluster<-data[which(data$clusternumber==i),]
            curr_cluster$id2<-sub("ncrna_", "", curr_cluster$id2) ## edit 2.0.5.1 // changed to id2
            curr_cluster <- curr_cluster[order(curr_cluster$id2),] ## edit 2.0.5.1 // changed to id2
            rownames(curr_cluster)<-NULL
            p <- curr_cluster$p.value
           
            # fix for problem with p-values == 1 // qnrom(1) returns Inf
            p[which(p==1)] <- 0.9999999999999999 
 
           if (length(curr_cluster$clusternumber) == orgcnt ) {

                # transform to normal distribution
                p <- qnorm(p, lower.tail = FALSE) 
            
                weight <- weightdata$weight
           
                # we make probits from our pvalues
                # then we combine the probits
                # then we use the combined probits to make pvalues again
                # i.e. joined pvalues

                # formula 2.1 from "A Note on Combning Dependent Tests of Significance" Joachim Hartung
                cp <- pnorm(sum(weight * p)/sqrt((1-rho)*sum(weight^2)+rho*(sum(weight)^2)), lower.tail = FALSE)
                if(is.nan(cp) == TRUE) { cp <- 0 } #############################  workaround for the NaN problem (usually for artifact predictions)
                all_pvalues <- c(all_pvalues,cp)

                line<-paste(cp,paste(curr_cluster$d1, collapse=";"), sep=";") ## edit 2.0.5.1 // changed to d1
                line<-paste(line,";",sep="")
                final.out_lines_temp <- c(final.out_lines_temp,line)

            } else {

                indices_missing <- which(! weightdata$V1 %in% curr_cluster$id2) ## edit 2.0.5.1 // changed to id2
                indices_available <- which(weightdata$V1 %in% curr_cluster$id2) ## edit 2.0.5.1 // changed to id2
     
                p <- qnorm(p, lower.tail = FALSE)
                cNresult <- condNormal(p, means, covmatrix, indices_available, indices_missing)
                
                p_full <- rep(NA, orgcnt)
                p_full[indices_available] <- p
                na_indices <- which(is.na(p_full))

                for (i in 1:length(na_indices)) {
                    # this is the actual sampling before we're only building the parameters
                    # NAs sometimes produced here due to very low energy values in several targets in the cluster <= -100
                    p_full[na_indices[i]] <- rnorm(1, mean=cNresult$condMean[i], sd=sqrt(cNresult$condVar[i,i]))
                }
                # formula 2.1 from "A Note on Combning Dependent Tests of Significance" Joachim Hartung
                cp <- pnorm(sum(weight * p_full)/sqrt((1-rho)*sum(weight^2)+rho*(sum(weight)^2)), lower.tail = FALSE)
                if(is.nan(cp) == TRUE) { cp <- 0 } #############################  workaround for the NaN problem (usually for artifact predictions)
                all_pvalues <- c(all_pvalues,cp)
               
                line<-paste(cp,paste(curr_cluster$d1, collapse=";"), sep=";") ## edit 2.0.5.1 // changed to d1
                line<-paste(line,";",sep="")
                final.out_lines_temp <- c(final.out_lines_temp,line)
            }
        } 
        h <- hist(all_pvalues, plot=FALSE, breaks=100)
        temp_err <- sum((1-h$intensities)**2)
        # check for error reduction
        if(temp_err < error) {
            error <- temp_err
            optimal_rho <- rho
            final.out_lines <- final.out_lines_temp
            output <- list(optimal_rho,final.out_lines)
        }
    } 
    return(output)
}

prelim_rho_output<-compute_rho(seq(0,1,by=0.1), data, weightdata, means, covmatrix)
prelim_rho <- prelim_rho_output[[1]]

if (prelim_rho == 0) {
    prelim_rho <- 0.1
}

## edit 2.0.6 // added this part to prepare cop2 combination input lines with reduced runtime
if (prep_for_cop2 == 1) {
    coprarna2_prep_lines <- prelim_rho_output[[2]]
    write.table(file="CopraRNA2_prep.csv", coprarna2_prep_lines, row.names=FALSE, col.names=FALSE, quote=FALSE)   
}

## edit 2.0.6 // removed rhodevelopment.txt

if (prep_for_cop2 == 0) {
    final_rho_output<-compute_rho(seq((prelim_rho-0.1),(prelim_rho+0.1),by=0.01), data, weightdata, means, covmatrix)
    final_rho <- final_rho_output[[1]]
    final_out_lines <- final_rho_output[[2]]
    write.table(file="CopraRNA1_with_pvsample.csv",final_out_lines, row.names=FALSE, col.names=FALSE, quote=FALSE)
}

