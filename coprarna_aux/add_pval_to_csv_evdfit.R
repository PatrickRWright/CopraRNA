## edit 2.0.4
# this script replaces the old
# add_pval_to_csv_evdfit.pl and 
# makes the Statistics::R dependency obsolete

# library(evir) ## edit 2.0.5 removed this dependency
# now sourcing gev and pgev locally here

# gev and pgev functions
# from the evir R package
# https://CRAN.R-project.org/package=evir
# ## edit 2.0.5 added this source locally 
# ## to be independent of evir

# start insert new functions

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

# end insert new functions

# read options ## edit 2.0.5.1
options <- read.table("CopraRNA_option_file.txt", sep=":")
region <- as.character(options$V2[4])

files_opt <- list.files(pattern="_opt.intarna.csv")
files_subopt <- list.files(pattern="_subopt.intarna.csv")

# sort to be 100% sure that the vectors 
# have the file for the same organism at
# the same place
files_opt <- sort(files_opt)
files_subopt <- sort(files_subopt)

for (i in 1:length(files_opt)) {

    curr_opt_file <- files_opt[i]
    curr_subopt_file <- files_subopt[i]

    # exception // I don't expect this to ever happen
    if (substr(curr_subopt_file, 1,20)==substr(curr_opt_file, 1,20)) { 
    } else { write("Files don't match in add_pval_to_csv_evdfit.R", file="err.log", append=T) }

    # opt
    data_opt <- read.table(curr_opt_file, sep=";", header=TRUE)
    energy_opt <- data_opt$E ## edit 2.0.5.1 // changed to E because headers in IntaRNA 2 are different
    # if the predicton is running on the full length transcript then use normalized energy
    if (region == "cds") {
        energy_opt <- data_opt$E_norm ## edit 2.0.5.1 // added E_norm for cds option
    }
    energy_opt <- energy_opt*(-1)
    # subopt
    data_subopt <- read.table(curr_subopt_file, sep=";", header=TRUE)
    energy_subopt <- data_subopt$E ## edit 2.0.5.1 changed to E because headers in IntaRNA 2 are different
    # if the predicton is running on the full length transcript then use normalized energy
    if (region == "cds") {
       energy_subopt <- data_subopt$E_norm ## edit 2.0.5.1 // added E_norm for cds option
    }
    energy_subopt <- energy_subopt*(-1)
    # fit evd and get parameters
    gevenergies <- gev(energy_opt) # fit
    xi <- gevenergies$par.ests[1] # read chi from fitting
    sigma <- gevenergies$par.ests[2] # read sigma
    mu <- gevenergies$par.ests[3] # read mu
    # get pvalues
    opt_pvals <- 1-(pgev(energy_opt, xi, mu, sigma))
    opt_pvals <- round(opt_pvals, digits=7)
    subopt_pvals <- 1-(pgev(energy_subopt, xi, mu, sigma))
    subopt_pvals <- round(subopt_pvals, digits=7)
    # add pvals to data frame
    data_opt <- cbind(data_opt, opt_pvals)
    data_subopt <- cbind(data_subopt, subopt_pvals)
    # sort by pvalue
    data_opt <- data_opt[order(data_opt$opt_pvals),]
    data_subopt <- data_subopt[order(data_subopt$subopt_pvals),]
    # change p-value column name 
    colnames(data_opt)[36]<-"p-value"    ## edit 2.0.5.1 // changed to 36 because there are now more columns
    colnames(data_subopt)[36]<-"p-value" ## edit 2.0.5.1 // changed to 36 because there are now more columns
    # write output
    # this overwrites the input files
    write.table(data_opt, file=curr_opt_file, sep=";", quote=F, row.names=F)
    write.table(data_subopt, file=curr_subopt_file, sep=";", quote=F, row.names=F)
}


