Hartung <- function(p, lambda=rep(1,length(p)), kappa=0.2, alpha=0.10)
#
# This function applies the modified inverse normal method for the combination of dependent p-values.
#
# Arguments:
# p:         vector of p-values.
# lambda:    vector of weights. It must be of the same length of p.
# kappa:     adjustment parameter. It can be either a positive value (0.2 is the default value) or "formula". If 
#            k == "formula", then it is computed as in Hartung, p. 853.
# alpha:     level for the 1-alpha confidence interval for rho (0.10 is the default).
#
# Value: 
# 	     a list of class ``htest'' containing
#	     statistic: the Ht test statistic
#	     parameter: the number of combined tests (p-values)
#	     p.value: the combined test p-value
#	     conf.int: the confidence interval for the estimated correlation
#	     estimate: the estimated correlation
#	     null.value: the specified hypothesized value under the null
#	     alternative: a character string describing the alternative hypothesis
#	     method: a character srting indicating the type of combination test (only Z-test is implemented)
#	     data.name: a character string giving the name of the data set
#
# Reference:
# Hartung, J. (1999): "A note on combining dependent tests of significance",
#                     Biometrical Journal, 41(7), 849--855.
#
# Author:        Claudio Lupi
# This version:  August 22, 2010.
#
{
  t       <- qnorm(p)
  n       <- length(p)
  avt     <- sum(t)/n
  q       <- sum((t - avt)^2)/(n-1)                          # Hartung, eqn. (2.2)
  rhohat  <- 1 - q
  rhostar <- max(-1/(n-1), rhohat)                           # Hartung, p. 851
  if (kappa=="formula") kappa <- (1 + 1/(n-1) - rhostar)/10  # Hartung, p. 853
  if (kappa=="formula2") kappa <- (1 + 1/(n-1) - rhostar)/5  # Hartung, p. 853

  # Hartung inverse normal corrected. See eqn. (2.4)
  Ht <- sum(lambda*t)/sqrt(sum(lambda^2)+((sum(lambda))^2-sum(lambda^2))*(rhostar+kappa*sqrt(2/(n-1))*(1-rhostar)))
  lower <- 1 - (n-1)/qchisq(alpha/2, (n-1)) * q
  upper <- 1 - (n-1)/qchisq((1-alpha/2), (n-1)) * q          # Hartung, eqn. (2.3)
  
  output <- list(statistic=c("Ht"=Ht),
	  parameter=c("N"=n),
	  p.value=pnorm(Ht),
	  conf.int=c(lower, upper),
	  estimate=c("average estimated correlation"=as.vector(rhohat)),
	  null.value=("Ht"=0),
	  alternative="less",
	  method="modified inverse normal combination",
	  data.name=deparse(substitute(p)))

  class(output) <- "htest"

  return(output)
}

