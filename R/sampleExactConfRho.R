#' Sample from exact sampling distribution for pearson correlation
#'
#' Given a pearson correlation coefficient \code{r} estimated from \code{n} observations, randomly sample \code{nSamps} from the sampling distribution.  
#'
#' The standard \code{\link{cor.test}} simply uses Fisher's transformation to approximate the sampling distribution.  But this breaks down for small n or strong correlation, exactly where the uncertainty matters.  This function implements the exact sampling distribution pdf.  As such, it will obviously be slower, but more accurate.
#'
#' @param r Observed pearson correlation.
#' @param n Number of paired observations from which r was derived.
#' @param nSamps Number of random samples to draw.
#' @return \code{nSamps} from the sampling distribution of the correlation coefficient.
#' @importFrom distr AbscontDistribution 
#' @importFrom gsl hyperg_2F1
sampleExactConfRho = function(r,n,nSamps=1000){
  #Given r - measured correlation, and n - number of observations
  vu = n-1
  #See here for definition of sampling distribution https://en.wikipedia.org/wiki/Pearson_correlation_coefficient#Using_the_exact_distribution
  d = function(x){
    logp = log(vu) + log(vu-1) + lgamma(vu-1)-lgamma(vu+.5)-.5*log(2*pi)+(vu-1)/2*log(1-r*r) +(vu-2)/2*log(1-x*x) + (1-2*vu)/2*log(1-x*r) + log(hyperg_2F1(3/2,-1/2,vu+0.5,(1+r*x)/2))
    return(exp(logp))
  }
  dist =AbscontDistribution(d=d,low1=-1,up1=1,low=-1,up=1)
  return(r(dist)(nSamps))
}


