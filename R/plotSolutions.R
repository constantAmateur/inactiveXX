#' Plot summary of solutions
#'
#' The EM algorithm used by \code{\link{inferInactiveX}} is rerun with many random starting points (1,000 by default).  This function plots a summary of these solutions, showing the log likelihood of each solution (y-axis) against the fraction of cells allocated to maternal (x-axis).  
#'
#' When the genotyping is not known beforehand the designation of maternal or paternal is arbitrary and you should expect the solutions to be evenly distributed either side of 0.5 on the x-axis.
#'
#' Jitter is added to points to allow the high density regions of the solution space to be visualised.
#'
#' @param fit The output of \code{\link{inferInactiveX}}.
#' @param jitter Fraction of plot range to jittery by in x,y direction, as a fraction of the dynamic range.
#' @return Nothing, but a plot is produced.
#' @importFrom stats rnorm
#' @importFrom graphics abline
#' @export
plotSolutions = function(fit,jitter=c(.01,.01)){
  if(length(fit$allFits)==1 && is.na(fit$allFits))
    stop("Fit didn't use EM, nothing to plot.")
  taus = sapply(fit$allFits, function(e) e$tau)
  Qs = sapply(fit$allFits,function(e) e$Q)
  d = data.frame(Q=Qs,tau=taus)
  plot(x=100*(taus+rnorm(length(taus),0,sd=abs(diff(range(taus)))*jitter[1])),
       y=Qs+rnorm(length(Qs),0,sd=abs(diff(range(Qs)))*jitter[2]),
       xlim=c(0,100),
       xlab='Xi skew',
       ylab='Log likelihood',
       main='Solution distribution',
       pch=19,
       cex=0.1)
  abline(v=50,col='darkred',lty=2)
  #gg = ggplot(d,aes(tau,Q)) +
  #  geom_hex() +
  #  geom_point(col = pointColour,size  = pointSize,alpha =  pointAlpha) +
  #  geom_vline(xintercept=0.5,col='red') +
  #  theme_bw() + 
  #  xlab('Fraction maternal') + 
  #  ylab('Log likelihood') +
  #  ggtitle('Solution distribution')
  #return(gg)
}
