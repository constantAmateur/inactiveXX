#' Plot summary of solutions
#'
#' The EM algorithm used by \code{\link{inferInactiveX}} is rerun with many random starting points (1,000 by default).  This function plots a summary of these solutions, showing the log likelihood of each solution (y-axis) against the fraction of cells allocated to maternal (x-axis).  
#'
#' When the genotyping is not known beforehand the designation of maternal or paternal is arbitrary and you should expect the solutions to be evenly distributed either side of 0.5 on the x-axis.
#'
#' @param fit The output of \code{\link{inferInactiveX}}.
#' @param pointColour Colour of points.
#' @param pointSize Size of points.
#' @param pointAlpha Transparency of points.
#' @param ... Passed to \code{\link{geom_hex}}.
#' @return A ggplot object containing the plot.
#' @importFrom ggplot2 ggplot geom_hex geom_point theme_bw xlab ylab ggtitle aes
#' @export
plotSolutions = function(fit,pointColour='darkgrey',pointSize=1.0,pointAlpha=1,...){
  if(length(fit$allFits)==1 && is.na(fit$allFits))
    stop("Fit didn't use EM, nothing to plot.")
  taus = sapply(fit$allFits, function(e) e$tau)
  Qs = sapply(fit$allFits,function(e) e$Q)
  d = data.frame(Q=Qs,tau=taus)
  gg = ggplot(d,aes(tau,Q)) +
    geom_hex() +
    geom_point(col = pointColour,size  = pointSize,alpha =  pointAlpha) +
    geom_vline(xintercept=0.5,col='red') +
    theme_bw() + 
    xlab('Fraction maternal') + 
    ylab('Log likelihood') +
    ggtitle('Solution distribution')
  return(gg)
}
