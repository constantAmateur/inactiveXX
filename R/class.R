#' Print summary of fit object
#'
#' Summary of features of fit returned by \code{\link{inferInactiveX}}.
#'
#' @param x Object to print.
#' @param ... Unused.
#' @return Nothing.  Print a message.
print.inactiveXX = function(x,...){
  cat("An object of class inactiveXX\n")
  cat(sprintf(" contains EM fit of X-inactivation state from %d random starts\n",length(x$allFits)))
  cat(sprintf(" %.02f%% cells express maternal X\n",x$tau*100))
  cat(sprintf(" %d cells with %d SNPs\n",length(x$states),length(x$genotype)))
  cat(sprintf(" %d cells (%.02f%%) called with high confidence\n",sum(x$cellSummary$highConfCall),sum(x$cellSummary$highConfCall)/length(x$states)*100))
  cat(sprintf(" %d SNPs (%.02f%%) called with high confidence\n",sum(x$snpSummary$highConfCall),sum(x$snpSummary$highConfCall)/length(x$genotype)*100))
  ###message(sprintf("Object of class inactiveXX, containing EM fit of X-inactivation state from %d random starts\n %d cells with %d SNPs\n %.02f%% cells expressing maternal X\n %d cells (%.02f%%) called with high confidence\n %d SNPs (%.02f%%) called with high confidence",
  ###                length(fit$allFits),
  ###                length(fit$states),
  ###                length(fit$genotype),
  ###                fit$tau*100,
  ###                sum(fit$cellSummary$highConfCall),
  ###                sum(fit$cellSummary$highConfCall)/length(fit$states)*100,
  ###                sum(fit$snpSummary$highConfCall),
  ###                sum(fit$snpSummary$highConfCall)/length(fit$genotype)*100
  ###                ))
}

