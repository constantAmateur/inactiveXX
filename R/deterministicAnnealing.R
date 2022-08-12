#' Deterministic annealing expectation maximisation
#'
#' Uses the deterministic annealing approach to wrap around \code{\link{inactiveXEM}}.  This smooths the likelihood function strongly to begin with, then successively reduces the smoothing in increments until it is zero.  This approach should allow the EM algorithm a better chance of converging to a global, rather than local, maximum. 
#'
#' Not intended for direct use, call \code{\link{inferInactiveX}} instead.
#'
#' @inheritParams inferInactiveX
#' @param dd A data.frame where rows are observations at a particular cell/SNP combination.  Columns must be \code{refCount}, \code{altCount}, \code{tot}, \code{cell}, and \code{SNP}.
#' @param ... Extra parameters passed to \code{\link{inactiveXEM}}
#' @return Fit summary for 
deterministicAnnealing = function(dd,betaStart,betaFac,tauInit,anchorCell,verbose=FALSE,...){
  #These will be randomly initialised by inactiveXEM below
  lambdas = NULL
  states = NULL
  tau = tauInit
  #Do settling in deterministic annealing to get good guess of initial parameters
  betas = betaStart*betaFac**seq(0,floor(-log(betaStart)/log(betaFac)))
  if(max(betas)!=1)
    betas = c(betas,1)
  i=0
  for(beta in betas){
    i=i+1
    if(verbose>2)
      message(sprintf('  [%d of %d] Running deterministic annealing with beta=%g.',i,length(betas),beta))
    tmp = inactiveXEM(dd,beta,tau,lambdas,states,...,verbose=verbose)
    lambdas = tmp$lambdas
    states = tmp$states
    tau = tmp$tau
  }
  #Make sure anchor cell is in state 1
  #if(states[anchorCell]<0.5){
  #  #Mirror around 0.5
  #  tmp$tau = 1-tau
  #  tmp$states = 1-states
  #  tmp$lambdas =1-lambdas
  #}
  return(tmp)
}


