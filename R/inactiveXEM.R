#' Inner EM loop
#'
#' Runs the EM loop to estimate the X chromosome genotype and inactivation state.  Not intended to be used directly.
#'
#' @inheritParams deterministicAnnealing
#' @param beta Deterministic annealing parameter for this fit. beta=1 solves the unperterbed likelihood function.
#' @param lambdas SNP state.  Maternal allele is alt when 1 and ref when 0.
#' @param states Cell states.  Maternal genotype when 1, (i.e. paternal X inactivated)
#' @param errRate The error rate to asume.
#' @param tol EM terminated when change in Q less than this value.
#' @param maxIter EM terminated after this many iterations, regardless of change in Q.
#' @return List of states, genotype, tau, iteration number, Q, dQ
#' @importFrom stats dbinom setNames
inactiveXEM = function(dd,beta,tauInit,lambdas,states,errRate,tol,maxIter,verbose=FALSE){
  #Model variable definitions:
  #lambda_i = The maternal allele is the alternate at SNP i when this is 1 and reference when 0
  #state_j = Cell j expresses the maternal genotype when this is 1 (i.e., paternal copy of X inactivated)
  #tau = Fraction of cells that are expressing the maternal genotype
  #Initialise
  if(is.null(lambdas)){
    lambdas = unique(dd$SNP)
    lambdas = setNames(sample(c(0,1),length(lambdas),replace=TRUE),lambdas)
  }
  if(is.null(states)){
    states = unique(dd$cell)
    states = setNames(sample(c(0,1),length(states),replace=TRUE),states)
  }
  cells = names(states)
  snps = names(lambdas)
  tau=tauInit
  Q = -Inf
  i=0
  while(TRUE){
    #E-step
    oldStates = states
    probs = lambdas[dd$SNP]*(2*errRate-1) + (1-errRate)
    #Calculate in log-space for numerical stability
    a = dbinom(dd$refCount,dd$tot,probs,log=TRUE) #Probability that maternal genotype is expressed
    b = dbinom(dd$altCount,dd$tot,probs,log=TRUE) #Probability that paternal genotype is expressed
    states = sapply(split(beta*(b-a),dd$cell),sum)
    #Back to linear space
    states = tau**beta/(tau**beta+(1-tau)**beta*exp(states))
    #Construct Q explicitly
    Qnew = states[dd$cell]*(log(tau)+dbinom(dd$refCount,dd$tot,probs,log=TRUE))+(1-states[dd$cell])*(log(1-tau)+dbinom(dd$altCount,dd$tot,probs,log=TRUE))
    Qnew = sum(Qnew)
    #M-step
    oldLamb = lambdas
    tau = mean(states)
    lambdas = (dd$refCount-dd$altCount)*(1-2*states[dd$cell])
    #Split and sum by snp
    lambdas = sapply(split(lambdas,dd$SNP),sum)
    lambdas = ifelse(lambdas>0,1,0)
    #Check termination condition
    i=i+1
    dQ = Qnew-Q
    Q=Qnew
    #Update us on the goings on
    if(verbose>2)
      message(sprintf("    [beta=%.02f i=%.03d] Q=%g, dQ = %g, %d SNPs were swapped, assignments moved by %g, and new tau is %g",beta,i,Q,dQ,sum(abs(oldLamb-lambdas)),sum(abs(oldStates-states)),tau))
    if(dQ< tol | i==maxIter)
      break
  }
  converged = i<maxIter
  if(i==maxIter)
    warning(sprintf("Run with beta=%g failed to converge in %d iterations.",beta,maxIter))
  return(list(states=states,lambdas=lambdas,tau=tau,Q=Q,dQ=dQ,i=i,converged=converged))
}
