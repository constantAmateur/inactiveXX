#' Estimates effective population bottleneck size
#'
#' Given Xi data from the same group of cells (cell type, cell state, tissue, etc.), estimate the size of the population bottleneck that would be needed to produce the observed spread in Xi skew across individuals.  This value will represent the smallest population bottleneck in the ancestory of the cell group investigated, which for many cells will be the "founder cell number"; the number of cells present when X-inactivation was established giving rise to this lineage.
#'
#' The estimation is performed using the model that the population size is related to the variance of the Xi skew across individuals as N = 1/(4*var), which is derived from the assumption that Xi skew for an individual is binomially distributed across the founder cells (with p=0.5).  This function takes into account two sources of uncertainty in estimating the effective population bottleneck size N.  Firstly, in each individual, the measurement of Xi skew has uncertainty, which decreases as the number of cells in each individual increases.  Secondly, the variance of the distribution of Xi skew from a finite number of individuals is uncertain.  Typically, most of the unceratainty is driven by the number of individuals, not the number of cells in each individual.
#' 
#' To estimate the uncertainty on the measured Xi skew in each individual, a beta prior with parameters alpha and beta are used.  The default values imply a completely uniform prior.
#'
#' By default, this function will calculate the population size for one group of cells.  Multiple groups can be estimated simultaneously if \code{cellsToUse} is given as a list of cellIDs, where then one estimate will be produced per list entry and named by the names of the \code{cellsToUse} list.
#'
#' @inheritParams cellGroupSummary
#' @param priorAlpha alpha pramater for beta distribution prior on individual Xi skew.
#' @param priorBeta beta pramater for beta distribution prior on individual Xi skew.
#' @param nSamps Number of random samples used to estimate uncertainty.
#' @param distQuants Which quantiles of the posterior distribution to provide?
#' @return A data.frame containing the estimate of the populaion size with uncertainty.
#' @importFrom stats rbeta rchisq var
#' @export
estPopSize = function(fits,cellsToUse=NULL,highConfOnly=FALSE,priorAlpha=1,priorBeta=1,nSamps=1000,distQuants=c(0.01,0.05,0.1,0.9,0.95,0.99)){
  dd = cellGroupSummary(fits,cellsToUse,highConfOnly)
  dd$postA = priorAlpha + dd$nMat
  dd$postB = priorAlpha + dd$nPat
  #Now simulate from each of them
  samps = matrix(rbeta(nSamps*nrow(dd),dd$postA,dd$postB),ncol=nSamps)
  #For each cell type, calculate the summary
  out=list()
  outN = list()
  for(cellGroup in unique(dd$cellGroup)){
    w = which(dd$cellGroup==cellGroup)
    if(is.na(cellGroup))
      w = which(is.na(dd$cellGroup))
    #For each random sample, estimate variance then incorporate the finite individuals variability
    sig2 = apply(samps[w,],2,var)
    sig2 = (length(w)-1)*sig2/rchisq(nSamps,length(w)-1)
    #Finally, convert these to Ns
    Ns =1/(4*sig2)
    #Prepare output summary
    t1 = data.frame(cellGroup = cellGroup,pointEst = 1/(4*var(dd$matFrac[w])))
    t2 = quantile(Ns,distQuants)
    names(t2) = paste0('q',gsub('%','',names(t2)))
    t1[,names(t2)]=t2
    out[[length(out)+1]] = t1
    if(is.na(cellGroup)){
      outN[[length(outN)+1]] = Ns
    }else{
      outN[[cellGroup]] = Ns
    }
  }
  out = do.call(rbind,out)
  attributes(out)$rawSamples = outN
  return(out)
}

