#' Estimate correlation of Xi skew between cell groups
#'
#' Given Xi skew estimates from a range of individuals, calculate the degree to which different groups of cells Xi skew are correlated.  This correlation can either be absolute, or relative to some estimate of the founder population skew (usually the mean or median Xi skew across cell types).
#'
#' As with \code{\link{estPopSize}}, there are two sources of uncertainty in these estimates, which depend on the number of cells of a given type within individuals, and the total number of individuals.  Both are taken into account in this estimation.
#'
#' Correlation coefficients are calculated for all pairwise combinations of cell types in \code{cellsToUse} that are listed in \code{groupsToCompare}.  
#'
#' When calculating correlation based on raw Xi skew, there is no need to include cell groupings other than the ones you are interested in in \code{cellsToUse}.  However, when \code{relativeTo} is set to something other than 'none', other cells types are needed.  This is because the comparison point is itself a derived quantity from the Xi skews, and as such should be derived independently from each random sample from the posterior. 
#'
#' Options for code{relativeTo} are: 'fixed', which correlates the Xi skew relative to \code{fixedVal} (0 by default), 'median', which correlates the difference of each Xi skew from the median value across all cell groupings, and 'tau', which correlates the difference from the average across all cells.
#'
#' When \code{exactConf} is \code{TRUE}, the exact sampling distribution for the correlation coefficient is used to estimate uncertainty due to number of individuals.  This is more accurate, particularly for low numbers of individuals, but slow (usually not so slow as to be a major issue).  When set to \code{FALSE}, the Fisher transformation is instead used.
#'
#' @inheritParams estPopSize
#' @param cellsToUse List containing at least two entries, named by cell types, giving the cellIDs for the cells of that type.
#' @param groupsToCompare Which of the cell groupings defined by `cellsToUse` to calculate correlation on.
#' @param relativeTo Should we calculate correlation on raw Xi skew, or Xi skew relative to something?
#' @param fixedVal Value to measure Xi skew differences from.  Ignored unless `relativeTo='fixed'`.
#' @param exactConf Should the exact distribution be used to calculate uncertainty (slower, but more accurate)?
#' @param verbose Level of verbosity, higher numbers is more.
#' @param nParallel How many threads in parallel
#' @param mcParams Overwrite default parameters controlling the call to \code{\link{mclapply}}.  Only used when \code{exactConf=TRUE} and \code{nParallel>1}
#' @return A data.frame giving point estimates and confidence distributions for each pairwise comparison.
#' @importFrom bettermc mclapply
#' @importFrom stats rbeta cor median
#' @export
estCorrelation = function(fits,cellsToUse,groupsToCompare = names(cellsToUse),relativeTo=c('fixed','median','tau'),fixedVal=0,exactConf=FALSE,highConfOnly=FALSE,priorAlpha=1,priorBeta=1,nSamps=1000,distQuants=c(0.01,0.05,0.1,0.9,0.95,0.99),verbose=0,nParallel=1,mcParams=list()){
  relativeTo = match.arg(relativeTo)
  if(!inherits(cellsToUse,'list') || length(cellsToUse)<2)
    stop('cellsToUse must be a list of at least 2 cell groupings')
  dd = cellGroupSummary(fits,cellsToUse,highConfOnly)
  dd$postA = priorAlpha + dd$nMat
  dd$postB = priorAlpha + dd$nPat
  #Now simulate from each of them
  samps = matrix(rbeta(nSamps*nrow(dd),dd$postA,dd$postB),ncol=nSamps)
  #Work out comparison point
  if(relativeTo=='fixed'){
    compVal = matrix(fixedVal,nrow=nrow(dd),ncol=nSamps)
  }
  if(relativeTo=='median'){
    compVal = do.call(rbind,lapply(split(seq(nrow(samps)),dd$individual),function(e) apply(samps[e,,drop=FALSE],2,median,na.rm=TRUE)))
  }
  if(relativeTo=='tau'){
    compVal = matrix(sapply(fits,function(e) e$tau)[dd$individual],nrow=nrow(dd),ncol=nSamps,byrow = FALSE)
  }
  #Adjust for both the sampling and point estimates
  relVal = samps-compVal
  dd$relFrac = switch(relativeTo,
                      'fixed' = dd$matFrac - fixedVal,
                      'median' = dd$matFrac - sapply(split(dd$matFrac,dd$individual),median,na.rm=TRUE)[dd$individual],
                      'tau' = dd$matFrac - sapply(fits,function(e) e$tau)[dd$individual])
  #Keep just the ones we want to compare
  w = which(dd$cellGroup %in% groupsToCompare)
  if(length(w)==0)
    stop("No comparison groups found.")
  dd = dd[w,]
  samps = samps[w,]
  relVal = relVal[w,]
  #For each random sample, calculate the correlation matrix
  nn = unique(dd$individual)
  cc = groupsToCompare
  idx = cbind(match(dd$cellGroup,cc),match(dd$individual,nn))
  corVals = list()
  for(i in seq(nSamps)){
    tt = matrix(NA,nrow=length(cc),ncol=length(nn),dimnames=list(cc,nn))
    tt[idx] = relVal[,i]
    corVals[[i]] = cor(t(tt),use='pairwise.complete.obs')
  }
  #Do point estimate as well
  tt = matrix(NA,nrow=length(cc),ncol=length(nn),dimnames=list(cc,nn))
  tt[idx] = dd$relFrac
  pe = cor(t(tt),use='pairwise.complete.obs')
  #We have our correlation values, up to variability due to finite cell numbers, the rest is to add variablity due to finite indivdual numbers.
  #Get number of datapoints for mutual pairs
  nPairs = matrix(NA,nrow=length(cc),ncol=length(cc))
  rownames(nPairs)=colnames(nPairs)=cc
  for(a in cc){
    for(b in cc){
      nPairs[match(a,cc),match(b,cc)] = sum(table(dd[dd$cellGroup %in% c(a,b),'individual'])==2)
    }
  }
  #Work out which indicies we need to run on.  As symmetric and diagonals trivial, only need to operate on the upper triangle.  Those without 4 obs can't have confidence interval, so ignore.
  w = which(nPairs>=4,arr.ind=TRUE)
  w = w[w[,2]>w[,1],,drop=FALSE]
  corVals = lapply(corVals,function(e) {tmp=e;tmp[seq_along(tmp)]=NA;tmp[w]=e[w];tmp})
  if(exactConf){
    #Use exact distribution
    defParams = list(X = corVals,
                        FUN = function(e){
                          for(i in seq(nrow(w))){
                            e[w[i,,drop=FALSE]] = sampleExactConfRho(e[w[i,,drop=FALSE]],nPairs[w[i,,drop=FALSE]],1)
                          }
                          return(e)},
                        mc.allow.error=TRUE,
                        mc.allow.fatal=TRUE,
                        mc.silent=verbose<1,
                        mc.retry.silent = verbose<1,
                        mc.dump.frames='no',
                        mc.shm.ipc=FALSE,
                        mc.share.copy=FALSE,
                        mc.share.vectors=FALSE,
                        mc.progress= interactive() && verbose,
                        mc.preschedule=TRUE,
                        mc.retry=10,
                        mc.cores=nParallel)
    for(nom in names(mcParams))
      defParams[[nom]] = mcParams[[nom]]
    corVals = do.call(mclapply,defParams)
  }else{
    #Use fisher transformation
    corVals = lapply(corVals,function(e) {
                       #Sample from Fisher transformed distribution
                       t1 = rnorm(nrow(w),mean=atanh(e[w]),sd=1/sqrt(nPairs[w]-3))
                       #Convert back to original space and save
                       e[w] = tanh(t1)
                       return(e)})
  }
  #Make back into symmetric
  corVals = lapply(corVals,function(e) {e[lower.tri(e)]=t(e)[lower.tri(e)];e})
  #Flatten and return
  out = data.frame(cellGroupA = rownames(pe)[w[,1]],
                   cellGroupB = colnames(pe)[w[,2]],
                   pointEst = pe[w])
  t1 = do.call(rbind,lapply(corVals,function(e) e[w]))
  for(q in distQuants)
    out[,paste0('q',100*q)] = apply(t1,2,quantile,q)
  #Save raw samples
  attributes(out)$rawSamples = corVals
  return(out)
}
