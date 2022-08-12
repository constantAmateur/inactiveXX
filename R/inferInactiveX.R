#' Infer which copy of X is inactivated
#'
#' Given counts on the X chromosome, use EM to simultaneously estimate the genotype of the X chromosome and which version is inactivated in each cell.  If \code{cnts} includes columns \code{matCount} and \code{patCount}, the genotypes are assumed to have already been phased and rather than estimating the genotype, the genotype is assumed. 
#'
#' To ensure consistency of results, one cell (the \code{anchorCell}) is assumed to be in the maternally active state.  If not specified, this is picked to be the cell with the highest X chromosome coverage.
#'
#' The summary stats provided at the cell/SNP level essentially assume the other part of the data is correct and then measure how discrepent each cell or SNP is.  For example, the cell summary stats assumes all SNPs have been correctly phased and then measures how many reads are in conflict with this phasing across each cell.
#'
#' For samples with few cells and an extreme fraction of cells in one X-Inactivation state, it may not be possible to accurately estimate the genotype and X-Inactivation state of each cell.  However, in this situation the best fits across all random starts form show a trend where more extreme \code{tau} values (fraction of cells in one X-Inactivation state) always result in a better fit.  By contrast, where the best fit is well captured the top fits all cluster strongly around the same value.  
#'
#' This difference can be seen visually with \code{\link{plotSolutions}}, but the code will attempt to automatically detect when the fit should be more extreme.  This is done by calculating the difference between the best fit value of \code{tau} and the the top \code{tauDiffFrac*maxIter} best fits (divided by the best fit \code{tau}) and reporting when this difference exceeds \code{tauDiffThresh}.  The intuition here is that when the top fits all cluster around one value the average difference will be basically 0.
#'
#' @param cnts Usually the output of \code{\link{getAllelicExpression}}.  A GRanges of counts with refCount/altCount (matCount/patCount if pre-phased) with one row for each cell and SNP.
#' @param errRate Assumed all cause error rate.  10\% seems to be about the baseline, but can vary
#' @param logitCut Cells with probabilities (on logit) scale greater than this value are marked high confidence cells.
#' @param pCut Cells with p-value of excess discrepent counts under this value cannot be high confidence cells.
#' @param tauInit Initial guess for fraction of cells expressing maternal allele.
#' @param betaStart Initial beta to use for deterministic annealing.
#' @param betaFac How much to increase beta by in each pass through the deterministic annealing loop.
#' @param anchorCell ID of cell assumed to be maternally active.  If NULL, one with highest coverage used.
#' @param nStarts How many times to run the fit starting at random locations.
#' @param tol EM terminated when change in Q less than this value.
#' @param maxIter Terminate EM if more than this many iterations.
#' @param tauDiffFrac Sanity check performed using the top \code{tauDiffFrac} of fits.
#' @param tauDiffThresh Sanity check failed when the fractional difference in tau exceeds this value. 
#' @param nParallel How many threads?
#' @param verbose Be verbose?  Levels are 0 (off), 1 (outer loop), 2 (beta loop), 3 (EM loop)
#' @return A list with \code{states}, indicating which allele is active, \code{genotype} indicating if the reference allele is maternal for each SNP, \code{tau} the fraction of cells with maternal active, \code{allFits} which contains results for each of the \code{nStarts} random initiations, \code{cellSummary} which contains summary stats for each cell, \code{snpSummary} which contains summary stats for each SNP, and \code{dd} which contains the raw data.
#' @importFrom Matrix colSums
#' @importFrom stats dbinom pbinom aggregate p.adjust
#' @importFrom alleleIntegrator buildCountMatricies
#' @importFrom bettermc mclapply
#' @export
inferInactiveX = function(cnts,errRate=0.10,logitCut=3,pCut=0.2,tauInit=0.5,betaStart=.01,betaFac=1.3,anchorCell=NULL,nStarts=1000,tol=1e-6,maxIter=1000,tauDiffFrac=0.1,tauDiffThresh=0.05,nParallel=1,verbose=1){
  #Are the inputs phased?  If so, can just return the answer right now.
  if(!is.null(cnts$matCount) && !is.null(cnts$patCount)){
    #If they are, no need for EM.
    Xmat = buildCountMatricies(cnts)
    mCnts = colSums(Xmat$matCount)
    pCnts = colSums(Xmat$patCount)
    tot = mCnts+pCnts
    a = dbinom(mCnts,tot,1-errRate,log=TRUE)
    b = dbinom(pCnts,tot,1-errRate,log=TRUE)
    states = 1/(1+exp(b-a))
    tau = mean(states)
    #Store genotype if in standard format
    if(!is.null(cnts$altIsMum)){
      gtype = setNames(as.numeric(cnts$altIsMum[!duplicated(names(cnts))]),names(cnts)[!duplicated(names(cnts))])
    }else{
      gtype=NA
    }
    #Store dd 
    Xmat = buildCountMatricies(cnts,assays=c('refCount','altCount'))
    tMat = as(Xmat[[1]],'dgTMatrix')
    tPat = as(Xmat[[2]],'dgTMatrix')
    tMat= data.frame(SNP=rownames(tMat)[tMat@i+1],
                     cell=colnames(tMat)[tMat@j+1],
                     refCount = tMat@x
                     )
    tPat= data.frame(SNP=rownames(tPat)[tPat@i+1],
                     cell=colnames(tPat)[tPat@j+1],
                     altCount = tPat@x
                     )
    dd = merge(tMat,tPat,by=c('SNP','cell'))
    dd$tot = dd$refCount + dd$altCount
    dd = dd[dd$tot>0,]
    fin = list(tau=tau,states=states,genotype=gtype,Q=NA,allFits=NA,dd=dd)
  }else{
    #Prepare data
    Xmat = buildCountMatricies(cnts,assays=c('refCount','altCount'))
    cells = colnames(Xmat[[1]])
    snps = rownames(Xmat[[1]])
    #Pick anchor cell
    if(is.null(anchorCell)){
      Xmat$tot = Xmat[[2]] + Xmat[[1]]
      anchorCell = colnames(Xmat$tot)[which.max(colSums(Xmat$tot>0))]
    }
    if(!anchorCell %in% cells)
      stop("Anchor cell not found.")
    tMat = as(Xmat[[1]],'dgTMatrix')
    tPat = as(Xmat[[2]],'dgTMatrix')
    tMat= data.frame(SNP=rownames(tMat)[tMat@i+1],
                     cell=colnames(tMat)[tMat@j+1],
                     refCount = tMat@x
                     )
    tPat= data.frame(SNP=rownames(tPat)[tPat@i+1],
                     cell=colnames(tPat)[tPat@j+1],
                     altCount = tPat@x
                     )
    dd = merge(tMat,tPat,by=c('SNP','cell'))
    dd$tot = dd$refCount + dd$altCount
    dd = dd[dd$tot>0,]
    if(verbose)
      message(sprintf("Running EM fit with %d random starts using %d threads",nStarts,nParallel))
    #Run model n times and pick best
    out = mclapply(seq(nStarts),function(e) {
                     tmp = deterministicAnnealing(dd,betaStart,betaFac,tauInit,anchorCell,verbose,errRate=errRate,tol=tol,maxIter=maxIter)
                     if(verbose>1)
                       message(sprintf('[%d of %d] Fit found with tau = %.02f and Q=%g.',e,nStarts,tmp$tau,tmp$Q))
                     return(tmp)},
                     mc.allow.fatal=TRUE,
                     mc.silent=verbose<1,
                     mc.progress= interactive() && verbose,
                     mc.retry=10,
                     mc.cores=nParallel)
    #Check if it's likely there is a more extreme solution
    taus = sapply(out,function(e) e$tau)
    Qs = sapply(out,function(e) e$Q)
    o = order(Qs)
    taus = taus[o]
    Qs = Qs[o]
    tau = tail(taus,n=1)
    tauDiff = pmin(1-taus,taus)-tau
    #Get the top N and calculate the relative difference from the best value
    w = tail(seq_along(taus),n=ceiling(length(taus)*tauDiffFrac))
    if((mean(tauDiff[w])/tau) > tauDiffThresh)
      warning(sprintf('Top %d fits highly discepent.  It is likely that a more extreme X-Inactivation fraction is true, but there is insufficient data to estimate it.',length(w)))




    maxIter*.1

    #Pick solution with best likelihood
    best = out[[which.max(sapply(out,function(e) e$Q))]]
    fin = list(tau=best$tau,states=best$states,genotype=best$lambdas,dd=dd,allFits=out)
  }
  #Add in some more details to dd
  #fin$dd$isMat = fin$states[fin$dd$cell]>0.5
  fin$dd$matCount = ifelse(fin$genotype[fin$dd$SNP],fin$dd$refCount,fin$dd$altCount)
  fin$dd$patCount = ifelse(!fin$genotype[fin$dd$SNP],fin$dd$refCount,fin$dd$altCount)
  #fin$dd$mainCount = ifelse(fin$dd$isMat,fin$dd$patCount,fin$dd$matCount)
  #fin$dd$offCount = ifelse(!fin$dd$isMat,fin$dd$patCount,fin$dd$matCount)
  #Call high confidence states
  tmp = aggregate(cbind(matCount,patCount,tot) ~ cell,data=fin$dd,FUN=sum)
  tmp$offCount = pmin(tmp$matCount,tmp$patCount)
  tmp$matFrac = tmp$matCount/tmp$tot
  tmp$badFrac = tmp$offCount/tmp$tot
  tmp = tmp[match(names(fin$states),tmp$cell),]
  #Convert to logit scale for convenience
  tmp$states = log(fin$states/(1-fin$states))
  tmp$highConfCall = abs(tmp$states)>logitCut
  #Is it inconsistent with the specified error rate?
  tmp$pVal = pbinom(tmp$offCount-1,tmp$tot,errRate,lower.tail=FALSE)
  tmp$qVal = NA
  tmp$qVal[tmp$highConfCall] = p.adjust(tmp$pVal[tmp$highConfCall],method='BH')
  #The qVal filter is too tolerant of letting through crap
  tmp$highConfCall = tmp$highConfCall & tmp$pVal>=pCut
  rownames(tmp) = tmp$cell
  fin$cellSummary = tmp
  #Now propogate the high confidence calls back to the raw data
  fin$dd$isMat = ifelse(tmp[fin$dd$cell,'highConfCall'],fin$states[fin$dd$cell]>0.5,NA)
  fin$dd$mainCount = ifelse(fin$dd$isMat,fin$dd$patCount,fin$dd$matCount)
  fin$dd$offCount = ifelse(!fin$dd$isMat,fin$dd$patCount,fin$dd$matCount)
  #Now do the same for SNPs
  tt1 = aggregate(cbind(refCount,altCount,matCount,patCount,tot) ~ SNP,data=fin$dd,FUN=sum)
  tt2 = aggregate(cbind(mainCount,offCount,tot) ~ SNP,data=fin$dd,FUN=sum,na.rm=TRUE)
  colnames(tt2)[4]='totOff'
  tmp = merge(tt1,tt2,by='SNP',all=TRUE)
  rownames(tmp) = tmp$SNP
  tmp$badFrac = tmp$offCount/tmp$totOff
  tmp$matFrac = tmp$matCount/tmp$tot
  tmp$pVal = pbinom(tmp$offCount-1,tmp$totOff,errRate,lower.tail=FALSE)
  tmp$qVal = p.adjust(tmp$pVal,method='BH')
  tmp$highConfCall = FALSE
  tmp$highConfCall[which(tmp$qVal>=pCut)]=TRUE
  tmp$genotype = fin$genotype[tmp$SNP]
  fin$snpSummary = tmp
  class(fin) = c('list','inactiveXX')
  return(fin)
}
