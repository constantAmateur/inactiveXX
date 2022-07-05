#' Infer heterozygous SNPs from scRNA
#'
#' Looks at locations of common SNPs in the human population (on chromosome X by default) and keeps as heterozygous any location where both alleles are detected in the data.  Optionally only consider those cells in \code{cellsToUse} when doing this. 
#'
#' If there are multiple 10X BAM files relating to the same individual, all should be passed to this function (via the \code{BAMs} parameter) at the same time to maximise the power to detect both alleles.  This depends on the alleleCount binary being present, or the vartrix binary if \code{useVartrix} is set.
#'
#' As well as returning locations of heterozygous SNPs, this function also stores the raw counts at each cell generated as part of the estimation procedure.  These can then be reused for estimation of the X inactivation state by passing the output of \code{\link{hetSNPsFromRNA}} to \code{\link{filterCountsX}} (see help for \code{\link{filterCountsX}} for more details).
#'
#' @param BAMs The 10X BAMs to use.
#' @param refGenome The reference genome used to map \code{BAMs}.
#' @param SNPs Which SNPs to consider.  Default are 1k genome SNPs on X with frequency greater than 1\%.  Note these have been lifted over to the GRCh38 reference.
#' @param cellsToUse Only look at cells with these barcodes.  Must be specified if using vartrix.
#' @param errRate When testing for presence of allele, use this error rate.
#' @param alpha Consider an allele detected at this significance threshold.
#' @param preserve Which metadata entries to preserve (if they're present in SNPs)
#' @param useVartrix Should we use vartrix to get allele counts?
#' @param ... Passed to \code{\link{getAllelicExpression}}
#' @return A GRanges object with the heterozygous SNPs.  Stored in the @metadata slot is also the cell specific counts for those snps to save future calculation.
#' @importFrom alleleIntegrator getAllelicExpression aggregateByClusters
#' @importFrom stats pbinom
#' @importFrom S4Vectors mcols
#' @importFrom GenomicRanges match
#' @export
hetSNPsFromRNA = function(BAMs,refGenome,SNPs=X1k,cellsToUse=NULL,errRate=0.05,alpha=0.05,preserve=c('REF','ALT','AF','strata','regionType','geneID','geneName'),useVartrix=FALSE,...){
  #Refine what to preserve
  toKeep = intersect(preserve,colnames(mcols(SNPs)))
  #Get the counts
  if(useVartrix){
    if(is.null(cellsToUse))
      stop("cellToUse is required for Vartix quantification of counts")
    cnts = vartrixCnts(SNPs,BAMs,cellsToUse,refGenome,...)
  }else{
    cnts = getAllelicExpression(SNPs,refGenome,BAMs,...)
  }
  #Filter them if we're going to.  Will be redundant for vartrix, but won't hurt.
  if(!is.null(cellsToUse))
    cnts = cnts[cnts$cellID %in% cellsToUse]
  #Group by SNP
  #tmp = aggregateByClusters(cnts,rep(1,length(cnts)))
  tmp = aggregateByClusters(cnts,rep(1,length(cnts)),assays=c('refCount','altCount','Tot'),toKeep)
  tmp$cellID=NULL
  #Do test for heterozygousness
  pValA = pbinom(tmp$refCount-1,tmp$Tot,errRate,lower.tail=FALSE)
  pValB = pbinom(tmp$altCount-1,tmp$Tot,errRate,lower.tail=FALSE)
  tmp$pVal = pmax(pValA,pValB)
  tmp$isHet = tmp$pVal<alpha
  #Do some final cleanup of object before return
  tmp$refCnt = tmp$refCount
  tmp$altCnt = tmp$altCount
  tmp$total = tmp$Tot
  tmp$BAF = tmp$altCnt/tmp$total
  #Drop ones that will conflict later if kept
  tmp$regionID = tmp$A = tmp$C = tmp$G = tmp$T = tmp$Tot = tmp$refCount = tmp$altCount = tmp$matCount = tmp$patCount = NULL
  #Store the cnts for later
  tmp@metadata$cntsAll = cnts
  tmp@metadata$cnts  = cnts[!is.na(match(cnts,tmp[tmp$isHet]))]
  return(tmp)
}
