#' Filter X counts
#'
#' Given allele and cell resolved, counts at heterozygous SNPs on X, filter them down to just the useful ones.  \code{errRate} and \code{pCut} have both been calibrated quite precisely, don't change them unless you know what you're doing.
#'
#' In some cases, it may be useful to set \code{dropBiallelic} to \code{FALSE}.  The logic behind dropping biallelic SNPs (that is those for which both alleles are detected in the same cell) is sound.  However, for SNPs with high coverage, this can occur just due to errors a non-trivial fraction of the time.  A simple statistical test is used to try and reduce this, but will not be perfect.  While discarding SNPs with evidence of both alleles is still the consernative thing to do, it may come at the cost of throwing out too much data.
#'
#' Usually the output of \code{\link{hetSNPsFromRNA}} will be the input for this function (to \code{cnts}).  However, if for some reason you have calculated cell specific allele counts (usually because you have supplied a list of known heterozygous SNPs) from running \code{\link[alleleIntegrator]{getAllelicExpression}} or similar, that can also be supplied as input.
#'
#' Filtering is done on both the evolutionary strata (excluding those known to escape X-inactivation) and the SNPs location relative to genes (e.g. intronic, exonic, or intergenic).  To do this, it is assumed that a column called \code{regionType} will encode the SNPs relationship to genes and \code{strata} will encode the evolutionary strata.  This are pre-calculated when running \code{\link{hetSNPsFromRNA}} with default parameters, but if other inputs are used they will not be present.  In this case the user may either manually create these collumns, or (assuming \code{autoSet==TRUE}) allow them to be automatically created assuming a GRCh38 reference from the \code{\link{X1k}} data object.
#'
#' @param cnts The counts to filter.  Usually the output of \code{\link{hetSNPsFromRNA}}, \code{\link[alleleIntegrator]{getAllelicExpression}}, or \code{\link{vartrixCnts}}.
#' @param cellsToUse Filter out any cells not in this list.  If NULL, no filtering on cell barcdodes.
#' @param minObs Drop any SNP that is not observed in at least this many cells.
#' @param dropBiallelic Drop any SNP that has evidence of escaping X-inactivation.
#' @param dropStrata Drop those in these stratas.  Usually these mark regions where the chance of escape of X-Inactivation is high.
#' @param regionsToUse Drop regions not off this type.
#' @param errRate Error rate to assume for finding biallelic SNPs
#' @param pCut P-value cut to use for finding biallelic SNPs.
#' @param autoSet If strata and/or regions not set, set them automatically using a GRCh38 reference.
#' @return Filtered counts.
#' @importFrom stats pbinom
#' @importFrom GenomeInfoDb seqnames dropSeqlevels seqlevels
#' @importFrom S4Vectors mcols queryHits subjectHits
#' @importFrom GenomicRanges findOverlaps
#' @export
filterCountsX = function(cnts,cellsToUse=NULL,minObs=2,dropBiallelic=TRUE,dropStrata=c('PAR1','PAR2','strata2','strata3'),regionsToUse=c('Exonic','Intronic'),errRate=0.02,pCut=0.1,autoSet=TRUE){
  #Check if it's the hetSNPsFromRNA option 
  if(!is.null(cnts@metadata$cnts))
    cnts = cnts@metadata$cnts
  #Drop cells
  if(!is.null(cellsToUse))
    cnts = cnts[cnts$cellID %in% cellsToUse]
  #Drop regions
  if(is.null(cnts$regionType)){
    if(autoSet){
      message('regionType not found, setting assuming data has been mapped to GRCh38')
      cnts$regionType = 'InterGenic'
      t1 = X1k@metadata$gtf
      regions = list(Genic=t1[as.character(t1$type)=='gene',],
                     Intronic = t1[as.character(t1$type)=='transcript',],
                     Exonic = t1[as.character(t1$type)=='exon',])
      for(nom in names(regions)){
        o = findOverlaps(cnts,regions[[nom]],ignore.strand=TRUE)
        cnts$regionType[queryHits(o)] = nom
      }
    }else{
      warning("regionType unavailable in cnts, results will be less reliable as we cannot filter out low accuracy regions")
    }
  }
  #Filter if not still null
  if(!is.null(cnts$regionType)){
    cnts = cnts[cnts$regionType %in% regionsToUse]
  }
  #Drop Strata
  if(is.null(cnts$strata)){
    if(autoSet){
      message('strata not found, setting assuming data has been mapped to GRCh38')
      cnts$strata = NA
      o = findOverlaps(cnts,X1k@metadata$strata)
      cnts$strata[queryHits(o)] = names(X1k@metadata$strata)[subjectHits(o)]
      cnts$strata = factor(cnts$strata,levels=names(X1k@metadata$strata))
    }else{
      warning("strata unavailable in cnts, results will be less reliable as we cannot filter out regions likely to escape X-Inactivation..")
    }
  }
  #Is it still null?
  if(!is.null(cnts$strata)){
    cnts = cnts[!cnts$strata %in% dropStrata] 
  }
  cnts = cnts[as.character(seqnames(cnts))=='X']
  cnts = dropSeqlevels(cnts,setdiff(seqlevels(cnts),'X'))
  #Drop mat/pat counts if empty
  if(all(is.na(cnts$matCount)))
    cnts$matCount=NULL
  if(all(is.na(cnts$patCount)))
    cnts$patCount=NULL
  #If phased, construct allele resolved
  phased = 'altIsMum' %in% colnames(mcols(cnts))
  if(phased){
    cnts = cnts[!is.na(cnts$altIsMum),]
    cnts$matAllele = ifelse(cnts$altIsMum,cnts$altCount,cnts$refCount)
    cnts$patAllele = ifelse(cnts$altIsMum,cnts$refCount,cnts$altCount)
  }
  #Make sure regionID is defined
  if(is.null(cnts$regionID))
    cnts$regionID = as.character(cnts)
  #Drop biallelic
  if(dropBiallelic){
    #The dumb version
    #biallelicSNP = unique(cnts$regionID[cnts$refCount>0 &  cnts$altCount>0])
    #The more permissive and more sensible version.  Only difference should be that it lets through the occasional 50 vs 1 examples 
    biallelicSNP = unique(cnts$regionID[pbinom(pmin(cnts$refCount,cnts$altCount)-1,cnts$Tot,errRate,lower.tail = FALSE)<pCut])
    cnts = cnts[!cnts$regionID %in% biallelicSNP]
  }
  #Drop any SNP that's not observed often enough to be useful.  Usually this should be 2
  tmp = table(names(cnts))
  cnts = cnts[names(cnts) %in% names(tmp[tmp>=minObs])]
  return(cnts)
}


