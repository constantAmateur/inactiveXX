#' Infer X inactivation state for all cells
#'
#' This is a wrapper script that will perform the standard workflow to estimate the X-Inactivation state for every cell in a 10X experiment from the BAM files.
#'
#' @param bams10X BAM files output by cellranger.  Should include every cellranger BAM file relating to one individual with at least two X chromosomes.
#' @param refGenome The reference genome used to map the BAMs provided.
#' @param verbose Be verbose?
#' @return The output of \code{\link{inferInactiveX}}.
#' @export
inactiveXX = function(bams10X,refGenome,verbose=TRUE){
  if(verbose)
    message(sprintf("Finding heterozygous SNPs from %d BAM files",length(bams10X)))
  snps = hetSNPsFromRNA(bams10X,refGenome)
  #There's no need to go back to the BAMs a second time, get XCnts from same object
  XCnts = filterCountsX(snps)
  fit = inferInactiveX(XCnts)
  return(fit)
}
