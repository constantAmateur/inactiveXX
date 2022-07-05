#' Get valid cell barcodes for 10X experiments
#'
#' Given a set of cellranger outputs from which the X-inactivation status is to be estimated, fetches the barcodes considered to have cells.  You can either specify paths as you would when running \code{\link{Read10X}} (which is what this function uses internally), or pass the \code{BAMs} argument to be used by \code{\link{hetSNPsFromRNA}} and the path will be guessed.
#'
#' Note that \code{dataPaths} must be a named vector, as the names will be used to construct the allowable cell IDs.
#'
#' As well as getting the valid barcodes, this function will also do a basic check to ensure that the samples provided look like they only contain X chromosomes.  This is based on determining if more cells express at least one of the \code{genesX} or \code{genesY}.
#'
#'
#' @param dataPaths Named vector giving paths containing cellranger filtered cells to load.  Names must be unique giving the name to use for each 10X channel.
#' @param checkSex Should we check if all samples contain only X chromosomes?
#' @param genesX Genes specific to X chromosome, used to check sex.
#' @param genesY Genes specific to Y chromosome, used to check sex.
#' @return A vector of cell barcodes to be passed to the \code{cellsToUse} param of \code{\link{hetSNPsFromRNA}}
#' @importFrom Seurat Read10X
#' @importFrom Matrix colSums
#' @export
getCellBarcodes = function(dataPaths,checkSex=TRUE,genesX='XIST',genesY=c('RPS4Y1','RPS4Y2')){
  if(is.null(names(dataPaths)))
    stop("No names specified on vector dataPaths.  dataPaths must be a named vector")
  if(any(duplicated(names(dataPaths))))
    stop("dataPaths names are not unique.")
  #Check if bam files
  if(all(grepl('\\.bam$',dataPaths))){
    #Assume cellranger 3+ by default
    guessPaths = file.path(dirname(dataPaths),'filtered_feature_bc_matrix')
    w = which(!dir.exists(guessPaths))
    if(length(w)>0){
      guessPaths[w] = file.path(dirname(dataPaths[w]),'filtered_gene_bc_matrices')
      #Have to get the next level down in the heirachy if earlier
      for(ww in w){
        tmp = list.files(guessPaths[ww])
        if(length(tmp)!=1)
          stop("Could not auto-construct paths")
        guessPaths[ww] = file.path(guessPaths[ww],tmp)
      }
    }
  }else{
    #If not, just use what was given
    guessPaths = dataPaths
  }
  #Check everything exists
  if(!all(dir.exists(guessPaths))){
    w = which(!dir.exists(guessPaths))
    stop(sprintf("Counts not found for channels %s",paste(names(dataPaths)[w],collapse=', ')))
  }
  names(guessPaths) = names(dataPaths)
  #Now actually do what we need to do
  srat = Read10X(guessPaths)
  cellsToUse = colnames(srat)
  if(checkSex){
    #Get the vector telling us which column is from which sample
    sampLab = factor(gsub('_[ACGT]+(-[0-9]+)?$','',colnames(srat)),levels=names(guessPaths))
    cntsX = sapply(split(colSums(srat[genesX,,drop=FALSE])>0,sampLab),mean)
    cntsY = sapply(split(colSums(srat[genesY,,drop=FALSE])>0,sampLab),mean)
    noY = cntsX>cntsY
    if(any(!noY)){
      w = which(!noY)
      stop(sprintf("The following channels appear to posses a Y chromosome: %s",paste(names(guessPaths)[w],collapse=', ')))
    }
  }
  return(cellsToUse)
}

