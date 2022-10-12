#' Produces summary of cell groups from fit objects
#'
#' Produces table summarising the Xi skew for each cell type in \code{cellsToUse} for each individual in \code{fits}.  Individiual and cell type names are taken from the names of \code{cellsToUse} and \code{fits} if present.
#' 
#' @param fits A list containing fit objects produced by \code{\link{inferInactiveX}}, giving Xi states from multiple individuals.
#' @param cellsToUse If NULL, all cells.  If a vector of cell IDs, pop size for that cell type.  If a list of cellIDs, pop size for all.
#' @param highConfOnly Use all estimates of Xi, or just the high confidence calls?
#' @return A data.table giving the Xi summary stats by cell group and individual from fit objects.
cellGroupSummary = function(fits,cellsToUse=NULL,highConfOnly=FALSE){
  if(is.null(cellsToUse))
    message("No cells given, using all cells to estimate population bottleneck size")
  if(is.null(names(fits)))
    names(fits) = sprintf('Individual%d',seq_along(fits))
  #If we have multiple cell types, run each separately!
  if(!is.null(cellsToUse) && inherits(cellsToUse,'list')){
    out = list()
    if(is.null(names(cellsToUse)))
      names(cellsToUse) = sprintf('CellGroup%d',seq_along(cellsToUse))
    for(i in seq_along(cellsToUse)){
      out[[i]] = cellGroupSummary(fits,cellsToUse=cellsToUse[[i]],highConfOnly=highConfOnly)
      out[[i]]$cellGroup = names(cellsToUse)[i]
    }
    #Merge table
    out = do.call(rbind,out)
    return(out)
  }
  #Get just the cells we want from each individual
  dat = list()
  for(i in seq_along(fits)){
    if(is.null(fits[[i]]$cellSummary))
      stop(sprintf("Fit object %d, named %s, is not a valid Xi skew fit object.",i,names(fits)[i]))
    dat[[i]] = fits[[i]]$cellSummary
    if(!is.null(cellsToUse))
      dat[[i]] = dat[[i]][dat[[i]]$cell %in% cellsToUse,,drop=FALSE]
    if(highConfOnly)
      dat[[i]] = dat[[i]][dat[[i]]$highConfCall,,drop=FALSE]
  }
  names(dat) = names(fits)
  if(any(sapply(dat,nrow)==0)){
    message(sprintf("%d of %d individuals have no matching cells and will not be used.",sum(sapply(dat,nrow)==0),length(dat)))
    if(all(sapply(dat,nrow)==0))
      stop("No matching cells found.  Check cellsToUse correctly specified.")
    dat = dat[sapply(dat,nrow)>0]
  }
  #Create a summary for each sample
  dd = data.frame(individual = names(dat),
                  cellGroup = NA,
                  nMat = sapply(dat,function(e) sum(e$stateProbs)),
                  nPat = sapply(dat,function(e) sum(1-e$stateProbs)),
                  nTot = sapply(dat,nrow))
  dd$matFrac = dd$nMat/dd$nTot
  rownames(dd)=NULL
  return(dd)
}
