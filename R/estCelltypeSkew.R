#' Infer groups of cells with significant Xi skew 
#'
#' Looks from regions of the kNN expression graph that show statistically significant deviations from the Xi skew specified (\code{tau})using the \code{miloR} library.  This requires both a count matrix (\code{toc}) and estimates of the Xi state for the same cells (\code{fit}).  
#'
#' The expression matrix can be input as either a sparse matrix, with rows being genes and columns being cell names, a \code{\link{Seurat}} object, or a \code{\link[miloR]{Milo}} object.  Xi states are expected to be given as an Xi fit object, the output of \code{\link{inferInactiveX}}.
#'
#' The milo significance testing framework usually uses a negative binomial model, that expects multiple biological replicates to calculate over-dispersion.  As this is not possible in the context of testing for Xi skew, a binomial model is instead used. 
#'
#' The output of this function contains objects that can be used with the standard \code{miloR} plotting functions, such as \code{\link[miloR]{plotNhoodGraphDA}} and \code{\link[miloR]{plotDAbeeswarm}}.
#'
#' @param toc Expression matrix to use.  Can be either a sparse matrix, Seurat object, or Milo object.
#' @param fit The output of \code{\link{inferInactiveX}}
#' @param tau Maternal fraction to test against.  If NULL, the overall value from \code{fit} is used.
#' @param resultsPassthrough If \code{toc} contains metadata, pass through these metadata columns to the results objects (if present).
#' @param prop What fraction of cells to randomly sample?
#' @param k How many neighbours to use in graph.
#' @param d How many PCs to use.
#' @param logFC If FALSE, uses linear difference instead of logFC (label remains logFC to not break milo).
#' @param useGraphFDR Use the graph FDR.
#' @return A list with two entries, the milo object and an annotated table of results.
#' @importFrom Seurat NormalizeData CreateSeuratObject as.SingleCellExperiment FindVariableFeatures ScaleData RunPCA
#' @importFrom miloR Milo buildGraph makeNhoods countCells calcNhoodDistance graph nhoods nhoodIndex nhoodDistances graphSpatialFDR buildNhoodGraph annotateNhoods
#' @importFrom stats binom.test
#' @export
estCelltypeSkew = function(toc,fit,tau=NULL,resultsPassthrough='seurat_clusters',prop=0.2,k=30,d=50,logFC=FALSE,useGraphFDR=TRUE){
  #Set tau if not set
  if(is.null(tau))
    tau = fit$tau
  #Is toc just a matrix
  if(inherits(toc,'Matrix') || inherits(toc,'matrix')){
    milo = NormalizeData(CreateSeuratObject(toc),verbose=FALSE)
    milo = FindVariableFeatures(milo,verbose=FALSE)
    milo = ScaleData(milo,verbose=FALSE)
    milo = RunPCA(milo,npcs=50,approx=FALSE,verbose=FALSE)
    milo = Milo(as.SingleCellExperiment(milo))
  }else if(inherits(toc,'Seurat')){
    milo = Milo(as.SingleCellExperiment(toc))
  }else if(inherits(toc,'Milo')){
    milo = toc
  }else{
    stop('toc is not in expected input format')
  }
  #Add in the entries for the Xi state.  Do checks for consistency between cell naming conventions
  m = match(colnames(milo),fit$cellSummary$cell)
  #If none of them match, fail
  if(all(is.na(m)))
    stop("Cell names in 'toc' and 'fit' don't match.")
  #If fewer than 20% match, warning
  if(mean(is.na(m))<0.2)
    warning(sprintf("Fit information found for only %.02f%% of cells in 'toc'.",100*mean(is.na(m))))
  #Don't use if not very certain (posterior<90%) or no information 
  milo@colData$stateXi = fit$cellSummary$stateXi[m]
  #Drop those that are not determined 
  w = which(milo@colData$stateXi!='Undetermined')
  milo = milo[,w]
  #Good to go, build graphs and proceed.
  milo = buildGraph(milo,k=k,d=d)
  milo = makeNhoods(milo,prop=prop,k=k,d=d,refined=TRUE)
  tmp = quantile(colSums(nhoods(milo)))
  message(sprintf('Neighbourhoods have median size %d (range %d to %d)',round(tmp[3]),tmp[1],tmp[5]))
  if(tmp[3]<20)
    stop("Neighbourhoods to small, increase prop or k or both")
  #Now generate counts in neighbourhoods
  milo = countCells(milo, meta.data = data.frame(milo@colData), samples="stateXi")
  milo = calcNhoodDistance(milo, d=d)
  #Do our own test because of reasons
  dd = data.frame(milo@nhoodCounts)
  dd$tot = rowSums(dd)
  dd$frac = dd$Maternal/dd$tot
  dd$pVal = sapply(seq(nrow(dd)),function(e) binom.test(dd$Maternal[e],dd$tot[e],tau)$p.value)
  #This needs to mimic the results format exactly
  res = data.frame(logFC = if(logFC){log(dd$frac)-log(tau)}else{dd$frac-tau},
                   logCPM = 14, #Dummy values
                   `F` = 0.2, #Dummy values
                   PValue = dd$pVal,
                   FDR = p.adjust(dd$pVal,method='BH'),
                   Nhood = as.numeric(rownames(dd)))
  if(!useGraphFDR){
    res$SpatialFDR = res$FDR
  }else{
    tmp = graphSpatialFDR(x.nhoods = nhoods(milo),
                          graph = graph(milo),
                          k=milo@.k,
                          pvalues = res$PValue[order(res$Nhood)],
                          indices = nhoodIndex(milo),
                          distances = nhoodDistances(milo),
                          reduced.dimensions = milo@int_colData$reducedDim$PCA)
    res$SpatialFDR[order(res$Nhood)] = tmp
  }
  #Make thing for visualising
  milo = buildNhoodGraph(milo)
  #Add in metadata
  for(tgt in intersect(resultsPassthrough,colnames(milo@colData)))
    res = annotateNhoods(milo, res, coldata_col = tgt)
  #Now return the two objects for you to with what you will
  return(list(milo=milo,res=res))
}


