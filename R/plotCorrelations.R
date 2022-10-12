#' Plot Xi skew correlations
#'
#' Visualise the output of \code{\link{estCorrelation}} as a heatmap.  The organisation of the heatmap is driven by the point estimates of the correlation coefficient for each pair of cell groupings, while each tile shows the distribution over the uncertainty in those estimates.
#'
#' The output is a \code{ComplexHeatmap} object and so can be controlled in the normal way for that package.
#'
#' If you encounter errors about NAs, try disabling hierarchical clustering by passing \code{cluster_rows=FALSE} and \code{cluster_columns=FALSE} to \code{\link[ComplexHeatmap]{Heatmap}}.
#'
#' @param rhos Output of \code{\link{estCorrelation}}.
#' @param nBins How many bins to use for histograms showing uncertainty distribution. 
#' @param ... Passed to \code{\link{Heatmap}}.
#' @return A heatmap object.
#' @importFrom grid grid.rect grid.lines gpar
#' @importFrom graphics hist
#' @importFrom ComplexHeatmap Heatmap
#' @importFrom circlize colorRamp2
#' @export
plotCorrelations = function(rhos,nBins=50,...){
  xBig = attributes(rhos)$rawSamples
  #Make the point estimates into a martix
  xAvg = xBig[[1]]
  xAvg[cbind(match(rhos$cellGroupA,rownames(xAvg)),match(rhos$cellGroupB,colnames(xAvg)))] = rhos$pointEst
  xAvg[cbind(match(rhos$cellGroupB,rownames(xAvg)),match(rhos$cellGroupA,colnames(xAvg)))] = rhos$pointEst
  #Define the magic function that will draw distributions in cell
  cell_fun = function(j,i,x,y,width,height,fill){
    #White out the area
    #grid.rect(x = x, y = y, width = width, height = height,gp = gpar(col = NA, fill = 'white'))
    #Grab the raw data
    dat = unlist(lapply(xBig,function(e) e[i,j]))
    #Make histogram
    brks = seq(-1,1,length.out=nBins)
    t1 = hist(dat,breaks=brks,plot=FALSE)
    heights = t1$counts/max(t1$counts)
    #Draw histogram
    for(ii in seq_along(t1$counts)){
      grid.rect(x=x - width/2 + width*((ii-1)/nBins),
                y=y-height/2,
                just=c(0,0),
                width= width/nBins,
                height=heights[ii]*height,
                gp = gpar(col=NA,fill='black'))
    }
    #Vertical line in middle
    grid.lines(x=c(x,x),y=c(y-height/2,y+height/2),gp=gpar(lty=2))
  }
  hmParams = list(matrix=xAvg,
                  name='Correlation',
                  cell_fun=cell_fun)
  theDots = list(...)
  for(nom in names(theDots))
    hmParams[[nom]] = theDots[[nom]]
  hm = do.call(Heatmap,hmParams)
  return(hm)
}


