#' Plot model heatmap
#'
#' Make a plot summarising the model.  Specifically plots a heatmap of SNPs by cells.  This is annotated by various summary metrics for each SNP and cell, which are calculated assuming the other data modality is true (i.e., number of bad reads at each SNP calculated assuming the cell assignment is correct).
#'
#'  By default, the heatmap shows the (usually very sparse) raw data.  However, if \code{summariseBy} is set to \code{"SNP"} or \code{"cell"}, the heatmap is distorted to remove the sparsity and emphasise the performance of either the cell or SNP (depending on the value of \code{summariseBy}).
#'
#' When summarising by SNP (or cell), all NA values are removed from each row (or column).  This allows visualisation of the overall data for each SNP (or cell), at the cost of making comparisons between rows (or columns) more difficult.
#'
#' @param fit The output of \code{\link{inferInactiveX}}.
#' @param summariseBy Should the heatmap summarise by SNP, cell, or do no manipulation of data (the default).
#' @param tol Binarise cells into maternal/paternal with this tolerance.
#' @param minCov Only plot those SNPs with at least this many reads.
#' @param badCutoff Any SNP where the 95\% confidence interval for counts going against the model exceeds this fraction is marked as "bad".
#' @param truncAt Truncate logit probabilities with absolute value above this.
#' @return Nothing.  A heatmap is drawn.
#' @importFrom ComplexHeatmap Heatmap rowAnnotation columnAnnotation draw
#' @importFrom circlize colorRamp2
#' @importFrom stats quantile
#' @export
plotModelHeatmap = function(fit,summariseBy=c('none','SNP','cell'),tol=0.1,minCov=0,badCutoff=0.10,truncAt=10){
  summariseBy = match.arg(summariseBy)
  #########################
  # Construct ordered data
  dd = fit$dd
  mat = matrix(NA,nrow=length(unique(dd$SNP)),ncol=length(unique(dd$cell)),dimnames=list(unique(dd$SNP),unique(dd$cell)))
  mat[cbind(match(dd$SNP,unique(dd$SNP)),match(dd$cell,unique(dd$cell)))] = dd$refCount/dd$tot
  #Define how to split and order things
  x = fit$cellSummary[colnames(mat),]
  ssCol = factor(ifelse(x$highConfCall,ifelse(x$matFrac>0.5,'Pat','Mat'),'???'),levels=c('Mat','???','Pat'))
  ooCol = colnames(mat)[order(ssCol,-fit$states[colnames(mat)])]
  ssCol = ssCol[match(ooCol,colnames(mat))]
  mat = mat[,ooCol]
  #And the same for rows
  ssRow = factor(ifelse(fit$genotype[rownames(mat)]==1,'Mat','Pat'),levels=c('Mat','Pat'))
  ooRow = rownames(mat)[order(ssRow,-fit$snpSummary[rownames(mat),'tot'])]
  ssRow = ssRow[match(ooRow,rownames(mat))]
  mat = mat[ooRow,]
  #########
  # Filter
  toKeep = fit$snpSummary[rownames(mat),'tot']>minCov
  mat = mat[toKeep,,drop=FALSE]
  ssRow = ssRow[toKeep]
  if(nrow(mat)<5){
    warning(sprintf("Fewer than 5 SNPs with coverage over %d.  Not plotting.",minCov))
    return(NULL)
  }
  #############
  # Annotation
  x = fit$cellSummary[colnames(mat),]
  colAnnot = columnAnnotation(stateProb = pmin(truncAt,pmax(-1*truncAt,x$states)),
                              matFrac = x$matFrac,
                              badFrac = x$badFrac,
                              covCell = x$tot,
                              col = list(stateProb=colorRamp2(c(-truncAt,0,truncAt),c('#f1a340','#f7f7f7','#998ec3')),
                                         matFrac = colorRamp2(c(0,0.5,1),c('#d8b365','#f5f5f5','#5ab4ac')),
                                         badFrac = colorRamp2(c(0,.25,.5),c('#fee6ce','#fdae6b','#e6550d')),
                                         covCell = colorRamp2(c(1,quantile(x$tot,c(0.99))),c('#f7f7f7','#525252'))
                                         ),
                           annotation_name_side='left')
  y = fit$snpSummary[rownames(mat),]
  rowAnnot = rowAnnotation(snpIsBad = y$highConfCall,
                           matFrac = y$matFrac,
                           badFrac = y$badFrac,
                           covSNP = y$totOff,
                           col = list(snpIsBad = c(`TRUE`='darkred',`FALSE`='darkblue'),
                                      matFrac = colorRamp2(c(0,0.5,1),c('#d8b365','#f5f5f5','#5ab4ac')),
                                      badFrac = colorRamp2(c(0,.25,.5),c('#fee6ce','#fdae6b','#e6550d')),
                                      covSNP = colorRamp2(c(1,quantile(x$tot,c(0.99))),c('#f7f7f7','#525252'))
                                      ))
  ######
  # Raw
  if(summariseBy=='none'){
    ht = Heatmap(mat,
            name='refFrac',
            col = colorRamp2(c(0,0.5,1),c('#b2182b','#f7f7f7','#2166ac')),
            cluster_rows=FALSE,
            cluster_columns=FALSE,
            show_column_names = FALSE,
            show_row_names=FALSE,
            column_split = ssCol,
            row_split = ssRow,
            top_annotation = colAnnot,
            right_annotation =rowAnnot 
            )
    draw(ht, row_title = "Ref allele is",column_title='Cell assigned to')
  }else{
    #Here's my crazy alternative.  mark boundaries between SNPs with an NA.  Fill each SNP out to show individual reads.  Normalise length to be constant.
    refCnt = altCnt = mat
    refCnt[cbind(match(dd$SNP,rownames(refCnt)),match(dd$cell,colnames(refCnt)))] = dd$refCount
    altCnt[cbind(match(dd$SNP,rownames(altCnt)),match(dd$cell,colnames(altCnt)))] = dd$altCount
    snpIdxs = split(colnames(refCnt),ssCol)
    cellIdxs =split(rownames(refCnt),ssRow)
    if(summariseBy=='SNP'){
      out = list()
      for(i in seq(nrow(refCnt))){
        out[[i]] =list()
        for(nom in names(snpIdxs)){
          snpIdx = snpIdxs[[nom]]
          r = refCnt[i,snpIdx]
          a = altCnt[i,snpIdx]
          w = which(!is.na(r) | !is.na(a))
          r=r[w]
          a=a[w]
          if(length(r)==0){
            out[[i]][[nom]] = NA
          }else{
            #Now intersperse, with blockers
            mask = rep(NA,length(r)*3)
            mask[seq(1,length(mask),3)]=r
            mask[seq(2,length(mask),3)]=a
            mask[seq(3,length(mask),3)]=0
            mask = rep(rep(c(1,0,NA),length(r)),mask)
            out[[i]][[nom]] = mask
          }
        }
      }
      #Now need to normalise each block
      sizes = do.call(rbind,lapply(out,lengths))
      tgtSizes = apply(sizes,2,max)
      snpMat = list()
      for(nom in names(snpIdxs)){
        x = lapply(out,function(e) e[[nom]])
        tgtSize = tgtSizes[nom]
        idxs = seq(tgtSize)
        for(i in seq_along(x)){
          f = tgtSize/length(x[[i]])
          x[[i]] = x[[i]][floor((idxs-1)/f)+1]
        }
        snpMat[[nom]] = do.call(rbind,x)
      }
      snpMat = do.call(cbind,snpMat)
      ht = Heatmap(snpMat,
                   name='refFrac',
                   col = colorRamp2(c(0,0.5,1),c('#b2182b','#f7f7f7','#2166ac')),
                   row_split = ssRow,
                   right_annotation = rowAnnot,
                   column_split = factor(rep(names(tgtSizes),tgtSizes),levels=levels(ssCol)),
                   show_row_names=FALSE,
                   show_column_names=FALSE,
                   cluster_rows=FALSE,
                   cluster_columns=FALSE)
      draw(ht, row_title = "Ref allele is",column_title='Cell assigned to')
    }
    if(summariseBy=='cell'){
      #Do the reverse version.
      out = list()
      for(i in seq(ncol(refCnt))){
        out[[i]] =list()
        for(nom in names(cellIdxs)){
          cellIdx = cellIdxs[[nom]]
          r = refCnt[cellIdx,i]
          a = altCnt[cellIdx,i]
          w = which(!is.na(r) | !is.na(a))
          r=r[w]
          a=a[w]
          if(length(r)==0){
            out[[i]][[nom]] = NA
          }else{
            #Now intersperse, with blockers
            mask = rep(NA,length(r)*3)
            mask[seq(1,length(mask),3)]=r
            mask[seq(2,length(mask),3)]=a
            mask[seq(3,length(mask),3)]=0
            mask = rep(rep(c(1,0,NA),length(r)),mask)
            out[[i]][[nom]] = mask
          }
        }
      }
      #Now need to normalise each block
      sizes = do.call(rbind,lapply(out,lengths))
      tgtSizes = apply(sizes,2,max)
      cellMat = list()
      for(nom in names(cellIdxs)){
        x = lapply(out,function(e) e[[nom]])
        tgtSize = tgtSizes[nom]
        idxs = seq(tgtSize)
        for(i in seq_along(x)){
          f = tgtSize/length(x[[i]])
          x[[i]] = x[[i]][floor((idxs-1)/f)+1]
        }
        cellMat[[nom]] = do.call(cbind,x)
      }
      cellMat = do.call(rbind,cellMat)
      ht = Heatmap(cellMat,
                   name='refFrac',
                   col = colorRamp2(c(0,0.5,1),c('#b2182b','#f7f7f7','#2166ac')),
                   row_split = factor(rep(names(tgtSizes),tgtSizes),levels=levels(ssCol)),
                   column_split = ssCol,
                   top_annotation = colAnnot,
                   show_row_names=FALSE,
                   show_column_names=FALSE,
                   cluster_rows=FALSE,
                   cluster_columns=FALSE)
      draw(ht, row_title = "Ref allele is",column_title='Cell assigned to')
    }
  }
}
