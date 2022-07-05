#' SNPs on X chromosome
#'
#' Starting point is the 1000 genomes SNPs with AF greater than 1% on the X-chromosome.  This is then lifted over to GRCh38 and various bits of meta-data are layered on top.  Most importantly, each SNP is allocated to a region, a gene (if not inter-genic), and an evolutionary strata.  Certain stratas have a high probability of escaping X-Inactivation, making them less useful for estimating X-Inactivation state and so are best filtered out
#'
#' The following gives a complete(ish) description of how \code{X1k} was generated.  It is only completish as it starts from a file, \code{X1kFile}, which contains the 1000 genomes SNPs on X (hg19) with population allele frequency greater than 1%.
#' \itemize{
#' \item X1k = readRDS(X1kFile)
#' \item X1k38 = changeGenomeVersion(X1k,liftChain)
#' \item X1k38 = X1k38[as.character(seqnames(X1k38))=='X']
#' \item X1k38 = dropSeqlevels(X1k38,setdiff(seqlevels(X1k38),'X'))
#' \item #Definition from https://academic.oup.com/gbe/article/5/10/1863/522116, which is hg19
#' \item strata = c(0,2.78,5.04,8.43,30.62,55.78,75.53,99.98,130.82,145.73,155.72)*1e6
#' \item #These are the boundaries lifted over and the PAR definitions from a hg38 source
#' \item strata = c(0,2.781479,5.12,8.46,30.60,55.75,76.31,100.73,131.69,146.65,155.701383)*1e6
#' \item strata = GRanges('X',IRanges(strata,width=c(diff(strata),1e9)))
#' \item names(strata) = paste0('strata',seq_along(strata))
#' \item names(strata)[1] = 'PAR1'
#' \item names(strata)[length(strata)]='PAR2'
#' \item #So we definitely need to drop PAR regions.  Almost certainly want to drop starta 2 and 3 as well.  starta 4 and 5 are the final bits of the p arm and are less high quality than the q arm, but not catestrophically so.  Whether it's best to remove or keep them depends on how much bad SNPs fuck up the model (my feeling is not very much).  Then there's a background of homozygous SNPs miscalled as heterozygous.  If I'm too strict with throwing out regions likely to escape X-Inactivation I'll end up with proportionally more of these miscalled het SNPs.
#' \item #PAR = GRanges('X',IRanges(c(1,155701383),c(2781479,1e9)))
#' \item #Add this in as metadata
#' \item o = findOverlaps(X1k38,strata)
#' \item X1k38$strata[queryHits(o)] = names(strata)[subjectHits(o)]
#' \item X1k38$strata = factor(X1k38$strata,levels=names(strata))
#' \item #And annotate them
#' \item X1k38 = annotateSNPs(X1k38,gtf=gtf)
#' \item X1k = X1k38
#' \item X1k@metadata$strata = strata
#' \item X1k@metadata$txdb = NULL #Needed as this references external file
#' \item save(X1k,file='X1k.RData',compress='xz')
#' }
#'
#' @docType data
#'
#' @usage data(X1k)
#'
#' @source \href{http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/}{1000GenomesFTP}
#' @source \href{https://academic.oup.com/gbe/article/5/10/1863/522116}{Strata definitions}
#'
#' @format An object of class \code{"GRanges"}; see \code{\link[GenomicRanges]{GRanges}}.
#'
#' @keywords datasets
"X1k"
