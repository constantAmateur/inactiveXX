% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estCelltypeSkew.R
\name{estCelltypeSkew}
\alias{estCelltypeSkew}
\title{Infer groups of cells with significant Xi skew}
\usage{
estCelltypeSkew(
  toc,
  fit,
  tau = NULL,
  resultsPassthrough = "seurat_clusters",
  prop = 0.2,
  k = 30,
  d = 50,
  logFC = FALSE,
  useGraphFDR = TRUE
)
}
\arguments{
\item{toc}{Expression matrix to use.  Can be either a sparse matrix, Seurat object, or Milo object.}

\item{fit}{The output of \code{\link{inferInactiveX}}}

\item{tau}{Maternal fraction to test against.  If NULL, the overall value from \code{fit} is used.}

\item{resultsPassthrough}{If \code{toc} contains metadata, pass through these metadata columns to the results objects (if present).}

\item{prop}{What fraction of cells to randomly sample?}

\item{k}{How many neighbours to use in graph.}

\item{d}{How many PCs to use.}

\item{logFC}{If FALSE, uses linear difference instead of logFC (label remains logFC to not break milo).}

\item{useGraphFDR}{Use the graph FDR.}
}
\value{
A list with two entries, the milo object and an annotated table of results.
}
\description{
Looks from regions of the kNN expression graph that show statistically significant deviations from the Xi skew specified (\code{tau})using the \code{miloR} library.  This requires both a count matrix (\code{toc}) and estimates of the Xi state for the same cells (\code{fit}).
}
\details{
The expression matrix can be input as either a sparse matrix, with rows being genes and columns being cell names, a \code{\link{Seurat}} object, or a \code{\link[miloR]{Milo}} object.  Xi states are expected to be given as an Xi fit object, the output of \code{\link{inferInactiveX}}.

The milo significance testing framework usually uses a negative binomial model, that expects multiple biological replicates to calculate over-dispersion.  As this is not possible in the context of testing for Xi skew, a binomial model is instead used. 

The output of this function contains objects that can be used with the standard \code{miloR} plotting functions, such as \code{\link[miloR]{plotNhoodGraphDA}} and \code{\link[miloR]{plotDAbeeswarm}}.
}
