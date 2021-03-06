% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/minPath.R
\name{gMST}
\alias{gMST}
\title{Identify local clusters of correlated genes using a minimum spanning tree.}
\usage{
gMST(pmat, cutoff = 0.05, min.size = 3)
}
\arguments{
\item{pmat}{A p-value matrix}

\item{cutoff}{The p-value below which an edge is removed}

\item{min.size}{The smallest cluster size to report in the returned vector}
}
\value{
A vector of cluster labels
}
\description{
Identify local clusters of correlated genes using a minimum spanning tree.
}
\examples{
pm <- pMatrix(foldChange, unsigned=TRUE, output="orig")
clusters <- gMST(pMatrix, cutoff=0.05/dim(foldChange)[1])
}
