% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/matrixCalcs.R
\name{forceUniform}
\alias{forceUniform}
\title{Force an affinity or dissimilarity measure to a uniform distribution for use
with gamma clustering}
\usage{
forceUniform(x, lower.tail = FALSE)
}
\arguments{
\item{x}{A square affinity or dissimilarity matrix}

\item{lower.tail}{Set "TRUE" for dissimilarity and "FALSE" for affinity}
}
\value{
A p-value matrix that has been forced to a uniform distribution
}
\description{
Force an affinity or dissimilarity measure to a uniform distribution for use
with gamma clustering
}
\examples{
db <- distBoot(foldChange, partial=TRUE)
}
