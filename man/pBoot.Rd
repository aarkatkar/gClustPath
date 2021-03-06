% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/minPath.R
\name{pBoot}
\alias{pBoot}
\title{Test a given path length connecting two traits based on a bootstrapped
path length distribution of path lengths between genes}
\usage{
pBoot(d, z)
}
\arguments{
\item{d}{A vector of path lengths obtained from distBoot}

\item{z}{The path length being tested}
}
\value{
The p-value testing the shortness of the provided path-length
}
\description{
Test a given path length connecting two traits based on a bootstrapped
path length distribution of path lengths between genes
}
\examples{
db <- distBoot(foldChange, partial=TRUE)
p <- pBoot(db, 0.016)
}
