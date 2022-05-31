#' Compute eigengenes following gamma clustering
#'
#' @param x A data matrix
#' @param m Membership probabilities
#' @return A matrix of eigengenes where rows correspond to each cluster.
#' @examples
#' eigenGenes(foldChange, clusters$membership)
eigenGenes <- function(x, m) {
  xm <- apply(x, 1, mean)
  xsd <- apply(x, 1, sd)
  xt <- (x - xm) / xsd
  egenes <- t(apply(m, 1, function(mi) {
    .getEigen(xt, mi)
  }))
}

.getEigen <- function(xt, m) {
  if (sum(m) == 0) {
    egene <- rep(0, dim(xt)[2])
  } else {
    egene <- prcomp(xt * m)$rotation[, 1]
    weightMean <- apply(xt, 2, function(b) {
      weighted.mean(b, m)
    })
    if (cor(egene, weightMean) < 0) {
      egene <- -1 * egene
    }
  }
  return(egene)
}

#' Computes a p-value matrix from the difference of two p-value matrices.
#'
#' @param pmat1 First p-value matrix
#' @param pmat2 Second p-value matrix
#' @param twoTailed Compute a two-tailed p-value
#' @param unifOut Force the output to a uniform distribution
#' @return A new p-value matrix testing the ratio of the first p-value matrix
#' over the second (or vice versa if twoTailed = TRUE).
#' @examples
#' pm1 <- pMatrix(foldChange)
#' pm2 <- pMatrix(foldChange, z=traitData["TCDD",])
#' lapMat <- dLaplace(pm1, pm2)
dLaplace <- function(pmat1,
                     pmat2,
                     twoTailed = FALSE,
                     unifOut = TRUE) {
  pmatDiff <- log(pmat2) - log(pmat1)
  diag(pmatDiff) <- 0
  p <- 0.5 * exp(-abs(pmatDiff))
  if (twoTailed) {
    pval <- 2 * p
  } else {
    pval <- ifelse(pmatDiff > 0, p, 1 - p)
  }
  N <- dim(pmat2)[1]
  if (unifOut) {
    pval <- forceUniform(pval, lower.tail = TRUE)
  }
  return(pval)
}
