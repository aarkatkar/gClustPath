#' Compute a p-value matrix
#'
#' @param features A data matrix of features. Rows are features.
#' @param traits A data matrix of traits. Rows are traits.
#' @param z A matrix or vector of covariates for conditioning correlations
#' @param unsigned If TRUE, ignore the sign of correlations
#' @param partial If TRUE, compute a partial correlation matrix
#' @param output One of "unif", "orig", or "cor" corresponding to forced-uniform,
#' original p-values, or raw correlations, respectively
#' @return A (N x N) matrix of p-values. The first rows/columns correspond to
#' traits if supplied
#' @examples
#' pm <- pMatrix(foldChange)
pMatrix <- function(features,
                    traits = NULL,
                    z = NULL,
                    unsigned=FALSE,
                    partial=FALSE,
                    output="unif"){
  if (output != "unif" & output != "orig" & output != "cor") {
    warning('output must be "unif", "orig", or "cor"')
    return(NA)
  }
  correlations <- .getCorrel(features, traits, z, partial)
  corMat <- correlations$corMat
  if (output == "cor") {
    return(corMat)
  }

  n_phen <- correlations$n_phen
  n_cov <- correlations$n_cov
  # Get p-values
  if (partial | output == "unif") {
    p <- .unifPvalue(corMat, n_phen, dim(features)[2], unsigned)
  } else {
    p <- .tdistPvalue(corMat, n_phen, n_cov, dim(features)[2], unsigned)
  }
  return(p)
}

.unifPvalue <- function(corMat, n_phen, n_samp, squared) {
  N <- dim(corMat)[1]
  if (n_phen == 0) {
    rxy <- corMat
  } else {
    rxy <- corMat[(n_phen + 1):N, (n_phen + 1):N]
    dm <- corMat[1:n_phen, 1:N]
    suppressWarnings(tapp <- dm * sqrt((n_samp - 2) / (1 - dm^2)))
    if (squared) {
      papp <- 2 * pt(-abs(tapp), n_samp - 2)
    } else {
      papp <- pt(-tapp, n_samp - 2)
    }
  }
  if (squared) {
    rxy <- rxy^2
  }
  p <- forceUniform(rxy)
  if (n_phen != 0) {
    p <- rbind(papp[1:n_phen, (n_phen + 1):N], p)
    p <- cbind(t(papp), p)
  }
  return(p)
}

.tdistPvalue <- function(corMat, n_phen, n_cov, n_samp, squared) {
  N <- dim(corMat)[1]
  if (n_phen == 0) {
    rxy <- corMat
  } else {
    rxy <- corMat[(n_phen + 1):N, (n_phen + 1):N]
    dm <- corMat[1:n_phen, 1:N]
    suppressWarnings(tapp <- dm * sqrt((n_samp - 2) / (1 - dm^2)))
    if (squared) {
      papp <- 2 * pt(-abs(tapp), n_samp - 2)
    } else {
      papp <- pt(-tapp, n_samp - 2)
    }
  }
  rxy <- ifelse(rxy >= 1, 1, rxy)
  suppressWarnings(t <-
                     rxy * sqrt((n_samp - 2 - n_cov) / (1 - rxy^2)))
  if (squared) {
    p <- 2 * pt(-abs(t), n_samp - 2 - n_cov)
  } else {
    p <- pt(-t, n_samp - 2 - n_cov)
  }
  if (n_phen != 0) {
    p <- rbind(papp[1:n_phen, (n_phen + 1):N], p)
    p <- cbind(t(papp), p)
  }
  p <- matrix(c(p), nrow = N)
  return(p)
}


.getCorrel <- function(x, phen, z, partial) {
  if (partial) {
    corMat <- corpcor::pcor.shrink(t(x), verbose = FALSE)
    n_cov <- 0
  } else if (is.null(z)) {
    corMat <- coop::pcor(t(x), use = "complete.obs")
    n_cov <- 0
  } else {
    if (is.null(dim(z))) {
      z <- t(z)
    }
    corMat <- .covCor(x, z)
    n_cov <- dim(z)[1]
  }
  if (!is.null(phen)) {
    if (is.null(dim(phen))) {
      phen <- t(phen)
    }
    corMat <- .addPhen(corMat, x, phen)
    n_phen <- dim(phen)[1]
  } else {
    n_phen <- 0
  }
  return(list(
    corMat = corMat,
    n_phen = n_phen,
    n_cov = n_cov
  ))
}

.covCor <- function(x, z) {
  N <- dim(x)[1]
  xtemp <- rbind(x, z)
  rxy <- coop::pcor(t(xtemp), use = "complete.obs")
  for (i in 1:dim(z)[1]) {
    rxz <- rxy[N + i, ]
    suppressWarnings(sqrxz <- ifelse(rxz > 1, 0, sqrt(1 - rxz^2)))
    denom <- sqrxz %*% t(sqrxz)
    rxy <- ifelse(denom == 0, 1, (rxy - rxz %*% t(rxz)) / denom)
  }
  rxy <- rxy[1:N, 1:N]
  return(rxy)
}

.addPhen <- function(rxy, x, phen) {
  N <- dim(x)[1]
  nphen <- dim(phen)[1]
  dm <- cor(t(phen), t(rbind(phen, x)), use = "complete.obs")
  rxy <- rbind(dm[1:nphen, (nphen + 1):(N + nphen)], rxy)
  rxy <- cbind(t(dm), rxy)
  return(rxy)
}

#' Force an affinity or dissimilarity measure to a uniform distribution for use
#' with gamma clustering
#'
#' @param x A square affinity or dissimilarity matrix
#' @param lower.tail Set "TRUE" for dissimilarity and "FALSE" for affinity
#' @return A p-value matrix that has been forced to a uniform distribution
#' @examples
#' db <- distBoot(foldChange, partial=TRUE)
forceUniform <- function(x, lower.tail = FALSE) {
  vals <- x[upper.tri(x)]
  N_ <- length(vals)
  if (lower.tail) {
    rankMat <- data.table::frank(vals) / (N_ + 1)
  } else {
    rankMat <- (N_ + 1 - data.table::frank(vals)) / (N_ + 1)
  }
  p <- diag(0, nrow = dim(x)[1])
  p[upper.tri(p)] <- rankMat
  return(p + t(p))
}
