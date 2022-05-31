#' Visualize clusters with UMAP
#'
#' @param pmat A p-value matrix
#' @param labels Cluster labels to visualize
#' @examples
#' pm <- pMatrix(foldChange)
#' clusters <- gammaCluster(pm)
#' gammaPlot(pm, clusters$labels)
gammaPlot <- function(pmat, labels) {
  colorList <- scales::hue_pal()(max(labels))
  u <- as.data.frame(umap::umap(-log1p(-pmat), input = "dist")$layout)
  colnames(u) <- c("UMAP1", "UMAP2")
  u$clusters <- as.factor(labels)
  ggplot2::ggplot() + ggplot2::coord_cartesian() +
    ggplot2::scale_x_continuous() +
    ggplot2::scale_y_continuous() +
    ggplot2::scale_color_hue() +
    ggplot2::layer(
      data = u,
      mapping = ggplot2::aes(x = UMAP1,
                             y = UMAP2,
                             color = clusters),
      stat = "identity",
      geom = "point",
      position = ggplot2::position_jitter()
    ) +
    ggplot2::theme_classic()
}

#' Perform gamma clustering to identify clusters that maximize the probability
#' of null inter-cluster correlations
#'
#' @param pmat A p-value matrix
#' @param div.maxiter Maximum number of iterations for the divisive phase
#' @param op.maxiter Maximum number of iterations for the optimization phase
#' @return A list containing cluster labels and membership probabilities
#' @examples
#' pm <- pMatrix(foldChange)
#' clusters <- gammaCluster(pm)
gammaCluster <- function(pmat, div.maxiter=100, op.maxiter=200){
  A <- -log(pmat)
  diag(A)<- 0
  D <- -log1p(-pmat)
  diag(D) <- 0

  m <- .divGamma(A, D, div.maxiter)
  m <- .opGamma(m, A, op.maxiter)
  m <- m[unique(Rfast::colMaxs(m)),]
  m <- t(t(m) / colSums(m))
  clust <- Rfast::colMaxs(m)
  return(list(labels=clust, membership=m))
}

.opGamma <- function(m, A, maxiter){
  clusters <- Rfast::colMaxs(m, value = FALSE)
  for(i in seq_len(maxiter)){
    a <- m %*% A
    total <- Rfast::rowsums(m) - m
    p <- -pgamma(a, total, lower.tail = FALSE, log.p = TRUE)
    maxCol <- Rfast::colMaxs(p, parallel = TRUE)
    upper <- t(exp(t(p) - maxCol))
    total <- Rfast::colsums(upper, parallel = TRUE)
    m <- t(t(upper) / total)
    clusters0 <- Rfast::colMaxs(m, value = FALSE)
    if (Rfast::all_equals(clusters0, clusters, fast_result = TRUE)) {
      converged <- 1
      break
    }
    clusters <- clusters0
  }
  return(m)
}

.divGamma <- function(A, D, maxiter){
  m <- matrix(rep(1, dim(A)[1]), nrow=1)
  nmiss <- Inf
  c0 <- c()

  for (i in seq_len(maxiter)){
    a <- m %*% A
    d <- m %*% D

    total <- Rfast::rowsums(m) - m
    p1 <- -pgamma(a, total, lower.tail=FALSE, log.p = TRUE)
    p2 <- pgamma(d, total, lower.tail=FALSE, log.p = TRUE)
    p <- rbind(p1+p2, rep(0, dim(A)[1]))
    maxCol <- Rfast::colMaxs(p, parallel = TRUE)
    upper <- t(exp(t(p) - maxCol))
    total <- Rfast::colsums(upper, parallel = TRUE)
    m <- t(t(upper) / total)
    nmiss0 <- length(which(Rfast::colMaxs(m)==dim(m)[1]))
    if(nmiss0==nmiss){
      m <- m[unique(Rfast::colMaxs(m)),]
      break
    }
    nmiss <- nmiss0
  }
  return(m)
}

#' Gamma sub-clustering performed on each individual cluster before a final
#' optimization step. The output can be re-entered as input for recursive
#' optimization.
#'
#' @param pmat A p-value matrix
#' @param labels A vector of previous cluster labels
#' @param div.maxiter Maximum number of initial splits of each cluster
#' @param op.maxiter Maximum number of iterations for final optimization
#' @return A list of cluster labels and membership probabilities
#' @examples
#' pm <- pMatrix(foldChange)
#' clusters <- gammaCluster(pm)
#' rcluster <- recGamma(pmat, clusters$labels)
#' rcluster <- recGamma(pmat, rcluster$labels)
recGamma <- function(pmat, labels, div.maxiter=10, op.maxiter=100){
  N <- length(labels)
  k <- max(labels)
  sclust <- rep(0, N)
  for (i in seq_len(k)){
    if (sum(labels==i)<=2){
      sclust[labels==i]<-rep(max(sclust)+1, sum(labels==i))
    }else{
      pmat_sub<- forceUniform(pmat[labels==i, labels==i],lower.tail=TRUE)
      sublabels <- gammaCluster(pmat_sub, div.maxiter, op.maxiter)$labels
      sclust[labels==i]<- sublabels+1+max(sclust)
    }
  }
  sclust <- data.table::frank(sclust, ties.method="dense")
  k <- max(sclust)
  m <- matrix(rep(0, N*k), nrow=k)
  m[cbind(sclust,seq_len(N))]<-1
  A <- -log(pmat)
  diag(A)<- 0
  m <- .opGamma(m, A, op.maxiter)
  sclust <- Rfast::colMaxs(m)
  m <- m[sort(unique(sclust)),]
  m <- t(t(m) / colSums(m))
  sclust <- data.table::frank(sclust, ties.method = "dense")
  return(list(labels=sclust, membership=m))
}
