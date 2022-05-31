#' Identify local clusters of correlated genes using a minimum spanning tree.
#'
#' @param pmat A p-value matrix
#' @param cutoff The p-value below which an edge is removed
#' @param min.size The smallest cluster size to report in the returned vector
#' @return A vector of cluster labels
#' @examples
#' pm <- pMatrix(foldChange, unsigned=TRUE, output="orig")
#' clusters <- gMST(pMatrix, cutoff=0.05/dim(foldChange)[1])
gMST <- function(pmat, cutoff = 0.05, min.size = 3) {
  g <- igraph::graph.adjacency(-log1p(-pmat), weighted = TRUE)
  mst <- igraph::mst(g)
  drop <- which(igraph::E(mst)$weight > cutoff)
  filtered <- igraph::delete_edges(mst, igraph::E(mst)[drop])
  clust <- igraph::components(filtered)
  clust.labels <- clust$membership
  clust.size <- clust$csize
  disconnected <- which(clust.size[clust.labels] < min.size)
  clust.labels[disconnected] <- 0
  clusterOld <- c(0, unique(clust.labels[clust.labels != 0]))
  clust.labels <- match(clust.labels, clusterOld) - 1
  return(clust.labels)
}

#' Calculate a path between a source and target with the minimum likelihood
#' of observing any links under the null hypothesis
#'
#' @param s Name (or index) of the source node
#' @param t Name (or index) of the target node
#' @param pmat A p-value matrix
#' @param r A regularization constant, the log of which will be added to all
#' distances
#' @param clusters A vector of cluster labels. If null, all nodes will be
#' the same color
#' @param nodeNames A vector of node names. If null, indices will be used instead
#' @param plotpath If TRUE, plot the resulting path
#' @return The minimal path connecting the source to target
#' @examples
#' pm <- pMatrix(foldChange, traitData, unsigned=TRUE, output="orig")
#' p <- minPath(1, 2, pm, r=1-1/dim(foldChange)[1])
minPath <- function(s, t, pmat,
                    r = 1,
                    clusters=NULL,
                    nodeNames = NULL,
                    plotpath = TRUE) {

  if (is.null(nodeNames)) {
    nodeNames <- seq_len(dim(pmat)[1])
  }
  if (is.null(clusters)) {
    clusters <- rep(1, dim(pmat)[1])
  }

  s <- which(nodeNames == s)
  t <- which(nodeNames == t)

  n_phen <- dim(pmat)[1] - length(clusters)
  odist <- -log1p(-pmat)
  temp_pmat <- odist - log(r)
  diag(temp_pmat) <- 0
  g <- igraph::graph_from_adjacency_matrix(temp_pmat, weighted = TRUE)
  sp <- igraph::shortest_paths(g, s, t, output = "both")

  pvals <- igraph::E(g)[sp$epath[[1]]]$weight
  sp_vert <- sp$vpath[[1]]
  nodeNames <- nodeNames[sp_vert]
  nodeClusters <- c(rep(0, n_phen), clusters)[sp_vert] + 1
  if (plotpath) {
    .plotPath(sp_vert, clusters, nodeNames, nodeClusters)
  }
  return(list(
    pvals = pvals, clusters = nodeClusters - 1,
    names = nodeNames
  ))
}

.plotPath <- function(sp_vert, clusters, nodeNames, nodeClusters) {
  gdf <- data.frame(from = seq_len(length(sp_vert) - 1),
                    to = c(2:length(sp_vert)))
  gdf <- igraph::graph.data.frame(gdf)


  hpal <- c("white", scales::hue_pal()(max(clusters)))
  spal <- c("none", rep("circle", max(clusters)))


  plot(gdf,
       layout = igraph::layout_as_tree, vertex.label = nodeNames,
       vertex.color = hpal[nodeClusters],
       edge.arrow.size = 0.5, vertex.shape = spal[nodeClusters]
  )
}
#' Compute a bootstrap distribution for testing minimal path lengths
#'
#' @param x A data matrix
#' @param partial Identical to the "partial" parameter of pMatrix
#' @param unsigned Identical to the "unsigned" parameter of pMatrix
#' @param output Identical to the "output" parameter of pMatrix
#' @return A vector of path lengths obtained from a permuted data matrix
#' @examples
#' db <- distBoot(foldChange, partial=TRUE)
distBoot <- function(x, partial = FALSE, unsigned = TRUE, output = "orig") {
  n <- dim(x)[2]
  xShuf <- t(apply(x, 1, sample))
  pmatShuf <- pMatrix(xShuf,
                    unsigned = unsigned, partial = partial, output = output)
  D <- -log1p(-pmatShuf)
  diag(D) <- 0

  d <- Rfast::floyd(D)
  du <- d[upper.tri(d)]
  return(du)
}

#' Test a given path length connecting two traits based on a bootstrapped
#' path length distribution of path lengths between genes
#'
#' @param d A vector of path lengths obtained from distBoot
#' @param z The path length being tested
#' @return The p-value testing the shortness of the provided path-length
#' @examples
#' db <- distBoot(foldChange, partial=TRUE)
#' p <- pBoot(db, 0.016)
pBoot <- function(d, z) {
  delta <- z - d
  p <- pgamma(delta[delta > 0], 2, lower.tail = FALSE, log.p = TRUE)
  log.p <- sum(p)
  return(1 - exp(log.p))
}


