## Compute adjacency matrix of spatial nearest neighbors using Delauney triangulation
## coords is an Nx2 matrix of coordinates 
#' @export 
getSpatialNeighbors <- function(coords, dist_thresh_quantile=.99) {
    ## do Delaunay triangulation to estimate neighbors
    triplets <- geometry::delaunayn(coords) 

    ## expand list of triplets into symmetric list of pairs
    pairs <- triplets_to_pairs(triplets) 

    ## prune very large distances
    dists <- sqrt(rowSums((coords[pairs[, 1], ] - coords[pairs[, 2], ])^2))
    idx_keep <- which(dists < quantile(dists, dist_thresh_quantile)) 
    pairs <- pairs[idx_keep, ]
    
    ## convert list of edges to sparse adjacency matrix
    adjmat <- Matrix::sparseMatrix(
        i = pairs[, 1], 
        j = pairs[, 2],
        x = rep(1, nrow(pairs))
    ) 

    ## accounts for duplicate pairs - faster than `unique(pairs)`
    adjmat@x <- pmin(adjmat@x, 1) 
    return(adjmat)        
}
