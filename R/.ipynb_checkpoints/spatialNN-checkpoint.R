## Compute adjacency matrix of spatial nearest neighbors using Delauney triangulation
## coords is an Nx2 matrix of coordinates 
#' @export 
getSpatialNeighbors <- function(coords, dist_thresh_quantile=.99, return_weights=FALSE) {
    ## do Delaunay triangulation to estimate neighbors
    triplets <- geometry::delaunayn(coords) 

    ## expand list of triplets into symmetric list of pairs
    pairs <- triplets_to_pairs(triplets) 
    pairs <- unique(pairs)

    ## prune very large distances
    dists <- sqrt(rowSums((coords[pairs[, 1], ] - coords[pairs[, 2], ])^2))
    if (is.infinite(dist_thresh_quantile)) {
        idx_keep <- seq_len(nrow(pairs))
    } else {
        idx_keep <- which(dists < quantile(dists, dist_thresh_quantile)) 
    }
    
    if (return_weights == FALSE) {
        dists <- rep(1, nrow(pairs))
    }
    
    ## convert list of edges to sparse adjacency matrix
    adjmat <- Matrix::sparseMatrix(
        i = pairs[idx_keep, 1], 
        j = pairs[idx_keep, 2],
        x = dists[idx_keep],
        dims = c(nrow(coords), nrow(coords))
    ) 

    # ## accounts for duplicate pairs - faster than `unique(pairs)`
    # adjmat@x <- pmin(adjmat@x, 1) 
    return(adjmat)        
}
