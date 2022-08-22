## Useful function for creating a rectangle from corner points
#' @export 
st_rectangle <- function(coords) {
    st_polygon(list(rbind(c(coords[1], coords[3]), c(coords[2], coords[3]), c(coords[2], coords[4]), c(coords[1], coords[4]), c(coords[1], coords[3]))))    
}

#' @export 
infer_polygon <- function(x, y, verbose=FALSE) {
    if (verbose) plot(cbind(x, y))
    
    ## (1) Kernel Density Estimation
    kde_res <- MASS::kde2d(x, y, n = 100)  
    # kde_res <- MASS::kde2d(x, y, n = 100, lims = c(min(x) * 1, max(x) * 2, min(y) * 1, max(y) * 1))  
    if (verbose) plot(st_as_stars(kde_res$z))

    ## Find MST with appropriate kernel density thresholding 
    ## NOTE: a threshold that is too high will give multiple connected components 
    ## NOTE: when cell has multiple compartments (i.e. densities), this should be OK! 

    ## Initial kde quantile estimate says that ~50% of the image should be in the cell
    ## ALTERANATIVE STRATEGY TO TRY: interpolate density for each transcript and choose area to capture 95% of transcripts 
    kde_quantile <- 0.5
    while (TRUE) {
        if (verbose) print(glue('KDE quantile threshold = {kde_quantile}'))
        ## (2) Choose kernel density for segmentation cutoff 
        kde_level <- quantile(as.numeric(kde_res$z), kde_quantile)
        # if (verbose) hist(as.numeric(kde_res$z))

        ## (3) Extract boundary pixels of density 
        # cell <- extract_cells(kde_res$z > kde_level, parallel = FALSE, TRUE, cbind(kde_res$x, kde_res$y))[[1]]
        # plot(cell)
        boundary_points <- get_idx(kde_res$z > kde_level, TRUE)[[1]]
        boundary_points[, 1] <- kde_res$x[boundary_points[, 1]]
        boundary_points[, 2] <- kde_res$y[boundary_points[, 2]]
        if (verbose) plot(boundary_points)

        ## (4) Minimum Spanning Tree
        if (nrow(boundary_points) < 11) {
            return(st_polygon())
        }
        k <- 10
        nn_res <- RANN::nn2(boundary_points, k = k + 1, eps = 0)
        nn_idx <- nn_res$nn.idx[, 2:(k + 1)]
        nn_dist <- nn_res$nn.dists[, 2:(k + 1)]
        N <- nrow(boundary_points)
        adj <- Matrix::sparseMatrix(i = rep(1:N, each = k), j = c(t(nn_idx)), x = c(t(nn_dist)))
        adj <- adj + Matrix::t(adj)
        g <- igraph::graph_from_adjacency_matrix(adj, weighted = TRUE, mode = "undirected")
        tree <- igraph::mst(g)
        
        if (igraph::is_connected(tree)) {
            ## If your tree is connected, you're done! 
            break 
        } else {
            ## Otherwise, relax the quantile threshold for density map 
            kde_quantile <- kde_quantile * .8
        }
    }
    
    ## (5) Longest Path in Tree (with Dijkstra's SSSP)
    ## Heuristic: only check paths among roots 
    leaves <- which(igraph::degree(tree) == 1)
    all_paths <- map(leaves, function(leaf) igraph::shortest_paths(tree, leaf, setdiff(leaves, leaf))$vpath)
    all_paths <- Reduce(c, c(all_paths))
    node_order <- as.integer(all_paths[[which.max(map_dbl(all_paths, length))]])
    node_order <- c(node_order, node_order[1])
    res <- sf::st_polygon(list(boundary_points[node_order, ]))
    res <- st_simplify(res)
    return(res)
}

