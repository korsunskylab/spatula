## Useful function for creating a bbox of type POLYGON
#' @export 
st_bbox_polygon <- function(pts) {
    bbox <- st_bbox(pts)
    coords <- bbox[c('xmin', 'xmax', 'ymin', 'ymax')]
    res <- st_polygon(list(
        rbind(
            c(coords[1], coords[3]),
            c(coords[2], coords[3]), 
            c(coords[2], coords[4]), 
            c(coords[1], coords[4]), 
            c(coords[1], coords[3]))
    ))
    return(res)
}

#' @export 
st_rectangle <- function(xmin, xmax, ymin, ymax) {
    res <- st_polygon(list(
        rbind(
            c(xmin, ymin),
            c(xmax, ymin), 
            c(xmax, ymax), 
            c(xmin, ymax), 
            c(xmin, ymin))
    ))
    return(res)    
}

#' @export 
infer_polygon <- function(x, y, max_iter=20, verbose=FALSE) {
    if (verbose) plot(cbind(x, y))
    
    ## (1) Kernel Density Estimation
    kde_res <- try({
        # R.utils::withTimeout({
        MASS::kde2d(x, y, n = 100)
        # }, timeout = max_time, onTimeout = "error")        
    })
    if (inherits(kde_res, 'try-error')) {
        return(sf::st_polygon()) ## return empty polygon             
    } 
    # kde_res <- MASS::kde2d(x, y, n = 100, lims = c(min(x) * 1, max(x) * 2, min(y) * 1, max(y) * 1))  
    # if (verbose) plot(st_as_stars(kde_res$z))

    ## Find MST with appropriate kernel density thresholding 
    ## NOTE: a threshold that is too high will give multiple connected components 
    ## NOTE: when cell has multiple compartments (i.e. densities), this should be OK! 

    ## Initial kde quantile estimate says that ~50% of the image should be in the cell
    ## ALTERANATIVE STRATEGY TO TRY: interpolate density for each transcript and choose area to capture 95% of transcripts 
    kde_quantile <- 0.5
    iter <- 1
    while (TRUE) {
        if (verbose) print(glue('KDE quantile threshold = {kde_quantile}'))
        ## (2) Choose kernel density for segmentation cutoff 
        kde_level <- quantile(as.numeric(kde_res$z), kde_quantile)

        ## (3) Extract boundary pixels of density 
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
        adj <- Matrix::sparseMatrix(
            i = rep(1:N, each = k), 
            j = c(t(nn_idx)), 
            x = c(t(nn_dist)),
            dims = c(N, N)
        )
        adj <- adj + Matrix::t(adj)
        g <- igraph::graph_from_adjacency_matrix(adj, weighted = TRUE, mode = "undirected")
        tree <- try({
            # R.utils::withTimeout({
             igraph::mst(g)  
            # }, timeout = max_time, onTimeout = "error")        
        })
        if (inherits(tree, 'try-error')) {
            return(sf::st_polygon()) ## return empty polygon  
        }
        if (igraph::is_connected(tree)) {
            break
        }
        else {
            kde_quantile <- kde_quantile * 0.8
        } 
        iter <- iter + 1
        if (iter == max_iter) {
            return(sf::st_polygon()) ## return empty polygon              
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
                    

#' @export 
#' @param tx Data frame with columns x, y, and cell (cell=0 encodes background) 
infer_polygons <- function (tx, colname_x, colname_y, min_tx_per_cell = 5, parallel = FALSE) {
    if (parallel) {
        future::plan(future::multicore)
        iter_fxn <- function(...) {
            furrr::future_imap(..., .options = furrr::furrr_options(seed = TRUE))
        }
    }
    else {
        iter_fxn <- purrr::imap
    }
    cells <- tx %>% subset(cell != 0) %>% 
        split(.$cell) %>% 
        iter_fxn(function(tx_cell, idx) {
        shape <- tryCatch({
            if (nrow(tx_cell) < min_tx_per_cell) {
                sf::st_polygon()
            }
            else {
                spatula::infer_polygon(tx_cell[[colname_x]], tx_cell[[colname_y]], 
                  max_iter = 10)
            }
        }, error = function(e) {
            sf::st_polygon()
        })
        return(list(id = idx, shape = shape))
    })
    res <- st_sf(id = unlist(map(cells, "id")), shape = st_sfc(map(cells, 
        function(x) x[["shape"]])))
    return(res)
}
                                                                   
#' @export 
st_make_neighborhoods <- function(shapes, dist) {
    stopifnot(is(shapes, 'sfc'))
    neighborhoods <- st_buffer(shapes, dist = dist) 
    neighborhoods <- map2(neighborhoods, shapes, function(x, y) {
        tryCatch({
            st_difference(x, y)
        }, error = function(e) {
            st_difference(st_make_valid(x), st_make_valid(y)) ## try to save this polygon 
        }, finally = function(e) {
            st_polygon() ## return empty (in case of self intersections, which are pretty rare)                
        })
    }) %>% 
        st_sfc()
    return(neighborhoods)
}

#' @export 
st_make_neighborhoods_multi <- function(shapes, dist, nsplit) {
    stopifnot(is(shapes, 'sfc'))
    if (nsplit == 1) 
        return(st_make_neighborhoods(shapes, dist))
    
    plan(multicore)
    split_id <- rep(1:nsplit, each = round(length(shapes) / nsplit) + 1)[1:length(shapes)]
    neighborhoods <- split(shapes, split_id) %>% 
        future_imap(function(shapes_tile, idx) {
            st_sf(
                shape = st_make_neighborhoods(shapes_tile, dist), 
                id = as.integer(idx) ## this is to maintain the original order of the shapes 
            )  
        }, .options = furrr_options(seed = 42L)) %>% 
        bind_rows() %>% 
        arrange(id) %>% 
        with(shape) %>% 
        st_sfc()
    return(neighborhoods)
}

                                                                   
                     