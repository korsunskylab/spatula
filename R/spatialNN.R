## Compute adjacency matrix of spatial nearest neighbors using Delauney triangulation
## coords is an Nx2 matrix of coordinates 
#' @export 
spatialNN <- function(coords, dist_thresh=Inf, do_soi=TRUE) {
    ## Initial sf object 
    coords_sf <- tibble(geom = map2(coords[, 1], coords[, 2], function(.r, .c) st_point(c(.r, .c)))) %>% st_as_sf()
                                    
    ## Do Delauney triangulation
    nb_tri <- coords_sf %>% st_geometry() %>% tri2nb()
                                    
    ## (Optionally) prune graph using Sphere of Influence method
    if (do_soi) {
        nb_tri <- soi.graph(
            tri.nb = nb_tri, 
            coords = coords_sf %>% st_geometry() %>% unlist() %>% matrix(ncol = 2, byrow = TRUE)
        ) %>% 
        graph2nb()        
    }
                                    
    ## Convert to dgCMatrix adjacency matrix
    wb_tri <- spdep::nb2WB(nb_tri)
    mat <- Matrix::sparseMatrix(
        i = wb_tri$adj,
        p = c(0, cumsum(wb_tri$num)), 
        x = unlist(nbdists(nb_tri, coords))
    )
                                    
    ## (Optionally) prune graph on max distance
    mat@x[mat@x > dist_thresh] <- 0
    mat <- Matrix::drop0(mat)
    mat@x <- rep(1, length(mat@x)) ## return unweighted matrix 
    rownames(mat) <- colnames(mat) <- rownames(coords)
    return(mat)
}
