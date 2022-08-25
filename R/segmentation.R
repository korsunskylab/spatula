#' @export 
initialize_grid <- function(img_fname, section_size=2048) {
    img_dim <- dim(read_stars(img_fname)) ## use any image to get coordinates
    img_dim
    ## CAUTION: grid is 0-indexed
    xmax <- img_dim['x']
    ymax <- img_dim['y']

    ## make sure that the grid covers the whole space! 
    ## without this, grid will leave off edges 
    ymax <- ymax + (section_size - ymax %% section_size)
    xmax <- xmax + (section_size - xmax %% section_size)

    ## define the bounding box of the full image 
    bbox <- st_sfc(st_polygon(list(rbind(c(0, 0), c(xmax, 0), c(xmax, ymax), c(0, ymax), c(0, 0)))))
    grid <- st_make_grid(bbox, cellsize = section_size) 
    
    ## Create overlaps between grids 
    overlap_size = round(0.10 * section_size)
    grid <- grid %>% 
        map(function(.grid) {
            .grid <- as.matrix(.grid)
            idx_x0 <- which(.grid[, 1] == min(.grid[, 1]))
            idx_x1 <- which(.grid[, 1] == max(.grid[, 1]))
            .grid[idx_x0, 1] <- pmax(.grid[idx_x0, 1] - overlap_size, 0)
            .grid[idx_x1, 1] <- pmin(.grid[idx_x1, 1] + overlap_size, xmax)
            idx_y0 <- which(.grid[, 2] == min(.grid[, 2]))
            idx_y1 <- which(.grid[, 2] == max(.grid[, 2]))
            .grid[idx_y0, 2] <- pmax(.grid[idx_y0, 2] - overlap_size, 0)
            .grid[idx_y1, 2] <- pmin(.grid[idx_y1, 2] + overlap_size, ymax)
            return(st_polygon(list(.grid)))
        }) %>% 
        sf::st_as_sfc()
    return(grid)
}

#' @export 
read_subimgs <- function(fnames, subimg_coords, projection_mode='max', do_parallel=TRUE) {
    if (do_parallel) {
        plan(multicore)
        iter_fxn <- function(...) furrr::future_map(..., .options = furrr_options(seed=TRUE))
    } else {
        iter_fxn <- purrr::map
    }
    res <- iter_fxn(fnames, function(fname) {
        x0 <- max(min(subimg_coords[, 1]), 1)
        y0 <- max(min(subimg_coords[, 2]), 1)
        img_dim <- dim(read_stars(fname))
        x1 <- min(max(subimg_coords[, 1]), img_dim[['x']])
        y1 <- min(max(subimg_coords[, 2]), img_dim[['y']])
        rasterio = list(
            nXOff = x0, 
            nYOff = y0, 
            nXSize = x1 - x0 + 1, 
            nYSize = y1 - y0 + 1
            
            ## NOTE: use code below for downsampling 
            # nBufXSize = ceiling((1/downsample_int) * (x1 - x0 + 1)), 
            # nBufYSize = ceiling((1/downsample_int) * (y1 - y0 + 1)) 
        ) 
        img <- read_stars(fname, RasterIO = rasterio, proxy = FALSE)
        return(img)
    }) 
    if (projection_mode == 'max') {
        res <- res %>% purrr::map(1) %>% purrr::reduce(pmax) 
    } else if (projection_mode == 'sum') {
        res <- purrr::reduce(res, `+`) 
    } else {
        stop('Invalid projection_mode')
    }
    return(res)
}


#' @export 
do_mesmer <- function(
    .grid, 
    mesmer_app, 
    fnames_nucleus, 
    fnames_cytoplasm, 
    img_scale, 
    projection_mode='max', 
    mesmer_mode=c('whole-cell', 'both')[1]
) {
    ## Read subimages and collapse them 
    # message('Read images')
    .grid <- as.matrix(.grid)
    subimg_nucleus <- read_subimgs(fnames_nucleus, .grid, projection_mode)
    subimg_cytoplasm <- read_subimgs(fnames_cytoplasm, .grid, projection_mode)
    
    ## Segment cells
    # mesmer_res <- do_mesmer(subimg_nucleus, subimg_cytoplasm, img_scale) 
    ## Prepare image with correct dimensions
    img_mat <- array(dim = c(1, nrow(subimg_nucleus), ncol(subimg_nucleus), 2))
    if (is(subimg_nucleus, 'stars')) subimg_nucleus <- img_nucleus[[1]]
    if (is(subimg_cytoplasm, 'stars')) subimg_cytoplasm <- subimg_cytoplasm[[1]]
    img_mat[1, , , 1] <- subimg_nucleus
    img_mat[1, , , 2] <- subimg_cytoplasm
    
    ## Run Mesmer prediction
    # message('Run mesmer')
    mesmer_res <- mesmer_app$predict(img_mat, compartment=mesmer_mode, image_mpp=img_scale)

    ## Extract cells and (corresponding) nuclei from Mesmer result 
    # message('Extract polygons')
    cells <- seg_img_to_sf(mesmer_res[1, , , 1], .grid)
    if (mesmer_mode == 'both') {
        nuclei <- seg_img_to_sf(mesmer_res[1, , , 2], .grid)
        nuclei <- align_cells_nuclei(cells, nuclei)         
        return(list(cells=cells, nuclei=nuclei, boundary=st_polygon(list(.grid))))
    } else {
        return(list(cells=cells, boundary=st_polygon(list(.grid))))
    }
    
}


#' @export 
get_idx <- function(img, only_boundary) {
    res <- get_idx_cpp(img, only_boundary)
    res <- purrr::map(res, matrix, ncol=2, byrow=TRUE)
    return(res)
}

#' @export 
trace_outline <- function(pixels_mat) {   
    if (nrow(pixels_mat) < 11) {
        return(st_polygon())
    }
    k <- 10
    nn_res <- RANN::nn2(pixels_mat, k = k+1, eps = 0)
    nn_idx <- nn_res$nn.idx[, 2:(k+1)]
    nn_dist <- nn_res$nn.dists[, 2:(k+1)]
    N <- nrow(pixels_mat)
    adj <- Matrix::sparseMatrix(i = rep(1:N, each = k), j = c(t(nn_idx)), x = c(t(nn_dist)))
    # adj <- Matrix::sparseMatrix(i = rep(1:N, each = k), j = c(t(nn_idx)), x = rep(1, N * k))
    idx_remove <- which(adj[1, ] > 0)[1]
    adj[1, idx_remove] <- 0
    adj[idx_remove, 1] <- 0
    adj <- Matrix::drop0(adj)
    adj <- adj + Matrix::t(adj) ## make symmetric
    g <- igraph::graph_from_adjacency_matrix(adj, weighted = TRUE, mode = 'undirected')
    tree <- igraph::mst(g)
    o <- as.integer(igraph::bfs(tree, root = 1)$order)
    o <- c(o, o[1])
    res <- sf::st_polygon(list(pixels_mat[o, ]))
    res <- st_simplify(res)
    return(res)
}



#' @export 
extract_cells <- function(img, parallel=FALSE, exclude_interior=TRUE) {
    ## Extract cell pixels lists from segmentation label matrix
    cells_idx <- get_idx(as.matrix(img), exclude_interior) 
    # cells_idx <- get_idx(as.matrix(fread('mesmer_res/grid_1.txt')), TRUE) 

    if (parallel) {
        plan(multicore)
        cells <- sf::st_sfc(future_map(cells_idx, trace_outline))        
    } else {
        cells <- sf::st_sfc(map(cells_idx, trace_outline))        
    }
    return(cells)
    # ## Assign colors to cells for plotting 
    # cells_sf <- sf::st_sf(ID = paste0('Cell', seq_len(length(cells))), geometry = cells)
    # return(cells_sf)
}


#' @export 
seg_img_to_sf <- function(segmentation_img, .grid) {
    cells <- extract_cells(segmentation_img)

    ## NOT NEEDED IF WE STAY IN R?    
#     ## Convert to global coordinates
#     ##   Rotate each by 90 degrees
#     # rotation_matrix <- matrix(c(0, -1, 1, 0), nrow=2, byrow=TRUE)
#     rotation_matrix <- matrix(c(1, 0, 0, 1), nrow=2, byrow=TRUE) ## Identity 
#     # rotation_matrix <- matrix(c(0, 1, 1, 0), nrow=2, byrow=TRUE)
#     cells <- cells %>% 
#         map(as.matrix) %>% 
#         map(`%*%`, rotation_matrix) %>% 
#         map(list) %>% map(sf::st_polygon) %>% sf::st_sfc()

#     ##   translate it back to left corner = c(0,0)
#     coords_to_shift <- cells %>% map(as.matrix) %>% map(apply, 2, min) %>% purrr::reduce(pmin)
#     cells <- cells %>% map(as.matrix) %>%
#         map(sweep, 2, coords_to_shift - 1, '-') %>% 
#         map(list) %>% map(sf::st_polygon) %>% sf::st_sfc() 

    ##  offset by global grid coordinates
    x_offset <- min(.grid[, 1])
    y_offset <- min(.grid[, 2])
    cells <- cells + c(x_offset, y_offset)

    return(cells)
}

#' @export 
align_cells_nuclei <- function(cells, nuclei) {
    sf_nuclei <- st_sf(ID_nucleus = paste0('N', 1:length(nuclei)), geometry = nuclei)
    sf_cells <- st_sf(ID_cell = paste0('C', 1:length(cells)), geometry = cells)
    suppressWarnings({
        sf_joint <- sf_cells %>% 
            ## Assign up to 1 nucleus to each cell
            st_join(sf_nuclei,left = TRUE, largest = TRUE) %>% 
            ## Make sure that each nucleus is assigned to only 1 cell
            left_join(
                as.data.frame(st_join(sf_nuclei, sf_cells, left = FALSE, largest = TRUE)),
                by = c('ID_cell')
            ) %>% 
            subset(
                is.na(ID_nucleus.x) | ID_nucleus.x == ID_nucleus.y
            ) %>% 
            dplyr::select(-ID_nucleus.y) %>% 
            dplyr::rename(ID_nucleus = ID_nucleus.x, geometry_cell = geometry.x, geometry_nucleus = geometry.y) %>% 
            identity()
    })
    
    ## Sometimes, nucleus is actually bigger than cell. 
    ## This is clearly a mistake that needs to be addressed but ...
    ## For now, just crop the nucleus to fit inside the cell
    # sf_joint$geometry_nucleus <- st_sfc(map2(sf_joint$geometry_nucleus, sf_joint$geometry_cell, st_intersection))
    # return(list(cells = sf_joint$geometry_cell, nuclei = sf_joint$geometry_nucleus))
    nuclei_aligned <- st_sfc(map2(sf_joint$geometry_nucleus, sf_joint$geometry_cell, st_intersection))
    return(nuclei_aligned)
}



#' @export 
do_segment <- function(fnames_nucleus, fnames_cytoplasm, img_scale, section_size=2048, projection_mode='max') {
    ## (1) Make grid
    grid <- initialize_grid(fnames_nucleus[[1]], section_size)
    
    ## (2) Segment cells in each grid section
    deepcell <- reticulate::import('deepcell')
    mesmer_app <- deepcell$applications$Mesmer()
    cells_list <- imap(
        grid, 
        # grid[c(1001:1002)], 
        function(.grid, i) {
            message(i) ## This is just to keep track of progress
            do_mesmer(.grid, mesmer_app, fnames_nucleus, fnames_cytoplasm, img_scale, 'max')
        }
    ) 

    ## (3) Stitch the segmented images together into one 
    boundary_size <- round(0.1 * section_size)
    cells_registered <- Reduce(function(x, y) register_labels(x, y, boundary_size), cells_list)
    return(list(cells_list=cells_list, cells_registered=cells_registered))    
}



#' @export 
register_labels <- function(cells1, cells2, boundary_size) {
    ## Registers two sets of geometries with overlapping 
    plan(multicore)
    
    ## Create sf tables 
    sf1 <- sf::st_sf(
        ID = paste0('Cell', seq_len(length(cells1$cells)), '_1'), 
        geometry_cell = cells1$cells
    )
    sf2 <- sf::st_sf(
        ID = paste0('Cell', seq_len(length(cells2$cells)), '_2'), 
        geometry_cell = cells2$cells
    )
    sf1$edge_score <- compute_edge_scores(cells1$cells, cells1$boundary, boundary_size)
    sf2$edge_score <- compute_edge_scores(cells2$cells, cells2$boundary, boundary_size)

    do_nuclei <- FALSE
    if (!is.null(cells1$nuclei) & !is.null(cells2$nuclei)) {
        sf1$geometry_nucleus <- cells1$nuclei
        sf2$geometry_nucleus <- cells2$nuclei
        do_nuclei <- TRUE
    }

    ## get all the pairs of overlapping cells 
    ## Limit overlap analysis in overlap of grid 
    pairs <- get_overlapping_cells(sf1, sf2, cells1$boundary, cells2$boundary)

    ## Check for case of no overlap, Return union of cells
    if (length(pairs) == 0) {
        new_cells <- list(
            cells = c(cells1$cells, cells2$cells),
            boundary = st_union(cells1$boundary, cells2$boundary)
        )
        return(new_cells) 
    }
        
    pairs <- cbind(
        dplyr::select(sf1[map_int(pairs, 2), ], ID1 = ID, geometry1 = geometry_cell, edge_score1 = edge_score),
        dplyr::select(sf2[map_int(pairs, 1), ], ID2 = ID, geometry2 = geometry_cell, edge_score2 = edge_score)
    ) %>% 
        data.frame() 

    ## Prune pairs to have some minimal overlap 
    pairs$pct_overlap <- future_map2_dbl(
        pairs$geometry1, pairs$geometry2, 
        function(g1, g2) {
            area_overlap <- st_intersection(g1, g2) %>% st_area()
            area_overlap / min(st_area(g1), st_area(g2))
        },
        .options = furrr_options(seed=TRUE)   
    )
    
    ## TODO: filter for overlap above. We shouldn't check for empty pairs twice! 
    pairs <- subset(pairs, pct_overlap > .1) ## 20% overlap 
    if (nrow(pairs) == 0) {
        new_cells <- list(
            cells = c(cells1$cells, cells2$cells),
            boundary = st_union(cells1$boundary, cells2$boundary)
        )
        return(new_cells) 
    }
    
    
    
    ## Divide the bipartite graph into connected components and deal with each one separately 
    g <- pairs %>% 
        with(rbind(ID1, ID2)) %>% 
        igraph::graph(directed = FALSE)
    connected_components <- components(g)
    connected_components <- split(names(connected_components$membership), connected_components$membership)

    ## one-to-one relationships
    ## this represents agreeing segmentations 
    ## so the action is to merge them into new cells (with potential overlaps)
    one_to_ones <- connected_components[as.integer(which(map_int(connected_components, length) == 2))]
    new_cells_merged <- future_map2(
        sf1$geometry_cell[as.integer(gsub('Cell(\\d+)_1', '\\1', map_chr(one_to_ones, 1)))],
        sf2$geometry_cell[as.integer(gsub('Cell(\\d+)_2', '\\1', map_chr(one_to_ones, 2)))],
        st_union,
        .options = furrr_options(seed=TRUE)   
    ) %>% st_sfc()
    
    if (do_nuclei) {
        new_nuclei_merged <- map(one_to_ones, function(.cells) {
            res <- st_union(
                subset(as.data.frame(sf1), ID == .cells[[1]])$geometry_nucleus,
                subset(as.data.frame(sf2), ID == .cells[[2]])$geometry_nucleus
            ) %>% st_coordinates()
            ## Because nuclei can be missing in cells, handle edge case with empty polygon
            if (nrow(res) > 0) {
                res <- st_polygon(list(res[, 1:2])) 
            } else {
                res <- st_polygon()
                # res <- st_geometrycollection()
            }
            return(res)
        }) %>% st_sfc() 
    }

    ## other cases (not one-to-one)
    ## this represents conflicting segmentations 
    ## so the course of action is to choose cells from one image over those from another 
    other_components <- connected_components[as.integer(which(map_int(connected_components, length) > 2))]
    cells_kept_sf <- map(other_components, function(.cells) {
        cells1 <- grep('_1$', .cells, value = TRUE)
        cells2 <- grep('_2$', .cells, value = TRUE)
        
        ## keep cells from the image where cells are less on the border
        if (
            mean(subset(as.data.frame(sf1), ID %in% cells1)$edge_score) < 
            mean(subset(as.data.frame(sf2), ID %in% cells2)$edge_score)
        ) {
            subset(sf1, ID %in% cells1)
        } else {
            subset(sf2, ID %in% cells2)
            # subset(as.data.frame(sf2), ID %in% cells2)
        }  
    }) %>% 
        dplyr::bind_rows() 

    ## With all boundary regions resolved
    cells <- c(
        subset(as.data.frame(sf1), !ID %in% pairs$ID1)$geometry_cell,
        subset(as.data.frame(sf2), !ID %in% pairs$ID2)$geometry_cell,
        cells_kept_sf$geometry_cell,
        new_cells_merged
    )
    
    if (do_nuclei) {
        nuclei <- c(
            subset(as.data.frame(sf1), !ID %in% pairs$ID1)$geometry_nucleus,
            subset(as.data.frame(sf2), !ID %in% pairs$ID2)$geometry_nucleus,
            cells_kept_sf$geometry_nucleus,
            new_nuclei_merged
        )
        return(list(cells=cells, nuclei=nuclei, boundary=st_union(cells1$boundary, cells2$boundary)))
    } else {
        return(list(cells=cells, boundary=st_union(cells1$boundary, cells2$boundary)))        
    }
}


#' @export 
compute_edge_scores <- function(cells, grids, boundary_size) {
    if (is(grids, 'POLYGON')) {
        .grid_lines <- st_multilinestring(grids) ## POLYGON    
    } else {
        .grid_lines <- st_sfc(map(grids, st_multilinestring)) ## MULTIPOLYGON
    }
    dists <- apply(sf::st_distance(cells, .grid_lines), 1, min)
    res <- exp(-as.numeric(dists) / boundary_size)
    return(res)
}





#' @export 
get_overlapping_cells <- function(sf1, sf2, grid1, grid2) {
    ## get all the pairs of overlapping cells 
    ## Limit overlap analysis in overlap of grid 
    grid_overlap <- st_intersection(grid1, grid2)
    if (st_is_empty(grid_overlap)) {
        return(list())
    }
    
    suppressWarnings({
        ids_1 <- sf1 %>% st_crop(grid_overlap) %>% rownames(1) %>% as.integer
        ids_2 <- sf2 %>% st_crop(grid_overlap) %>% rownames(1) %>% as.integer
    })
    
    overlaps <- st_intersects(
        sf1$geometry_cell[ids_1], 
        sf2$geometry_cell[ids_2], 
        sparse = TRUE
    )
    ## Remove empty overlaps 
    overlaps <- overlaps[map_int(overlaps, length) > 0]    
    if (length(overlaps) == 0) {
        return(list())
    }

    pairs <- which(map_int(overlaps, length) > 0) %>% 
        map(function(i) {
            map(ids_2[overlaps[[i]]], c, ids_1[i])
            # map(overlaps[[i]], c, i)
        }) %>% 
        purrr::reduce(append)
    return(pairs)
}

#' @export 
register_labels_multi <- function(.cells_list, boundary_size) {
    ## Function to decide order of pairwise registration of adjacent images
    while (length(.cells_list) > 1) {
        ##   Choose any two images of smallest size 
        ncells <- .cells_list %>% map('cells') %>% map_dbl(length)
        idx <- order(ncells)[1:2]
        new_cells <- register_labels(.cells_list[[idx[[1]]]], .cells_list[[idx[[2]]]], boundary_size)
        
        .cells_list <- .cells_list[setdiff(seq_len(length(.cells_list)), idx)]
        .cells_list[[length(.cells_list) + 1]] <- new_cells
    }    
    return(.cells_list[[1]])
}


