## Functions for quantifying things about cells 
## Analogous to regionprops in skimage 

#' @export 
st_assign_pt <- function(pts, ...) {
    UseMethod('st_assign_pt')
}

#' @export 
st_assign_pt.default <- function(pts, shapes_sf, shape_id_col) {
    stop(paste0('st_assign_pt not supported for pts of type ', class(pts)))
}

#' @export 
st_assign_pt.data.frame <- function(pts, ...) {
    st_assign_pt(data.table(pts), ...)
}


#' @export 
st_assign_pt.sf <- function(pts, shapes_sf, shape_id_col) {
    pts$cell <- st_assign_pt.sfc(st_geometry(pts), shapes_sf, shape_id_col)
    return(pts)
}

#' @export 
st_assign_pt.sfc <- function (pts, shapes_sf, shape_id_col) {
    stopifnot(is(shapes_sf, "sf"))
    stopifnot(shape_id_col %in% colnames(shapes_sf))
    if (nrow(shapes_sf) == 0) return(rep(0L, length(pts)))        
    
    suppressWarnings({
        ## NOTE: use st_buffer to handle edge case of bbox with no area 
        bbox <- st_intersection(
            st_buffer(st_as_sfc(st_bbox(pts)), .1), 
            st_buffer(st_as_sfc(st_bbox(shapes_sf)), .1)
        )
        if (length(bbox) == 0) return(rep(0L, length(pts)))        
        # bbox <- st_buffer(bbox, 1e-8) ## needed? 
        
        ## for efficiency, 
        # pts <- st_crop(pts, bbox) ## bad idea! tx in must == tx out 
        tx_in <- st_contains(bbox, pts)[[1]] ## transcripts inside bbox 
        tx_out <- setdiff(seq_len(length(pts)), tx_in)

        if (length(tx_in) == 0) return(rep(0L, length(pts)))        
        
        ## Crop is expensive, intersects much cheaper 
        # shapes_sf <- st_crop(shapes_sf, bbox)
        shapes_sf <- shapes_sf[st_intersects(bbox, shapes_sf)[[1]], ]        
        
        cell_ids <- st_intersects(pts[tx_in], st_geometry(shapes_sf)) %>% 
            map(head, 1) %>% as.integer()
        res <- rep(0L, length(pts))
        res[tx_in] <- cell_ids
        res[tx_in] <- shapes_sf[[shape_id_col]][res[tx_in]]
        res[is.na(res)] <- 0L
    })
    return(res)
}

#' @export 
st_assign_pt.data.table <- function(
    tx_dt, shapes_sf, colname_x, colname_y, gridn=10, shape_id_col='cell', parallel=FALSE, verbose=FALSE
) {
    stopifnot(is(tx_dt, 'data.table'))
    stopifnot(colname_x %in% colnames(tx_dt) & colname_y %in% colnames(tx_dt))
    tx_dt[, ORDER := 1:.N] ## to preserve order 
    
    ## (1) Before we start, limit to common area shared by transcripts and cells
    if (verbose) message('(1) Make bounding boxes')
    bbox_tx <- st_rectangle(
        min(tx_dt[[colname_x]]), max(tx_dt[[colname_x]]), 
        min(tx_dt[[colname_y]]), max(tx_dt[[colname_y]])
    )
    bbox_cells <- st_as_sfc(st_bbox(shapes_sf))
    bbox <- st_intersection(bbox_tx, bbox_cells)
    bbox <- st_buffer(bbox, 1e-8) ## needed? 
    bbox <- st_as_sfc(st_bbox(bbox)) ## simplify bbox 
    bbox_coords <- st_bbox(bbox)

    if (st_is_empty(bbox)) {
        tx_dt$cell <- 0L
        return(tx_dt)
    }

    # Do we really need to do this cropping? grid should do it for us below: 
    # CAUTION: cropping tx directly is dangerous! We want to just append a cell column to the tx
    # Here, we just append a column by reference 
    if (verbose) message('(2.1) crop transcripts')
    tx_dt[
        , INBBOX := (
            eval(parse(text = colname_x)) >= bbox_coords[['xmin']] & 
            eval(parse(text = colname_x)) < bbox_coords[['xmax']] & 
            eval(parse(text = colname_y)) >= bbox_coords[['ymin']] & 
            eval(parse(text = colname_y)) < bbox_coords[['ymax']]         
        )
    ] 
    suppressWarnings({
        ## NOTE: cropping is expensive (b/c of intersection)
        ##       intersect is much cheaper 
        if (verbose) message('(2.2) crop cells')
        shapes_sf <- shapes_sf[st_intersects(bbox, shapes_sf)[[1]], ] 
        # shapes_sf <- st_crop(shapes_sf, bbox)
        grid <- st_make_grid(bbox, n = gridn) 
    })

    ## (2) Split data 
    if (verbose) message('(3.1) Split transcripts')
    tx_list <- map(grid, function(g) {
        bbox_g <- st_bbox(g)
        tx_dt[INBBOX == TRUE][
            eval(parse(text = colname_x)) >= bbox_g[['xmin']] & 
            eval(parse(text = colname_x)) < bbox_g[['xmax']] & 
            eval(parse(text = colname_y)) >= bbox_g[['ymin']] & 
            eval(parse(text = colname_y)) < bbox_g[['ymax']]  
        ]
    })
    
    grid_use <- which(map_int(tx_list, nrow) > 0)
    grid <- grid[grid_use]
    tx_list <- tx_list[grid_use]    
    tx_list <- map(tx_list, st_as_sf, coords = c(colname_x, colname_y))        
    stopifnot(sum(map_int(tx_list, nrow)) == sum(tx_dt$INBBOX)) 
    if (verbose) message('(3.2) Split cells')    
    tile_cells <- st_intersects(grid, shapes_sf)
    shapes_list <- map(tile_cells, function(i) shapes_sf[i, ]) 
    
    ## (3) Do the assignment and put it back together 
    if (verbose) message('(4) do assignment')
    if (parallel == TRUE) {
        plan(multicore)
        res <- future_map2(tx_list, shapes_list, function(tx_g, shapes_g) {
            st_assign_pt(tx_g, shapes_g, shape_id_col) %>% 
                dplyr::select(ORDER, cell)
        }, .options = furrr_options(seed = TRUE))
    } else {
        res <- map2(tx_list, shapes_list, function(tx_g, shapes_g) {
            st_assign_pt(tx_g, shapes_g, shape_id_col) %>% 
                dplyr::select(ORDER, cell)
        }) 
    } 
    if (verbose) message('(5) put everything back together')
    ## NOTE: because this is a data.table, all operations below are in memory! 
    tx_dt[data.table(bind_rows(res)), on = 'ORDER', cell := i.cell]
    setorder(tx_dt, ORDER)
    tx_dt[, `:=`(ORDER = NULL, INBBOX = NULL)]
    setnafill(tx_dt, fill = 0L, cols = 'cell')
    return(tx_dt)
}


#' @export 
tx_to_counts <- function(genes, cells, remove_bg = TRUE) {
    if (remove_bg) {
        idx <- which(cells != 0)
        cells <- cells[idx]
        genes <- genes[idx]
    }
    genes <- factor(genes)
    cells <- factor(cells)
    counts <- Matrix::sparseMatrix(
        i = as.integer(genes), 
        j = as.integer(cells), 
        x = rep(1, length(genes)),
        dims = c(length(levels(genes)), length(levels(cells)))
    )
    rownames(counts) <- levels(genes)
    colnames(counts) <- levels(cells)
    return(counts)
}

#' @export 
st_aggregate_pts_shapes <- function(
    tx_dt, shapes_sf, colname_x, colname_y, colname_ptname, colname_shapename, gridn, verbose=FALSE,
    return_type=c('mat', 'list')[1]
) {
    ## Note to self: this is huge! Chop it up! 
    ## Maybe also a simpler function to just to the whole thing without tiling? 
    
    ## (1) crop to shared area 
    ##       don't care about removing points here 
    ##       
    bbox_tx <- st_rectangle(
            min(tx_dt[[colname_x]]), max(tx_dt[[colname_x]]), 
            min(tx_dt[[colname_y]]), max(tx_dt[[colname_y]])
        )
    bbox_cells <- st_as_sfc(st_bbox(shapes_sf))
    bbox <- st_intersection(bbox_tx, bbox_cells)
    bbox <- st_buffer(bbox, 1e-8) ## needed? 
    bbox <- st_as_sfc(st_bbox(bbox)) ## simplify bbox 
    bbox_coords <- st_bbox(bbox)

    if (st_is_empty(bbox)) {
        ## do something? 
        stop('No overlap between points and shapes')
    }

    tx_dt <- tx_dt[
        eval(parse(text = colname_x)) >= bbox_coords[['xmin']] & 
        eval(parse(text = colname_x)) < bbox_coords[['xmax']] & 
        eval(parse(text = colname_y)) >= bbox_coords[['ymin']] & 
        eval(parse(text = colname_y)) < bbox_coords[['ymax']], 
    ] 

    suppressWarnings({
        shapes_sf <- shapes_sf[st_intersects(bbox, shapes_sf)[[1]], ] 
        grid <- st_make_grid(bbox, n = gridn) 
    })

    ## (2) Split data 
    if (verbose) message('(3.1) Split transcripts')
    tx_list <- map(grid, function(g) {
        bbox_g <- st_bbox(g)
        tx_dt[
            eval(parse(text = colname_x)) >= bbox_g[['xmin']] & 
            eval(parse(text = colname_x)) < bbox_g[['xmax']] & 
            eval(parse(text = colname_y)) >= bbox_g[['ymin']] & 
            eval(parse(text = colname_y)) < bbox_g[['ymax']]  
        ]
    })
    grid_use <- which(map_int(tx_list, nrow) > 0)
    grid <- grid[grid_use]
    tx_list <- tx_list[grid_use]    
    tx_list <- map(tx_list, st_as_sf, coords = c(colname_x, colname_y))        
    if (verbose) message('(3.2) Split cells')    
    tile_cells <- st_intersects(grid, shapes_sf)
    shapes_list <- map(tile_cells, function(i) shapes_sf[i, ]) 


    ## (3) do point counting for each shape in each tile 
    if (verbose) message('(4) count transcripts')
    res_list <- map2(tx_list, shapes_list, function(.tx_dt, .shapes_sf) {
        pt_names <- .tx_dt[[colname_ptname]]
        overlap_res <- st_contains(
            st_geometry(.shapes_sf), ## shapes sfc 
            st_geometry(st_as_sf(.tx_dt, coords = c(colname_x, colname_y))) ## points sfc
        )
        overlap_res <- map(overlap_res, function(i) pt_names[i]) ## faster in C++? 
        data.table(
            shape_id = as.character(.shapes_sf[[colname_shapename]]), ## when integer, code below fails 
            pt_names = overlap_res,
            key = 'shape_id'
        )
    })

    ## (4) combine over tiles, taking into account shapes that appear in two tiles 
    if (verbose) message('(5) combine tiles')
    foo <- function(x, y) {
        cells_common <- intersect(x$shape_id, y$shape_id)
        if (length(cells_common) == 0) {
            res <- bind_rows(x, y)
        } else {
            res <- x[cells_common][
                y[cells_common], on = 'shape_id' ## join cells
            ][
                , .(pt_names = list(c(pt_names[[1]], i.pt_names[[1]]))), by = shape_id ## merge pts 
            ][] %>% 
                rbind(x[setdiff(x$shape_id, cells_common)]) %>% 
                rbind(y[setdiff(y$shape_id, cells_common)])
       }
        data.table::setkey(res, 'shape_id')
        return(res)    
    }
    res <- Reduce(foo, res_list)

    if (return_type == 'list') {
        return(res)
    } else if (return_type == 'mat') {
        ## (5) Collapse into matrix
        if (verbose) message('(6) collapse into matrix')
        return(tx_to_counts2(res, 'pt_names', 'shape_id'))
    } else {
        stop('invalid return_type')
    }
}
                     

#' @export 
## this version uses ipx format for dgCMatrix instead of ijx 
tx_to_counts2 <- function(cell_lists, colname_ptname, colname_shapename) {
    pt_names <- unlist(cell_lists[[colname_ptname]])
    pt_names <- factor(pt_names)
    shape_names <- cell_lists[[colname_shapename]]
    counts <- Matrix::sparseMatrix(
        i = as.integer(pt_names), ## genes 
        p = cumsum(c(0, map_int(cell_lists[[colname_ptname]], length))), ## ntx per cell 
        x = rep(1, length(pt_names)),
        dims = c(nlevels(pt_names), length(shape_names)),
        dimnames = list(levels(pt_names), shape_names)
    )    
    return(counts)
}
