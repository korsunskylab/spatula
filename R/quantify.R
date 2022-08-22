## Functions for quantifying things about cells 
## Analogous to regionprops in skimage 

## NOTES: 
##   (1) assumes column names: global_x, global_y, gene
#' @export 
assign_tx_to_cell <- function(tx, cells, micron_to_pixel_mat, do_grid, boundary=NULL, verbose=FALSE) {
    if (do_grid == TRUE & !is.null(boundary)) {
        stop('TODO: when boundary not defined, computing bounding box of cells')
    }

    if (verbose) message('Transform microns to pixel coordinates')
    ## To directly compare image coordinates to transcript coordinates, need to transform to same space 
    ## Mosaic image transformation matrix provided by Vizgen 
    # rot_mat <- fread('micron_to_mosaic_pixel_transform.csv')
    # rot_mat <- fread('/mnt/efs/fs1/vizgen/vizgen_analyses/HuColonCa_FFPE_PH1_CellBoundary_QC_V8_LH_10-13-2021/region_0/micron_to_mosaic_pixel_transform.csv')
    rot_mat <- t(as.matrix(micron_to_pixel_mat))

    ## NOTE: the rotation matrix assumes that third column is intercept (not z!)
    tx_coords_transformed <- as.matrix(cbind(tx[, .(global_x, global_y)], 1)) %*% rot_mat
    tx_coords_transformed <- data.frame(tx_coords_transformed[, 1:2])
    colnames(tx_coords_transformed) <- c('X1', 'X2')


    ## Sanity check: these must be zero. Otherwise, transformation is wrong 
    stopifnot(all(tx_coords_transformed$X1 >= 0) & all(tx_coords_transformed$X2 >= 0))

    tx <- tx %>% 
        cbind(data.frame(tx_coords_transformed)) %>% 
        dplyr::select(x=X1, y=X2, gene) %>% 
        tibble::rowid_to_column('rowID')
    head(tx)


    if (verbose) message('Filter cells to be within transcript bbox')
    ## coords: xmin, xmax, ymin, ymax
    bbox_tx = st_rectangle(c(min(tx$x), max(tx$x), min(tx$y), max(tx$y)))
    cells_sf <- st_sf(
        ID = paste0('Cell', 1:length(cells)),
        cell = cells
    )
    cell_names <- cells_sf %>% st_crop(bbox_tx) %>% with(ID) 
    cells_sf <- subset(cells_sf, ID %in% cell_names)
    nrow(cells_sf)

    if (verbose) message('Make grid')
    if (do_grid) {
        boundary <- st_crop(st_sfc(boundary), bbox_tx)[[1]]
        grid <- st_intersection(st_make_grid(boundary, n = 10), boundary)
    } else {
        ## if no grid, make grid one point 
        ## NOTE: st_bbox was not working here, so here is a manual version: 
        bbox_cells <- st_rectangle(c(
            min(st_coordinates(cells_sf$cell)[, 'X']), 
            max(st_coordinates(cells_sf$cell)[, 'X']), 
            min(st_coordinates(cells_sf$cell)[, 'Y']), 
            max(st_coordinates(cells_sf$cell)[, 'Y'])
        ))
        grid <- st_sfc(bbox_cells)
    }

    if (verbose) message('Aggregate transcripts within grid sections')
    ## TODO: move processing inside one grid to its own function to not repeat code 
    if (do_grid) {
        ## Need to remember cell names to deal with cells on the border 
        ## NOTE: most of the time is spent subsetting cells with st_crop
        ##       if we knew that cells were sorted spatially, we could compute bbox for those cells 
        ## CAUTION: future_imap passes full copy of tx table to each thread, 
        ##          so this is very memory intensive
        ##          let's consider explicitely splitting tx first by grid location and then processing each piece of parallel 
        parallel = FALSE
        if (parallel) {
            options(future.globals.maxSize=1e10)
            plan(multicore)
            iter_fxn <- function(...) {
                furrr::future_imap(..., .options = furrr_options(seed = TRUE))
            }
        } else {
            iter_fxn <- purrr::imap
        }

        tx_segmented <- iter_fxn(grid, function(.grid, i) {
            message(i)
            ## Subset data to this part of grid 
            .tx <- tx[
                between(x, st_bbox(.grid)[['xmin']], st_bbox(.grid)[['xmax']]) & 
                between(y, st_bbox(.grid)[['ymin']], st_bbox(.grid)[['ymax']])
            ]
            if (nrow(.tx) == 0) return(.tx)

            ## get IDs of cells to avoiding cropping cell at border 
            suppressWarnings({
                cell_names <- cells_sf %>% st_crop(.grid) %>% with(ID)        
                if (length(cell_names) == 0) return(.tx)
                .cells <- subset(cells_sf, ID %in% cell_names)$cell

                ## intersect the points and polygons
                tx_points <- .tx %>% 
                    dplyr::select(x, y) %>% 
                    as.matrix() %>% 
                    st_multipoint() %>% 
                    st_sfc() %>% 
                    st_cast('POINT') 
                overlap_res <- st_contains(.cells, tx_points, sparse = TRUE)
            })

            ## assign each transcript to a cell (if inside)
            .tx$cell <- st_intersects(tx_points, .cells) %>% 
                map(head, 1) %>% ## avoid double assignment
                as.integer()

            ## by convention, set background to cell 0 
            .tx$cell[is.na(.tx$cell)] <- 0 
            .tx$cell <- c('extracellular', cell_names)[.tx$cell + 1]        
            return(.tx)
        }) %>% 
            bind_rows() %>% 
            dplyr::arrange(rowID) ## return in same order as original tx     
        return(tx_segmented)
    } else {
        tx_points <- tx %>% 
            dplyr::select(x, y) %>% 
            as.matrix() %>% 
            st_multipoint() %>% 
            st_sfc() %>% 
            st_cast('POINT') 
        overlap_res <- st_contains(cells_sf$cell, tx_points, sparse = TRUE)
        tx$cell <- st_intersects(tx_points, cells_sf$cell) %>% 
            map(head, 1) %>% ## avoid double assignment
            as.integer()

        ## by convention, set background to cell 0 
        tx$cell[is.na(tx$cell)] <- 0 
        tx$cell <- c('extracellular', cell_names)[tx$cell + 1]        
        return(tx)
    }
    
}
    