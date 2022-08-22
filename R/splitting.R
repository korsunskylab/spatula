#' @export 
split_grid <- function(grid_point) {
    bbox <- grid_point$bbox
    xmid <- mean(bbox[1:2])
    ymid <- mean(bbox[3:4])
    ## NOTE: the 1e-10 term avoids duplicates 
    bboxes_new <- list(
        c(bbox[1], xmid - 1e-10, bbox[3], ymid - 1e-10),
        c(xmid, bbox[2], bbox[3], ymid - 1e-10),
        c(xmid, bbox[2], ymid, bbox[4]),
        c(bbox[1], xmid - 1e-10, ymid, bbox[4])
    )
    
    ## Make four new grid points 
    res <- map(bboxes_new, function(bbox_test) {
        res <- list(
            tx = grid_point$tx[between(global_x, bbox_test[1], bbox_test[2]) & between(global_y, bbox_test[3], bbox_test[4])], 
            bbox = bbox_test,
            bbox_geom = st_rectangle(bbox_test)
        )
        res$n <- nrow(res$tx)
        return(res)
    })
}


#' @export 
split_tx <- function(tx, max_tx, max_voxels) {
    ## Initialize grid with all transcripts
    grid <- list(
        list(
            tx = tx, 
            bbox = c(min(tx$global_x), max(tx$global_x), min(tx$global_y), max(tx$global_y))
        )
    )
    grid[[1]]$n <- nrow(grid[[1]]$tx)
    grid[[1]]$bbox_geom <- st_rectangle(grid[[1]]$bbox)

    ## Keep splitting grid points until each has at most max_tx transcripts
    .i <- 0
    while (TRUE) {
        .i <- .i + 1
        grids_split <- which(map_int(grid, 'n') > max_tx)
        if (length(grid) >= max_voxels) {
            break
        } else if (length(grids_split) > 0) {
            i <- grids_split[1]
            grid <- append(grid, split_grid(grid[[i]]))
            grid[[i]] <- NULL
            grids_split <- which(map_int(grid, 'n') > max_tx)
        } else {
            break
        }
    }
    ## For QC purposes, compute the transcript density of each region 
    for (i in seq_len(length(grid))) {
        grid[[i]]$density <- grid[[i]]$n / st_area(grid[[i]]$'bbox_geom')
    }    
    return(grid)
}
