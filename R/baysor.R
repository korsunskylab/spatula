## Utility functions for running Baysor with tiling, from R 

#' @export 
baysor.split_tx_files <- function(output_dir, max_tx_per_grid) {
    ## Check if dir exists
    if (!dir.exists(output_dir)) {
        dir.create(output_dir)
    }

    ## Check if dir is writeable 
    if (file.access(output_dir, mode = 2) == -1L) {
        stop('Specified output_dir does not have write permissions.')
    }

    grid <- split_tx(
        tx, 
        max_tx = max_tx_per_grid, 
        max_voxels = round(nrow(tx) / max_tx_per_grid) * 5 ## to prevent worst case behavior of too much splitting for 
    ) 


    err_status <- grid %>%
        map('tx') %>% 
        imap(function(tx_grid, grid_id) {
            dirname <- file.path(output_dir, paste0('g', grid_id))        
            if (!dir.exists(dirname)) dir.create(dirname)
            fname <- as.character(glue('{dirname}/tx_baysor.csv'))
            fwrite(tx_grid, fname, sep = ',') 
        })    
    return(length(grid))
}


#' @export 
baysor.collect_tx <- function(dir, testing=FALSE, mintx=10) {
    ## Read cell assignment from each tile 
    # plan(multicore)
    tx <- list.files(dir, pattern = 'segmentation.csv$', recursive = TRUE, full.names = TRUE) %>% 
        map(fread) %>% 
        # future_map(fread) %>% 
        data.table::rbindlist(idcol = 'tile')

    if (testing) {
        ncells <- nrow(unique(tx[cell != 0, .(tile, cell)]))
    }

    ## renumber cells to be consecutive integers 
    tx <- unique(tx[cell != 0, .(tile, cell)])[
        , .N, by = tile
    ][
        , offset := cumsum(N)
    ][
        , offset := c(0, offset[1:(.N-1)])
    ][
        tx, on = 'tile'
    ][
        cell != 0, cell := cell + offset    
    ][
        , `:=`(tile = NULL, N = NULL, offset = NULL) ## remove this house-keeping information 
    ][] 

    if (testing) {
        stopifnot(ncells == max(tx$cell))
    }    
    
    ## Release transcripts from low-count cells back to the wild 
    tx[
        , N := sum(!grepl('Blank', gene)), by = cell 
    ][
        , cell := case_when(
            N >= mintx ~ cell, 
            TRUE ~ 0L
        )
    ]
    return(tx)
}

#' @export 
baysor.collect_cells <- function(tx, min_tx_per_cell) {
    cells <- infer_polygons(tx, 'x', 'y', min_tx_per_cell, parallel = TRUE)
    cells <- cbind(cells, st_coordinates(st_centroid(cells$shape))) %>% 
        dplyr::rename(x = X, y = Y)
    cells$area <- st_area(cells$shape)
    cells <- cells %>% 
        left_join(
            tx[
                cell != 0, 
                .(
                    n_transcripts = .N, 
                    avg_confidence = mean(assignment_confidence), 
                    cluster = .SD[, .N, by = cluster][order(-N)][1, cluster]
                ), 
                by = cell
            ][
                , cell := as.character(cell)
            ], 
            by = c('id' = 'cell')
        )
    ## Tableau classic 20
    pal <- c(
        '#1f77b4', '#aec7e8', '#ff7f0e', '#ffbb78', '#2ca02c', '#98df8a', 
        '#d62728', '#ff9896', '#9467bd', '#c5b0d5', '#8c564b', '#c49c94', '#e377c2', '#f7b6d2', 
        '#7f7f7f', '#c7c7c7', '#bcbd22', '#dbdb8d', '#17becf', '#9edae5'
    )
    cells$rand_color <- sample(pal, nrow(cells), TRUE)
    return(cells)
}



#' @export 
baysor.prepare_cmds <- function(ntiles, output_dir, fname_baysor_config, max_parallel_jobs) {
    output_dir <- normalizePath(output_dir)
    # dir <- normalizePath(params_pipeline$output_dir)
    file.copy(
        from = fname_baysor_config, 
        to = file.path(output_dir, 'config.toml')
    )

    ## Create list of baysor commands to run 
    ## Write them to file - scheduler will read from here 
    i <- seq_len(ntiles)
    baysor_cmds <- glue('nohup baysor run -c datadir/config.toml -o datadir/g{i} datadir/g{i}/tx_baysor.csv :cell 1> datadir/g{i}/out.out 2> datadir/g{i}/err.err &') 
    writeLines(baysor_cmds, file.path(output_dir, 'Baysor_cmd_list.txt'))

    ## Write scheduler script 
    ## Because this is run inside Docker, directory is simply /root/datadir
    scheduler_code <- paste(c(
    '#!/bin/bash -e
    readarray -t cmds < datadir/Baysor_cmd_list.txt
    declare -i maxRunningJobs=',
    # params_pipeline$max_parallel_jobs,
    max_parallel_jobs, 
    '\ndeclare -i totalJobs=${#cmds[@]}
    declare -i jobsSent=0
    while [ $jobsSent -lt $totalJobs ]
    do
        if [ `jobs -r | wc -l | tr -d " "` -lt $maxRunningJobs ]
        then 
            echo "${cmds[jobsSent]}"
            eval "${cmds[jobsSent]}"
            jobsSent=$jobsSent+1
        fi 
        sleep .1
    done
    wait
    '
    ), collapse = '')
    writeLines(scheduler_code, file.path(output_dir, 'scheduler.sh'))
    cmd <- as.character(glue('chmod a+x {output_dir}/scheduler.sh'))
    system(cmd, wait = TRUE)
}

#' @export 
baysor.run_baysor_tiles_docker <- function(ntiles, output_dir, fname_baysor_config) {
    output_dir <- normalizePath(output_dir)
    file.copy(
        from = fname_baysor_config, 
        to = file.path(output_dir, 'config.toml')
    )

    ## Open one instance of Docker
    cmd <- as.character(glue('docker run --name MY_INSTANCE_NAME -t -d -v {output_dir}:/root/datadir vpetukhov/baysor:master'))
    system(cmd, wait = TRUE)

    ## Run scheduler inside Docker 
    cmd <- 'docker exec MY_INSTANCE_NAME datadir/scheduler.sh'
    system(cmd, wait = TRUE)

    ## Once schedule is finished, stop and kill the docker instance 
    system('docker stop MY_INSTANCE_NAME', wait = TRUE)
    system('docker rm MY_INSTANCE_NAME', wait = TRUE)    
}

# baysor.run_baysor_tiles_binary <- function(ntiles, output_dir, fname_baysor_config) {

# }


#' @export 
baysor.finish <- function() {
    ## Remove temporary files 
    if (opts$remove_temp_files) {
        for (fname in list.files(output_dir, full.names = TRUE)) {
            unlink(fname, recursive = TRUE)
        }    
    }
    

    ## Write out all main output files 
    fwrite(tx, file.path(output_dir, 'transcripts.csv'))
    writeLines(rownames(counts), file.path(output_dir, 'genes.tsv'))
    spatula::writeMM(counts, file.path(output_dir, 'counts.mtx'))
    ## parquet files are readable in R and python 
    st_write_parquet(dplyr::select(cells, id, shape), file.path(output_dir, 'shapes.parquet'))    
    fwrite(st_drop_geometry(cells), file.path(output_dir, 'cells.csv'), sep = ',')
    
}

#' @export 
baysor.run <- function(
    tx, 
    fname_baysor_config, 
    output_dir, 
    opts = list(
        baysor_mode = c('docker', 'binary')[1], ## baysor can be run through docker container or natively installed binary 
        max_tx_per_grid = 2e7, ## maximum number of transcripts per tile 
        compute_shapes = TRUE, 
        max_parallel_jobs = 10, ## maximum number of concurrent Baysor jobs running 
        min_tx_per_cell = 10, 
        remove_temp_files = TRUE ## set to FALSE for debugging 
    )
) {
    ## Check inputs 
    message('WARNING: for now, tx needs to have columns: x, y, gene, cell')
    stopifnot(colnames(tx) == c('x', 'y', 'gene', 'cell'))
    stopifnot(file.exists(fname_baysor_config))
    if (opts$baysor_mode == 'docker') {
        if (!'docker' %in% strsplit(system('groups', intern = TRUE), ' ')[[1]]) {
            stop('You do not belong to the docker group. Please follow the instructions here to run docker without root privileges: https://docs.docker.com/engine/install/linux-postinstall/')
        }
    }
    
    ntiles <- baysor.split_tx_files(output_dir, opts$max_tx_per_grid) 
    baysor.prepare_cmds(ntiles, output_dir, fname_baysor_config, opts$max_parallel_jobs)
    if (opts$baysor_mode == 'docker') {
        baysor.run_baysor_tiles_docker(ntiles, output_dir, fname_baysor_config)    
    } else if (opts$baysor_mode == 'binary') {
        stop('binary mode not implemented yet. Please use docker mode.')
        # baysor.run_baysor_tiles_binary(ntiles, output_dir, ..., fname_baysor_config) 
    } else {
        stop('invalid baysor_mode')
    }
    
    tx <- baysor.collect_tx(output_dir, testing=FALSE, mintx=opts$min_tx_per_cell) 
    counts <- tx_to_counts(tx$gene, tx$cell, remove_bg = TRUE)
    cells <- baysor.collect_cells(tx, opts$min_tx_per_cell)
    environment(baysor.finish) <- environment() 
    baysor.finish() 
}
