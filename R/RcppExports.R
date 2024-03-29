# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

triplets_to_pairs <- function(triplets) {
    .Call('_spatula_triplets_to_pairs', PACKAGE = 'spatula', triplets)
}

dgCMat_to_pairs_cpp <- function(ivec, pvec, xvec) {
    .Call('_spatula_dgCMat_to_pairs_cpp', PACKAGE = 'spatula', ivec, pvec, xvec)
}

get_idx_cpp <- function(X, only_boundary = TRUE) {
    .Call('_spatula_get_idx_cpp', PACKAGE = 'spatula', X, only_boundary)
}

