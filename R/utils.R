#' @export 
dgCMat_to_pairs <- function(mat, weighted=FALSE) {
    res <- dgCMat_to_pairs_cpp(mat@i, tail(mat@p, -1), mat@x)
    if (weighted) {
        return(res)
    } else {
        return(res[, 1:2])
    }
}
