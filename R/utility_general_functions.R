#' @useDynLib bambu, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

#' @noRd
.onUnload <- function(libpath) {
    library.dynam.unload("bambu", libpath)
}
