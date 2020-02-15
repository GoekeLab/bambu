#' @useDynLib bamboo, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL


.onUnload <- function (libpath) { library.dynam.unload("bamboo", libpath) }

