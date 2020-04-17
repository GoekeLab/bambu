#' @useDynLib bamboo, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

#' @noRd
.onUnload <- function (libpath) { library.dynam.unload("bamboo", libpath) }

