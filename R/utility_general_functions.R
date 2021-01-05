#' @useDynLib bambu, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

#' @noRd
.onUnload <- function(libpath) {
    library.dynam.unload("bambu", libpath)
}


#' @noRd
helpFun <- function(chr, chrRanges, bamFile) {
    return(GenomicAlignments::grglist(GenomicAlignments::readGAlignments(
        file = bamFile,
        param = Rsamtools::ScanBamParam(
            flag = Rsamtools::scanBamFlag(isSecondaryAlignment = FALSE),
            which = chrRanges[chr]),
        use.names = FALSE)))
}


