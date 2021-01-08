#' Function to prepare reads for processing from bam file
#' @param bamFile bamFile
#' @inheritParams bambu
#' @noRd
prepareDataFromBam <- function(bamFile, yieldSize = NULL, verbose = FALSE, 
                               seqlevels = NULL) {
    if (methods::is(bamFile, "BamFile")) {
        if (!is.null(yieldSize)) {
            Rsamtools::yieldSize(bamFile) <- yieldSize
        } else {
            yieldSize <- Rsamtools::yieldSize(bamFile)
        }
    } else if (!grepl(".bam", bamFile)) {
        stop("Bam file is missing from arguments.")
    } else {
        if (is.null(yieldSize)) {
            yieldSize <- NA
        }
        bamFile <- Rsamtools::BamFile(bamFile, yieldSize = yieldSize)
    }
    bf <- open(bamFile)
    readGrgList <- list()
    counter <- 1
    while (Rsamtools::isIncomplete(bf)) {
        readGrgList[[counter]] <-
            GenomicAlignments::grglist(GenomicAlignments::readGAlignments(bf,
            param = Rsamtools::ScanBamParam(flag =
                Rsamtools::scanBamFlag(isSecondaryAlignment = FALSE)),
            use.names = FALSE))
        # readGrgList<-c(readGrgList,GenomicAlignments::grglist(reads))
        if (verbose) show(min(length(readGrgList),
            counter * yieldSize, na.rm = TRUE))
        counter <- counter + 1
    }
    on.exit(close(bf))
    if (length(readGrgList) > 1) {
        readGrgList <- do.call(c, readGrgList)
    } else {
        readGrgList <- readGrgList[[1]]
    }
    # remove microexons of width 1bp from list
    readGrgList <- readGrgList[GenomicRanges::width(readGrgList) > 1]
    if(!is.null(seqlevels)) { # only keep ranges from seqlevels provided
        readGrgList <- GenomeInfoDb::keepSeqlevels(readGrgList,
                                          value = seqlevels,
                                          pruning.mode = "coarse")
    }
    mcols(readGrgList)$id <- seq_along(readGrgList)
    return(readGrgList)
}

