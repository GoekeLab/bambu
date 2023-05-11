#' Function to prepare reads for processing from bam file
#' @param bamFile bamFile
#' @inheritParams bambu
#' @importFrom methods is
#' @importFrom Rsamtools yieldSize yieldSize<- BamFile isIncomplete 
#'     ScanBamParam scanBamFlag
#' @importFrom GenomicAlignments grglist readGAlignments
#' @importFrom GenomicRanges width
#' @noRd
prepareDataFromBam <- function(bamFile, yieldSize = NULL, verbose = FALSE, use.names = FALSE) {
    if (is(bamFile, "BamFile")) {
        if (!is.null(yieldSize)) {
            yieldSize(bamFile) <- yieldSize
        } else {
            yieldSize <- yieldSize(bamFile)
        }
    } else if (!grepl(".bam", bamFile)) {
        stop("Bam file is missing from arguments.")
    } else {
        if (is.null(yieldSize)) {
            yieldSize <- NA
        }
        bamFile <- BamFile(bamFile, yieldSize = yieldSize)
    }
    bf <- open(bamFile)
    readGrgList <- list()
    counter <- 1
    while (isIncomplete(bf)) {
        reads = readGAlignments(bf,
            param = ScanBamParam(flag =
                scanBamFlag(isSecondaryAlignment = FALSE), what = "cigar"),
            use.names = use.names)
        readGrgList[[counter]] = grglist(reads)
        softClip5Prime <-pmax(0,as.numeric(gsub('^(\\d*)[S].*','\\1',mcols(reads)$cigar)), na.rm=T)
        softClip3Prime <-pmax(0,as.numeric(gsub('.*\\D(\\d*)[S]$','\\1',mcols(reads)$cigar)), na.rm=T)
        hardClip5Prime <-pmax(0,as.numeric(gsub('^(\\d*)[H].*','\\1',mcols(reads)$cigar)), na.rm=T)
        hardClip3Prime <-pmax(0,as.numeric(gsub('.*\\D(\\d*)[H]$','\\1',mcols(reads)$cigar)), na.rm=T)
        mcols(readGrgList[[counter]])$softClip5Prime = softClip5Prime
        mcols(readGrgList[[counter]])$softClip3Prime = softClip3Prime
        mcols(readGrgList[[counter]])$hardClip5Prime = hardClip5Prime
        mcols(readGrgList[[counter]])$hardClip3Prime = hardClip3Prime
        counter <- counter + 1
    }
    on.exit(close(bf))
    if (length(readGrgList) > 1) {
        readGrgList <- do.call(c, readGrgList)
    } else {
        readGrgList <- readGrgList[[1]]
    }
    # remove microexons of width 1bp from list
    readGrgList <- readGrgList[width(readGrgList) > 1]
    mcols(readGrgList)$id <- seq_along(readGrgList)
    return(readGrgList)
}

