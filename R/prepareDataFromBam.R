#' Function to prepare reads for processing from bam file
#' @param bamFile bamFile
#' @inheritParams bambu
#' @importFrom methods is
#' @importFrom Rsamtools yieldSize yieldSize<- BamFile isIncomplete 
#'     ScanBamParam scanBamFlag
#' @importFrom GenomicAlignments grglist readGAlignments
#' @importFrom GenomicRanges width
#' @noRd
prepareDataFromBam <- function(bamFile, yieldSize = NULL, verbose = FALSE, use.names = FALSE, demultiplexed = NULL) {
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
    cells <- c()
    umi <- c()
    
    while (isIncomplete(bf)) {
        readGrgList[[counter]] <-
            grglist(readGAlignments(bf,
            param = ScanBamParam(flag =
                scanBamFlag(isSecondaryAlignment = FALSE)),
            use.names = use.names))
        
        ### add ### 
        if (isTRUE(demultiplexed)){
          mcols(readGrgList[[counter]])$CB <-  substr(names(readGrgList[[counter]]), 1, 16)
          mcols(readGrgList[[counter]])$UMI <- substr(names(readGrgList[[counter]]), 18, 29)
          names(readGrgList[[counter]]) <- NULL 
          cells <- unique(c(cells, mcols(readGrgList[[counter]])$CB))
          mcols(readGrgList[[counter]])$CB <- factor(mcols(readGrgList[[counter]])$CB, levels = cells)
          
          umi <- unique(c(umi, mcols(readGrgList[[counter]])$UMI))
          mcols(readGrgList[[counter]])$UMI <- factor(mcols(readGrgList[[counter]])$UMI, levels = umi)
        }
        ### add ### 
        
        counter <- counter + 1
    }
    on.exit(close(bf))
    rm(cells)
    rm(umi)
    if (length(readGrgList) > 1) {
        readGrgList <- do.call(c, readGrgList)
    } else {
        readGrgList <- readGrgList[[1]]
    }
    # remove microexons of width 1bp from list
    readGrgList <- readGrgList[width(readGrgList) > 1]
    
    return(readGrgList)
}

