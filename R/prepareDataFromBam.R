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
    if(demultiplexed) use.names = TRUE
    while (isIncomplete(bf)) {
        ### add ###
        alignmentInfo <- readGAlignments(bf, param = ScanBamParam(tag = c("BC", "UG"), 
                                         flag = scanBamFlag(isSecondaryAlignment = FALSE)), 
                                         use.names = use.names)
        ### add ###
        readGrgList[[counter]] <-grglist(alignmentInfo)
        ### add ### 
        if (isTRUE(demultiplexed)){
            mcols(readGrgList[[counter]])$CB <- ifelse(!is.na(mcols(alignmentInfo)$BC), mcols(alignmentInfo)$BC, 
                                                       substr(names(readGrgList[[counter]]), 1, 16))
            
            mcols(readGrgList[[counter]])$UMI <- ifelse(!is.na(mcols(alignmentInfo)$UG), mcols(alignmentInfo)$UG, 
                                                       substr(names(readGrgList[[counter]]), 18, 29))
            
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

