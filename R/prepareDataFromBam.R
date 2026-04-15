#' Function to prepare reads for processing from bam file
#' @param bamFile bamFile
#' @inheritParams bambu
#' @importFrom methods is
#' @importFrom Rsamtools yieldSize yieldSize<- BamFile isIncomplete 
#'     ScanBamParam scanBamFlag
#' @importFrom GenomicAlignments grglist readGAlignments
#' @importFrom GenomicRanges width
#' @noRd
prepareDataFromBam <- function(bamFile, yieldSize = NULL, verbose = FALSE,
        use.names = FALSE, extractBarcodeUMI = FALSE, dedupUMI = FALSE) {
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
    use.names.OG <- use.names
    if(extractBarcodeUMI) use.names <- TRUE
    while (isIncomplete(bf)) {
        alignmentInfo <- readGAlignments(bf, param = ScanBamParam(tag = c("CB", "UB"), 
                                         flag = scanBamFlag(isSecondaryAlignment = FALSE)), 
                                         use.names = use.names)
        readGrgList[[counter]] <-grglist(alignmentInfo)
        if(extractBarcodeUMI){
            # parse CB and UMI from the bam file, either from CB/UB tags or read names.
            # read name format: CB_UMI#READNAME (CB & UMI cannot have '_', otherwise parsing fails)
            mcols(readGrgList[[counter]])$CB <- case_when(
                !is.na(mcols(alignmentInfo)$CB) ~ mcols(alignmentInfo)$CB,
                grepl("^[^_]+_[^#]+#", names(readGrgList[[counter]]), perl = TRUE) ~ sub("_.*", "", names(readGrgList[[counter]])),
                TRUE ~ NA
            )

            mcols(readGrgList[[counter]])$UMI <- case_when(
                !is.na(mcols(alignmentInfo)$UB) ~ mcols(alignmentInfo)$UB,
                grepl("^[^_]+_[^#]+#", names(readGrgList[[counter]]), perl = TRUE) ~ sub("^[^_]+_([^#]+)#.*$", "\\1", names(readGrgList[[counter]])),
                TRUE ~ NA
            )
        }
        
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

    if (extractBarcodeUMI){
        mcols(readGrgList)$CB <- factor(mcols(readGrgList)$CB, levels = sort(unique(mcols(readGrgList)$CB)))
    }

    # remove microexons of width 1bp from list
    readGrgList <- readGrgList <- readGrgList[sum(width(readGrgList)) > 1]
    numNoCBs <- sum(is.na(mcols(readGrgList)$CB))
    if(numNoCBs > 0){
        message("Removing ", numNoCBs, " reads that were not assigned barcodes. If this is unexpected check the barcode map input")
        readGrgList <- readGrgList[!is.na(mcols(readGrgList)$CB)]
    }
    if(dedupUMI){
        #UMI deduplication by barcode
        start.ptm <- proc.time()
        numUMIs <- length(na.omit(unique(mcols(readGrgList)$UMI))) # remove NA UMIs
        if(numUMIs > 100){
            df <- data.frame(umi = mcols(readGrgList)$UMI, 
                barcode = mcols(readGrgList)$CB,
                lengths = sum(width(readGrgList)))
            df <- df %>% mutate(id = row_number()) %>% group_by(barcode, umi) %>% summarise(primary.id = id[which.max(lengths)])
            readGrgList <- readGrgList[df$primary.id]
            end.ptm <- proc.time()
            message("UMI deduplication time ", round((end.ptm - start.ptm)[3] / 60, 3), " mins.")
        } else {
            message("Only ", numUMIs, " detected. Not performing UMI deduplication. If this is unexpected, double check the --chemistry argument")
        }

 
        #readGrgList = c(readGrgList.filt, unname(readGrgList.keep))
    }
    if(!use.names.OG) names(readGrgList) <- NULL
    seqlevels(readGrgList) <- as.character(unique(getChrFromGrList(readGrgList)))
    return(readGrgList)
}


#' Function to clip sequences 
#' @noRd
clipFunction <- function(cigarData, grep_pattern, replace_pattern){
    return(suppressWarnings(pmax(0,as.numeric(gsub(grep_pattern,replace_pattern,
                                                   cigarData)), na.rm=T)))
}
