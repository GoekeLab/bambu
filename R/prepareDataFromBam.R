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
        use.names = FALSE, demultiplexed = FALSE, cleanReads = FALSE) {
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
    use.names.OG = use.names
    if(demultiplexed | cleanReads) use.names = TRUE
    while (isIncomplete(bf)) {
        alignmentInfo <- readGAlignments(bf, param = ScanBamParam(tag = c("BC", "UG"), 
                                         flag = scanBamFlag(isSecondaryAlignment = FALSE)), 
                                         use.names = use.names)
        readGrgList[[counter]] <-grglist(alignmentInfo)
        if (isTRUE(demultiplexed)){
            mcols(readGrgList[[counter]])$CB <- ifelse(!is.na(mcols(alignmentInfo)$BC), mcols(alignmentInfo)$BC, 
                                                       substr(names(readGrgList[[counter]]), 1, 16))
            
            mcols(readGrgList[[counter]])$UMI <- ifelse(!is.na(mcols(alignmentInfo)$UG), mcols(alignmentInfo)$UG, 
                                                       substr(names(readGrgList[[counter]]), 18, 29))
            
            cells <- unique(c(cells, mcols(readGrgList[[counter]])$CB))
            mcols(readGrgList[[counter]])$CB <- factor(mcols(readGrgList[[counter]])$CB, levels = cells)
            umi <- unique(c(umi, mcols(readGrgList[[counter]])$UMI))
            mcols(readGrgList[[counter]])$UMI <- factor(mcols(readGrgList[[counter]])$UMI, levels = umi)
        }
        if(cleanReads){
            softClip5Prime <-pmax(0,as.numeric(gsub('^(\\d*)[S].*','\\1',GenomicAlignments::cigar(alignmentInfo))), na.rm=T)
            softClip3Prime <-pmax(0,as.numeric(gsub('.*\\D(\\d*)[S]$','\\1',GenomicAlignments::cigar(alignmentInfo))), na.rm=T)
            hardClip5Prime <-pmax(0,as.numeric(gsub('^(\\d*)[H].*','\\1',GenomicAlignments::cigar(alignmentInfo))), na.rm=T)
            hardClip3Prime <-pmax(0,as.numeric(gsub('.*\\D(\\d*)[H]$','\\1',GenomicAlignments::cigar(alignmentInfo))), na.rm=T)
            mcols(readGrgList[[counter]])$clip5Prime = pmax(softClip5Prime, hardClip5Prime)
            mcols(readGrgList[[counter]])$clip3Prime = pmax(softClip3Prime, hardClip3Prime)
            rev = as.vector(strand(alignmentInfo) == '-')
            rev2 = grepl("_-.+of", names(alignmentInfo))
            temp = mcols(readGrgList[[counter]])$clip5Prime
            mcols(readGrgList[[counter]])$clip5Prime[rev != rev2] = mcols(readGrgList[[counter]])$clip3Prime[rev != rev2]
            mcols(readGrgList[[counter]])$clip3Prime[rev != rev2] = temp[rev != rev2]
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
    # remove microexons of width 1bp from list
    readGrgList <- readGrgList[width(readGrgList) > 1]
    if(cleanReads){
        #extract duplicated reads from flexiplex to clean
        #leave other reads alone as supplimental alignments maybe fusion transcripts
        #commented out because it takes awhile
        # dt = data.table(name = substr(readNames2,0,nchar(readNames2[1])-2), 
        #          strand = substr(readNames2,nchar(readNames2[1]),nchar(readNames2[1])))
        # dt[, id := .I]
        # dt <- dt[, .(ids = list(id), toFilt = any(strand == "+" & strand == "-")), by = name]
        # readGrgList.keep = readGrgList[c(dt$ids[!dt$toFilt])]
        # readGrgList.filt = readGrgList[c(dt$ids[dt$toFilt])]

        df = data.frame(name = names(readGrgList), 
            clip5 = mcols(readGrgList)$clip5Prime)
        df = df %>% mutate(id = row_number()) %>% group_by(name) %>% summarise(primary.id = id[which.min(clip5)])
        readGrgList = unname(readGrgList[df$primary.id])
        #readGrgList = c(readGrgList.filt, unname(readGrgList.keep))
    }
    if(!use.names.OG) {names(readGrgList) <- NULL }
    return(readGrgList)
}

