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
        use.names = FALSE, demultiplexed = FALSE, cleanReads = TRUE, dedupUMI = FALSE) {
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
    if(!isFALSE(demultiplexed) | cleanReads) use.names = TRUE
    if(grepl(".[ct]sv$",demultiplexed)){
        readMap = read.table(demultiplexed, 
            sep = ifelse(grepl(".tsv$",demultiplexed), "\t", ","), header = FALSE)
    }
    while (isIncomplete(bf)) {
        alignmentInfo <- readGAlignments(bf, param = ScanBamParam(tag = c("BC", "UG"), 
                                         flag = scanBamFlag(isSecondaryAlignment = FALSE)), 
                                         use.names = use.names)
        readGrgList[[counter]] <-grglist(alignmentInfo)
        if (!isFALSE(demultiplexed)){
            if(isTRUE(demultiplexed)){
                mcols(readGrgList[[counter]])$BC <- ifelse(!is.na(mcols(alignmentInfo)$BC), mcols(alignmentInfo)$BC, 
                                                        gsub("(^[GACT]+(?=_)).*", '\\1', names(readGrgList[[counter]]), perl = TRUE))
                mcols(readGrgList[[counter]])$UMI <- ifelse(!is.na(mcols(alignmentInfo)$UG), mcols(alignmentInfo)$UG, 
                                                        gsub(".*((?<=_)[GACT]*(?=#)).*", '\\1', names(readGrgList[[counter]]), perl = TRUE))
            } else{
                mcols(readGrgList[[counter]])$BC = NA
                mcols(readGrgList[[counter]])$UMI = "NA"
                mcols(readGrgList[[counter]])$BC = readMap[,2][match(names(readGrgList[[counter]]),readMap[,1])]
                if(ncol(readMap)>2){
                    mcols(readGrgList[[counter]])$UMI = readMap[,3][match(names(readGrgList[[counter]]),readMap[,1])]
                }
            }
            cells <- unique(c(cells, mcols(readGrgList[[counter]])$BC))
            mcols(readGrgList[[counter]])$BC <- factor(mcols(readGrgList[[counter]])$BC, levels = cells)
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
    numNoBCs = sum(is.na(mcols(readGrgList)$BC))
    if(numNoBCs > 0){
        message("Removing ", numNoBCs, " reads that were not assigned barcodes. If this is unexpected check the barcode map input")
    readGrgList = readGrgList[!is.na(mcols(readGrgList)$BC)]
    }
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

       #select alignments closest to barcode
        start.ptm <- proc.time()
        df = data.frame(name = names(readGrgList), 
            clip5 = mcols(readGrgList)$clip5Prime)
        df = df %>% mutate(id = row_number()) %>% group_by(name) %>% summarise(primary.id = id[which.min(clip5)])
        readGrgList = unname(readGrgList[df$primary.id])
        end.ptm <- proc.time()
        message("Primary alignment selection time ", round((end.ptm - start.ptm)[3] / 60, 3), " mins.")
        
    }
    if(dedupUMI){
        #UMI deduplication by barcode
        start.ptm <- proc.time()
        numUMIs = length(unique(mcols(readGrgList)$UMI))
        if(numUMIs > 100){
            df = data.frame(umi = mcols(readGrgList)$BC, 
                barcode = mcols(readGrgList)$UMI,
                lengths = sum(width(readGrgList)))
            df = df %>% mutate(id = row_number()) %>% group_by(barcode, umi) %>% summarise(primary.id = id[which.max(lengths)])
            readGrgList = readGrgList[df$primary.id]
            end.ptm <- proc.time()
            message("UMI deduplication time ", round((end.ptm - start.ptm)[3] / 60, 3), " mins.")
        } else {
            message("Only ", numUMIs, " detected. Not performing UMI deduplication. If this is unexpected, double check the --chemistry argument")
        }

 
        #readGrgList = c(readGrgList.filt, unname(readGrgList.keep))
    }
    if(!use.names.OG) {names(readGrgList) <- NULL }
    return(readGrgList)
}

