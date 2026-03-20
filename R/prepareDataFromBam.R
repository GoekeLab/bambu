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
    use.names.OG <- use.names
    if(!isFALSE(demultiplexed) | cleanReads) use.names <- TRUE
    # grepl(".[ct]sv$",demultiplexed)
    if(!is.logical(demultiplexed)){ # if demultiplexed argument is not logical value  
        if(file.exists(demultiplexed)){ # check if file path exists, it has to be a regular delimited file, can be compressed 
            readMap <- fread(demultiplexed,header = FALSE, data.table = FALSE) # changed function to more efficiently read in data and also correct error of "no lines available as input" for csv format input
        }else{
            stop("Provided barcode to map file does not exists! Please provide the correct path to demultiplex argument!")
        }
    }
    while (isIncomplete(bf)) {
        alignmentInfo <- readGAlignments(bf, param = ScanBamParam(tag = c("CB", "UB"), 
                                         flag = scanBamFlag(isSecondaryAlignment = FALSE)), 
                                         use.names = use.names)
        readGrgList[[counter]] <-grglist(alignmentInfo)
        if (!isFALSE(demultiplexed)){ # if demultiplexed is TRUE or a string path 
            if(isTRUE(demultiplexed)){ # if demultiplexed is TRUE
      
                mcols(readGrgList[[counter]])$CB <- case_when(!is.na(mcols(alignmentInfo)$CB) ~ mcols(alignmentInfo)$CB, 
                                                              grepl("^[^_]+_[^#]+#", names(readGrgList[[counter]]), perl = TRUE) ~ sub("_.*", "", names(readGrgList[[counter]])), # a checkpoint to see whether CB is contained in the name, with specific format CB_UMI#READNAME 
                                                              TRUE ~ NA) 

                mcols(readGrgList[[counter]])$UMI <- case_when(!is.na(mcols(alignmentInfo)$UB) ~ mcols(alignmentInfo)$UB,
                                                               grepl("^[^_]+_[^#]+#", names(readGrgList[[counter]]), perl = TRUE) ~ sub("^[^_]+_([^#]+)#.*$", "\\1", names(readGrgList[[counter]])), # a checkpoint to see whether UMI is contained in the name, with specific format CB_UMI#READNAME, 
                                                               TRUE ~ NA) 
            } else{ # if demultiplexed is a string path
                mcols(readGrgList[[counter]])$CB <- NA
                mcols(readGrgList[[counter]])$UMI <- NA
                mcols(readGrgList[[counter]])$CB <- readMap[,2][match(names(readGrgList[[counter]]),readMap[,1])]
                if(ncol(readMap)>2){
                    mcols(readGrgList[[counter]])$UMI <- readMap[,3][match(names(readGrgList[[counter]]),readMap[,1])]
                }
            }
            cells <- unique(c(cells, mcols(readGrgList[[counter]])$CB))
            mcols(readGrgList[[counter]])$CB <- factor(mcols(readGrgList[[counter]])$CB, levels = cells)
            umi <- unique(c(umi, mcols(readGrgList[[counter]])$UMI))
            mcols(readGrgList[[counter]])$UMI <- factor(mcols(readGrgList[[counter]])$UMI, levels = umi)
        }
        if(cleanReads){
            softClip5Prime <- clipFunction(cigarData = GenomicAlignments::cigar(alignmentInfo), grep_pattern = '^(\\d*)[S].*', replace_pattern = '\\1')
            softClip3Prime <- clipFunction(cigarData = GenomicAlignments::cigar(alignmentInfo), grep_pattern = '.*\\D(\\d*)[S]$', replace_pattern = '\\1')
            hardClip5Prime <- clipFunction(cigarData = GenomicAlignments::cigar(alignmentInfo), grep_pattern = '^(\\d*)[H].*', replace_pattern = '\\1')
            hardClip3Prime <- clipFunction(cigarData = GenomicAlignments::cigar(alignmentInfo), grep_pattern = '.*\\D(\\d*)[H]$', replace_pattern = '\\1')
            # softClip5Prime <-suppressWarnings(pmax(0,as.numeric(gsub('^(\\d*)[S].*','\\1',GenomicAlignments::cigar(alignmentInfo))), na.rm=T))
            # softClip3Prime <-suppressWarnings(pmax(0,as.numeric(gsub('.*\\D(\\d*)[S]$','\\1',GenomicAlignments::cigar(alignmentInfo))), na.rm=T))
            # hardClip5Prime <-suppressWarnings(pmax(0,as.numeric(gsub('^(\\d*)[H].*','\\1',GenomicAlignments::cigar(alignmentInfo))), na.rm=T))
            # hardClip3Prime <-suppressWarnings(pmax(0,as.numeric(gsub('.*\\D(\\d*)[H]$','\\1',GenomicAlignments::cigar(alignmentInfo))), na.rm=T))
            mcols(readGrgList[[counter]])$clip5Prime <- pmax(softClip5Prime, hardClip5Prime)
            mcols(readGrgList[[counter]])$clip3Prime <- pmax(softClip3Prime, hardClip3Prime)
            rev <- as.vector(strand(alignmentInfo) == '-')
            rev2 <- grepl("_-.+of", names(alignmentInfo))
            temp <- mcols(readGrgList[[counter]])$clip5Prime
            mcols(readGrgList[[counter]])$clip5Prime[rev != rev2] <- mcols(readGrgList[[counter]])$clip3Prime[rev != rev2]
            mcols(readGrgList[[counter]])$clip3Prime[rev != rev2] <- temp[rev != rev2]
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
    readGrgList <- readGrgList <- readGrgList[sum(width(readGrgList)) > 1]
    numNoCBs <- sum(is.na(mcols(readGrgList)$CB))
    if(numNoCBs > 0){
        message("Removing ", numNoCBs, " reads that were not assigned barcodes. If this is unexpected check the barcode map input")
        readGrgList <- readGrgList[!is.na(mcols(readGrgList)$CB)]
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
        df <- data.frame(name = names(readGrgList), 
            clip5 = mcols(readGrgList)$clip5Prime)
        df <- df %>% mutate(id = row_number()) %>% group_by(name) %>% summarise(primary.id = id[which.min(clip5)])
        readGrgList <- readGrgList[df$primary.id]
        end.ptm <- proc.time()
        message("Primary alignment selection time ", round((end.ptm - start.ptm)[3] / 60, 3), " mins.")
        
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
