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
    #polyA variables
    searchLengthUnaligned = 20
    mat <- nucleotideSubstitutionMatrix(match = 1, mismatch = -2, baseOnly = TRUE)
    polyAPattern=DNAStringSet(paste0(rep('T',searchLengthUnaligned),collapse=''))  
    polyAPatternLong=DNAStringSet(paste0(rep('T',2000),collapse=''))  
    polyTPattern=DNAStringSet(paste0(rep('A',searchLengthUnaligned),collapse=''))  
    polyTPatternLong=DNAStringSet(paste0(rep('A',2000),collapse=''))  
    ###
    while (isIncomplete(bf)) {
        reads = readGAlignments(bf,
            param = ScanBamParam(flag =
                scanBamFlag(isSecondaryAlignment = FALSE), what = "seq"),
            use.names = use.names)
        readGrgList[[counter]] = grglist(reads)
        softClip5Prime <-pmax(0,as.numeric(gsub('^(\\d*)[S].*','\\1',GenomicAlignments::cigar(reads))), na.rm=T)
        softClip3Prime <-pmax(0,as.numeric(gsub('.*\\D(\\d*)[S]$','\\1',GenomicAlignments::cigar(reads))), na.rm=T)
        hardClip5Prime <-pmax(0,as.numeric(gsub('^(\\d*)[H].*','\\1',GenomicAlignments::cigar(reads))), na.rm=T)
        hardClip3Prime <-pmax(0,as.numeric(gsub('.*\\D(\\d*)[H]$','\\1',GenomicAlignments::cigar(reads))), na.rm=T)
        mcols(readGrgList[[counter]])$softClip5Prime = softClip5Prime
        mcols(readGrgList[[counter]])$softClip3Prime = softClip3Prime
        mcols(readGrgList[[counter]])$hardClip5Prime = hardClip5Prime
        mcols(readGrgList[[counter]])$hardClip3Prime = hardClip3Prime

        #todo: seperate this by alignment strand
        mcols(readGrgList[[counter]])$polyA5 = findPolyATail(mcols(reads)$seq, softClip5Prime, "start", polyAPattern, polyAPatternLong, mat)
        mcols(readGrgList[[counter]])$polyA3 = findPolyATail(mcols(reads)$seq, softClip3Prime, "end", polyTPattern, polyTPatternLong, mat)
        # mcols(readGrgList[[counter]])$polyA = rep(NA, length(reads))
        # mcols(readGrgList[[counter]])$polyA[strand(reads)=="+"] = findPolyATail(reads[strand(reads)=="+"], softClip5Prime, polyAPattern, polyAPatternLong, mat)
        # mcols(readGrgList[[counter]])$polyA[strand(reads)=="-"] = findPolyATail(reads[strand(reads)=="-"], softClip3Prime, polyTPattern, polyTPatternLong, mat)
        
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

barcodeAlignmentExtended <-  function(pattern, subject,type='local-global',...)  ## should replace function above, contains percent identity (pid)
{
  data <- matrix(NA,ncol=4,nrow=length(subject))
  colnames(data) <- c('score','pid','start', 'end')
  #  seq.align <- pairwiseAlignment(pattern=rep(DNAStringSet(pattern),length(subject)),subject=subject,type='global-local',...)
  seq.align <- pairwiseAlignment(pattern=subject,subject=pattern,type=type,...)
  data[,'score'] <- score(seq.align)
  data[,'pid'] <- pid(seq.align)
  data[,'start'] <- start(pattern(seq.align))
  data[,'end'] <- end(pattern(seq.align))
  return(data)
}

