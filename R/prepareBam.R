#' Function to prepare reads for processing from bam file
#'@title PREPAREDATAFROMBAM
#'@param bamFile
#'@param yieldSize
#'@export
prepareDataFromBam <- function(bamFile, yieldSize=NULL) {
  ## Todo: optimise which data is essential
  ## Todo: don't use names for reads, use index instead? Names are require a lot of memory


  if(class(bamFile)=='BamFile') {
    if(!is.null(yieldSize)) {
      yieldSize(bamFile) <- yieldSize
    } else {
      yieldSize <- yieldSize(bamFile)
    }

  }else if(!grepl('.bam',bamFile)){
     stop("Bam file is missing from arguments.")

    }else{
    if(is.null(yieldSize)) {
      yieldSize <- NA
    }
    bamFile <- Rsamtools::BamFile(bamFile, yieldSize = yieldSize)
  }
  bf <- open(bamFile)
  #readJunctions <- GRangesList()
  readGrglist <- GenomicRanges::GRangesList()
  # startRanges <- GRanges()
  # endRanges <- GRanges()
  # readAlignmentStrand <- NULL
  # readIdToName <- NULL
  counter <- 1

  while(Rsamtools::isIncomplete(bf)){
    # reads=readGAlignments(bamFile, param=ScanBamParam(flag=scanBamFlag(isSecondaryAlignment=FALSE),what = scanBamWhat()), use.names=T)
    reads=GenomicAlignments::readGAlignments(bf,
                          param=Rsamtools::ScanBamParam(flag=Rsamtools::scanBamFlag(isSecondaryAlignment=FALSE)), ## whether supplementary alignments would be considered as well?? To confirm with Jonathan.
                          use.names=TRUE)
    #readJunctions <- c(readJunctions,junctions(reads))
    readGrglist<-c(readGrglist,GenomicAlignments::grglist(reads))
    # startRanges <- c(startRanges,GRanges(seqnames=seqnames(reads),ranges=IRanges(start=start(reads),end=start(reads))))  ## can be extracted from readGrglist, remove here? range(readGrglist[1:10])
    # endRanges <- c(endRanges,GRanges(seqnames=seqnames(reads),ranges=IRanges(start=end(reads),end=end(reads))))## can be extracted from readGrglist, remove here? range(readGrglist[1:10])
    # readIdToName <- c(readIdToName,names(reads))## replace with integer Ids or use real names
    # readAlignmentStrand <- c(readAlignmentStrand,as.character(strand(reads)))
    show(counter* yieldSize)
    counter <- counter + 1
  }

  close(bf)
  readGrglist <- readGrglist[GenomicRanges::width(readGrglist)>1]  # remove microexons Of width 1bp from list
  readNames <- names(readGrglist)
  names(readGrglist) <- 1:length(readGrglist)  # names needed to be replaced as some reads are multiple times mapped (distinct parts of the read which are compatible)

  # readJunctions <- myGaps(readGrglist)
  # names(readJunctions) <- paste('read',1:length(readJunctions),sep='.')  ## later: remove names (they blow up memory usage!), maybe hashmap? hm <- hashmap(1:nrow(data), data$names)
  # names(startRanges) <- names(readJunctions)## later: remove names (they blow up memory usage!) remove this object and replace with range(readGrglist) whenever needed
  # names(endRanges) <- names(readJunctions)## later: remove names (they blow up memory usage!)remove this object and replace with range(readGrglist) whenever needed
  # names(readIdToName) <- names(readJunctions)  ## requires a lot of memory, replace later
  # names(readAlignmentStrand) <- names(readJunctions)

  return(list(readGrglist=readGrglist, readNames=readNames)) # startRanges=startRanges, endRanges=endRanges))
}

