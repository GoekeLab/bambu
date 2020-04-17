#' Function to prepare reads for processing from bam file
#' @noRd
prepareDataFromBam <- function(bamFile, yieldSize=NULL, verbose = FALSE) {
  ## don't use names for reads, use index instead? Names are require a lot of memory

  if(class(bamFile)=='BamFile') {
    if(!is.null(yieldSize)) {
      Rsamtools::yieldSize(bamFile) <- yieldSize
    } else {
      yieldSize <- Rsamtools::yieldSize(bamFile)
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
  readGrglist <- GenomicRanges::GRangesList()
  counter <- 1

  while(Rsamtools::isIncomplete(bf)){
    # reads=readGAlignments(bamFile, param=ScanBamParam(flag=scanBamFlag(isSecondaryAlignment=FALSE),what = scanBamWhat()), use.names=T)
    reads=GenomicAlignments::readGAlignments(bf,
                          param=Rsamtools::ScanBamParam(flag=Rsamtools::scanBamFlag(isSecondaryAlignment=FALSE)),
                          use.names=TRUE)
    readGrglist<-c(readGrglist,GenomicAlignments::grglist(reads))
    if(verbose) show(min(length(reads), counter* yieldSize, na.rm=T))
    counter <- counter + 1
  }

  close(bf)
  readGrglist <- readGrglist[GenomicRanges::width(readGrglist)>1]  # remove microexons of width 1bp from list
  mcols(readGrglist)$qname = names(readGrglist)
  names(readGrglist) <- 1:length(readGrglist)  # names needed to be replaced as some reads are multiple times mapped (distinct parts of the read which are compatible)
  return(readGrglist)
}

