#' Function to prepare reads for processing from bam file
#' @param bamFile bamFile
#' @inheritParams bambu
#' @noRd
prepareDataFromBam <- function(bamFile, yieldSize=NULL, verbose = FALSE, ncore = 1) {
  
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
  
  # parallel processing of single files by reading chromosomes separately
  if(ncore>1){ 
    
    chr <- scanBamHeader(bamFile)[[1]]
    chrRanges <- GRanges(seqnames=names(chr), ranges=IRanges(start=1, end=chr))
    bpParameters <- BiocParallel::bpparam()
    #===# set parallel options: If more CPUs than samples available, use parallel computing on each sample, otherwise use parallel to distribute samples (more efficient)
    bpParameters$workers <- ncore
    readGrgList <- BiocParallel::bplapply(1:length(chr),
                            helpFun,
                            chrRanges=chrRanges,
                            bamFile=bamFile,
                            BPPARAM = bpParameters)
   } else {
  
  
  bf <- open(bamFile)
  #readGrgList <- GenomicRanges::GRangesList()
  
  readGrgList <- list()
  counter <- 1
  
  while(Rsamtools::isIncomplete(bf)){
    readGrgList[[counter]] <- GenomicAlignments::grglist(GenomicAlignments::readGAlignments(bf,
                                                                                            param=Rsamtools::ScanBamParam(flag=Rsamtools::scanBamFlag(isSecondaryAlignment=FALSE)),
                                                                                            use.names=FALSE))
    
    # readGrgList<-c(readGrgList,GenomicAlignments::grglist(reads))
    if(verbose) show(min(length(readGrgList), counter* yieldSize, na.rm=T))
    counter <- counter + 1
  }
  
  close(bf)
  }
  if(length(readGrgList)>1) {
    readGrgList <-  do.call(c, readGrgList)
  } else   {
    readGrgList <- readGrgList[[1]]
  }
  readGrgList <- readGrgList[GenomicRanges::width(readGrgList)>1]  # remove microexons of width 1bp from list
  mcols(readGrgList)$id = 1:length(readGrgList) 
  return(readGrgList)
}

 helpFun<- function(chr, chrRanges, bamFile){
  return(GenomicAlignments::grglist(readGAlignments(file=bamFile,
                                                    param=Rsamtools::ScanBamParam(flag=Rsamtools::scanBamFlag(isSecondaryAlignment=FALSE),
                                                                                  which=chrRanges[chr]),
                                                    use.names=FALSE)))
  }