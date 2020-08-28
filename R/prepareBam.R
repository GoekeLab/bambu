#' Function to prepare reads for processing from bam file
#' @param bamFile bamFile
#' @inheritParams bambu
#' @noRd
prepareDataFromBam <- function(bamFile, yieldSize=NULL, verbose = FALSE, ncore = 1) {
  
  if(is(bamFile,'BamFile')) {
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
  #readGrgList <- GenomicRanges::GRangesList()
  
  readGrgList <- list()
  counter <- 1
  
  while(Rsamtools::isIncomplete(bf)){
    readGrgList[[counter]] <- GenomicAlignments::grglist(GenomicAlignments::readGAlignments(bf,
                                                                                            param=Rsamtools::ScanBamParam(flag=Rsamtools::scanBamFlag(isSecondaryAlignment=FALSE)),
                                                                                            use.names=FALSE))
    
    # readGrgList<-c(readGrgList,GenomicAlignments::grglist(reads))
    if(verbose) show(min(length(readGrgList), counter* yieldSize, na.rm=TRUE))
    counter <- counter + 1
  }
  
  close(bf)
 
  if(length(readGrgList)>1) {
    readGrgList <-  do.call(c, readGrgList)
  } else   {
    readGrgList <- readGrgList[[1]]
  }
  readGrgList <- readGrgList[GenomicRanges::width(readGrgList)>1]  # remove microexons of width 1bp from list
  mcols(readGrgList)$id = seq_along(readGrgList) 
  return(readGrgList)
}

 helpFun<- function(chr, chrRanges, bamFile){
  return(GenomicAlignments::grglist(readGAlignments(file=bamFile,
                                                    param=Rsamtools::ScanBamParam(flag=Rsamtools::scanBamFlag(isSecondaryAlignment=FALSE),
                                                                                  which=chrRanges[chr]),
                                                    use.names=FALSE)))
  }
