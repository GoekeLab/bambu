######
## Todo (2019-10-25 JG):
## combine transcripts from multiple runs to make the transcript IDs comparable
## assign tx Ids only after filtering
## clean up code

## 1) assign read classes to genes (read classes might be assigned to a gene but not to a transcript, for example very low quality reads)
## 2) among read classes assigned to genes but not to transcripts, identify new transcripts and annotated by type (e.g. new first exons, new last exon, alternative splicing), filter for high quality read classes only
## 3) among read classes not assigned to any gene, filter high quality read classes and add as new genes
## 4) extend annotations with new transcripts/genes
## 5) recalculate distance table for combined annotations+new transcripts/genes (optional?)
## 6) summarise statistics: how many reads assigned to known transcripts, how many assigned to new transcripts, how many new transcripts, ... some other meaningful statistics?
## general) clean up code

## todo: create annotation object to avoid multiple copying and processing of the same data (eg resorting exonsByTx etc)
## todo: add option to output reads corrected and stranded as bed file
## todo: check code for unncecceary commands, columns that are not needed, ways to improve speed remove redundancy
## include test cases to measure accuracy against default??
## better strand prediction for transcripts without clear strand from junction
## include no junction reads
## assign all reads to all transcripts (CY)


#####
#' Isoform reconstruction using genomic alignments
#'@title ISOFORM RECONSTRUCTION
#'@title ISOFORM RECONSTRUCTION
#'@description
#'
#'@param bamFile A string variable that indicates the path to genome bam file.
#'@param txdb txdb object.
#'@param genomeFAFile A string variable that indicates the path to genome annotation .fa file.
#'@param stranded A logical variable that indicates whether the experiment is a stranded run or non-stranded run (default to non-stranded).
#'@param protocol A string variable indicates the sequencing protocol used in this experiment, one of the following values: not used at the momen (default to NULL).
#'@param minimumReadSupport A integer value that indicates that least number of reads needed to support, default to 2.
#'@param prefix prefix for new gene names (for combining multiple runs gene Ids can be made unique and merged later)
#'@param minimumTxFraction A numeric value indicates the minimum transcript coverage fraction needed to support, default to 0.02, i.e., 2%.
#'@param yieldSize A numeric value indicates the yieldSize.
#'@export
#'@examples
#' \dontrun{
#'  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
#'  samplefile <- system.file("extdata","Homo_sapiens.GRCh38.91.annotations-txdb.sqlite", package = "bamboo")
#'  txdb <- loadDb(samplefile) #
#'  test.bam <- system.file("extdata", "GIS_HepG2_cDNAStranded_Rep5-Run4_chr9_108865774_109014097.bam", package = "bamboo")
#'  standardJunctionModelFile <- system.file("extdata", "standardJunctionModel-temp.rds", package = "bamboo")
#'  fa.file <- system.file("extdata", "Homo_sapiens.GRCh38.dna_sm.primary_assembly_chr9.fa.gz", package = "bamboo")
#'  isore(bamFile = test.bam,  txdb = txdb,genomeFA = FaFile(fa.file))
#'  }
isore <- function(bamFile,
                  txdb=NULL, ##CY: this should all be based on R objects in memory, not files (if possible)
                  txdbTablesList=NULL, ## optional  ## will save a lot of time for multiple data sets
                  genomeDB=NULL, ## to be deleted   ## is required to avoid providing a fasta file with genome sequence, helpful for most users
                  genomeFA=NULL, ## genome FA file, should be in .fa format
                  stranded=FALSE,
                  protocol=NULL,
                  prefix='',  ## prefix for new gene names (for combining multiple runs gene Ids can be made unique and merged later)
                  minimumReadSupport=2,
                  minimumTxFraction = 0.02,
                  yieldSize = NULL,
                  quickMode = FALSE)
{
  show('deprecated - use bamboo directly')
}


isore.preprocessBam <- function(bamFile, yieldSize = NULL){

  cat('### load data ### \n')
  start.ptm <- proc.time()
  ## create BamFile object from character ##
  if(class(bamFile)=='BamFile') {
    if(!is.null(yieldSize)) {
      yieldSize(bamFile) <- yieldSize
    } else {
      yieldSize <- Rsamtools::yieldSize(bamFile)
    }
  }else if(!grepl('.bam',bamFile)){
    stop("Bam file is missing from arguments.")
  }else{
    if(is.null(yieldSize)) {
      yieldSize <- NA
    }

  }

  readData <- prepareDataFromBam(bamFile)
  end.ptm <- proc.time()
  cat(paste0('Finished loading data in ', round((end.ptm-start.ptm)[3]/60,1), ' mins. \n'))
  return(readData)
}

isore.constructReadClasses <- function(readGrgList,
                                       runName='sample1',
                                       txdbTablesList, ## has to be provided (function should be called through bamboo, so is optional through that main function)
                                       genomeDB=NULL, ## is required to avoid providing a fasta file with genome sequence, helpful for most users
                                       genomeFA=NULL, ## genome FA file, should be in .fa format
                                       stranded=FALSE,
                                       protocol=NULL,
                                       prefix='',  ## prefix for new gene names (for combining multiple runs gene Ids can be made unique and merged later)
                                       minimumReadSupport=2,
                                       minimumTxFraction = 0.02,
                                       yieldSize = NULL,
                                       quickMode = FALSE){



 ## todo: which preprocessed junction correction model to use?

  if(is.null(genomeFA)){
    stop("GenomeFA file is missing.")
  }else if(class(genomeFA) != 'FaFile'){
    if(!grepl('.fa',genomeFA)){
      stop("GenomeFA file is missing.")
    }else{
      genomeFA <- Rsamtools::FaFile(genomeFA)
    }
  }


  # cat('### load data ### \n')
  # start.ptm <- proc.time()
  # ## create BamFile object from character ##
  # if(class(bamFile)=='BamFile') {
  #   if(!is.null(yieldSize)) {
  #     yieldSize(bamFile) <- yieldSize
  #   } else {
  #     yieldSize <- Rsamtools::yieldSize(bamFile)
  #   }
  # }else if(!grepl('.bam',bamFile)){
  #   stop("Bam file is missing from arguments.")
  # }else{
  #   if(is.null(yieldSize)) {
  #     yieldSize <- NA
  #   }
  #   bamFile <- Rsamtools::BamFile(bamFile, yieldSize = yieldSize)
  # }
  #
  # readGrglist <- prepareDataFromBam(bamFile)
  # end.ptm <- proc.time()
  # cat(paste0('Finished loading data in ', round((end.ptm-start.ptm)[3]/60,1), ' mins. \n'))


  unlisted_junctions <- unlist(myGaps(readGrgList))
  cat('### create junction list with splice motif ### \n')
  start.ptm <- proc.time()
  uniqueJunctions <- createJunctionTable(unlisted_junctions, genomeDB = genomeDB, genomeFA=genomeFA)
  end.ptm <- proc.time()
  cat(paste0('Finished creating junction list with splice motif in ', round((end.ptm-start.ptm)[3]/60,1), ' mins. \n'))



  cat('### infer strand/strand correction of junctions ### \n')
  junctionTables <- junctionStrandCorrection(uniqueJunctions, unlisted_junctions, txdbTablesList[['intronsByTxEns']], stranded=stranded)
  uniqueJunctions <- junctionTables[[1]]
  unlisted_junctions <- junctionTables[[2]]
  rm(junctionTables)
  gc()

  cat('### find annotated introns ### \n')
  uniqueJunctions$annotatedJunction <- (!is.na(GenomicRanges::match(uniqueJunctions, unique(txdbTablesList[['unlisted_introns']]))))

  # Indicator: is the junction start annotated as a intron start?
  annotatedStart <- tapply(uniqueJunctions$annotatedJunction,  uniqueJunctions$junctionStartName,sum)>0
  uniqueJunctions$annotatedStart <- annotatedStart[uniqueJunctions$junctionStartName]
  rm(annotatedStart)
  gc()

  # Indicator: is the junction end annotated as a intron end?
  annotatedEnd <- tapply(uniqueJunctions$annotatedJunction, uniqueJunctions$junctionEndName,sum)>0
  uniqueJunctions$annotatedEnd <- annotatedEnd[uniqueJunctions$junctionEndName]
  rm(annotatedEnd)
  gc()

  cat('### build model to predict true splice sites ### \n')
  start.ptm <- proc.time()
  if(sum(uniqueJunctions$annotatedJunction)>4500 &sum(!uniqueJunctions$annotatedJunction)>5000){  ## note: these are thresholds that should be adjusted, or changed. Also can look into the code to find out what is the issue, probably number of training data per strand?
    predictSpliceSites <- predictSpliceJunctions(uniqueJunctions,junctionModel = NULL)
    uniqueJunctions=predictSpliceSites[[1]]
    junctionModel=predictSpliceSites[[2]]
  } else {
    junctionModel = standardJunctionModels_temp
    predictSpliceSites <- predictSpliceJunctions(uniqueJunctions,junctionModel = junctionModel)
    uniqueJunctions=predictSpliceSites[[1]]
    show('Warning: junction correction with not enough data, precalculated model is used')
  }
  rm(predictSpliceSites)  # clean up should be done more efficiently
  gc()
  end.ptm <- proc.time()
  cat(paste0('Model to predict true splice sites built in ', round((end.ptm-start.ptm)[3]/60,1), ' mins. \n'))


  cat('### correct junctions based on set of high confidence junctions ### \n')
  start.ptm <- proc.time()
  uniqueJunctions <- findHighConfidenceJunctions(uniqueJunctions, junctionModel)
  uniqueJunctions$mergedHighConfJunctionIdAll_noNA <- uniqueJunctions$mergedHighConfJunctionId
  uniqueJunctions$mergedHighConfJunctionIdAll_noNA[is.na(uniqueJunctions$mergedHighConfJunctionId)] <- names(uniqueJunctions[is.na(uniqueJunctions$mergedHighConfJunctionId)])
  uniqueJunctions$strand.mergedHighConfJunction <- as.character(strand(uniqueJunctions[uniqueJunctions$mergedHighConfJunctionIdAll_noNA]))
  end.ptm <- proc.time()
  cat(paste0('Finished correcting junction based on set of high confidence junctions in ', round((end.ptm-start.ptm)[3]/60,1), ' mins. \n'))
  rm(junctionModel)
  gc()

  cat('### create transcript models (read classes) from spliced reads ### \n')
  start.ptm <- proc.time()
  readClassListSpliced <- constructSplicedReadClassTables(uniqueJunctions, unlisted_junctions, readGrgList, mcols(readGrgList)$qname, quickMode = quickMode)  ## speed up this function ##
  end.ptm <- proc.time()
  cat(paste0('Finished create transcript models (read classes) for reads with spliced junctions in ', round((end.ptm-start.ptm)[3]/60,1), ' mins. \n'))
  rm(list = c('uniqueJunctions','unlisted_junctions'))
  gc()

  cat('### create single exon transcript models (read classes) ### \n')
  start.ptm <- proc.time()

  singleExonReads <- unlist(readGrgList[elementNROWS(readGrgList)==1])
  referenceExons<-unique(c(granges(unlist(readClassListSpliced[['exonsByReadClass']][readClassListSpliced$readClassTable$confidenceType=='highConfidenceJunctionReads' & readClassListSpliced$readClassTable$strand!='*'])), granges(unlist(txdbTablesList[['exonsByTx']]))))
  readClassListUnsplicedWithAnnotation <- constructUnsplicedReadClasses(singleExonReads,referenceExons, mcols(readGrgList)$qname, confidenceType = 'unsplicedWithin', prefix='unsplicedWithin',stranded=stranded) ### change txId OK

  singleExonReadsOutside <- singleExonReads[!(mcols(readGrgList)$qname[as.integer(names(singleExonReads))] %in% readClassListUnsplicedWithAnnotation$readTable$readId)]
  rm(list=c('singleExonReads'))
  gc()

  combinedSingleExonRanges <- reduce(singleExonReadsOutside,ignore.strand=!stranded)
  readClassListUnsplicedReduced <- constructUnsplicedReadClasses(singleExonReadsOutside,combinedSingleExonRanges, mcols(readGrgList)$qname, confidenceType = 'unsplicedNew', prefix='unsplicedNew',stranded=stranded)
  rm(list = c('singleExonReadsOutside','combinedSingleExonRanges','readGrgList'))
  gc()

  end.ptm <- proc.time()
  cat(paste0('Finished create single exon transcript models (read classes) in ', round((end.ptm-start.ptm)[3]/60,1), ' mins. \n'))


  exonsByReadClass = c(readClassListSpliced$exonsByReadClass, readClassListUnsplicedWithAnnotation$exonsByReadClass, readClassListUnsplicedReduced$exonsByReadClass)
  readClassTable = rbind(readClassListSpliced$readClassTable, readClassListUnsplicedWithAnnotation$readClassTable, readClassListUnsplicedReduced$readClassTable)
  #readTable = rbind(readClassListSpliced$readTable, readClassListUnsplicedWithAnnotation$readTable, readClassListUnsplicedReduced$readTable)
  rm(list=c('readClassListSpliced','readClassListUnsplicedWithAnnotation','readClassListUnsplicedReduced'))
  gc()


  #bamFile.basename <- tools::file_path_sans_ext(basename(BiocGenerics::path(bamFile)))
  counts <- matrix(readClassTable$readCount, dimnames = list(names(exonsByReadClass), runName))
  colDataDf <- DataFrame(name=runName, row.names=runName)
  mcols(exonsByReadClass) <- select(readClassTable, chr.rc = chr, strand.rc=strand, intronStarts, intronEnds, confidenceType)
  # readTable is currently not returned
  se <- SummarizedExperiment::SummarizedExperiment(assays=SimpleList(counts=counts),
                                                   rowRanges = exonsByReadClass,
                                                   colData = colDataDf)


  rm(list=c('counts','exonsByReadClass','readClassTable'))
  return(se)
}


isore.combineTranscriptCandidates <- function(readClassSe, readClassSeRef=NULL, stranded=FALSE){
# seListTMP=seList[[3]]
# readClassSe <- seListTMP[grepl('highConfidenceJunctionReads|unsplicedNew', rowData(seListTMP)$confidenceType),]
  show('combine new transcript candidates')
  if(is.null(readClassSeRef)){  #if no reerence object is given, create one from a readClassSe object
    counts <- assays(readClassSe)$counts
    start <- matrix(min(start(rowRanges(readClassSe))), dimnames = dimnames(counts))
    end <- matrix(max(end(rowRanges(readClassSe))), dimnames = dimnames(counts))
    rowData <- as_tibble(rowData(readClassSe))
    rowData$start <- rowMins(start)
    rowData$end <- rowMaxs(end)
    rowData <- rowData %>% select(chr=chr.rc, start, end, strand=strand.rc, intronStarts, intronEnds, confidenceType)

    readClassSeRef <- SummarizedExperiment::SummarizedExperiment(assays=SimpleList(counts=counts, start=start, end=end),
                                                                 rowData = rowData,
                                                                 colData = colData(readClassSe))

  }else {
    colDataCombined <- rbind(colData(readClassSeRef), colData(readClassSe))

    readClassSeRefTBL <- as_tibble(rowData(readClassSeRef), rownames='id')
    readClassSeTBL <- as_tibble(rowData(readClassSe), rownames='id') %>%
                                mutate(start=min(start(rowRanges(readClassSe))),
                                       end=max(end(rowRanges(readClassSe))))

  rowData.spliced <- full_join(filter(readClassSeRefTBL, confidenceType=='highConfidenceJunctionReads'),
                          filter(readClassSeTBL, confidenceType=='highConfidenceJunctionReads'),
                          by=c('chr'='chr.rc','strand'='strand.rc','intronStarts', 'intronEnds'), suffix=c('.ref','.new'))

  #create first SE object for spliced Tx


  counts.splicedRef <- matrix(0, dimnames=list(1:nrow(rowData.spliced),rownames(colData(readClassSeRef))), ncol=nrow(colData(readClassSeRef)), nrow = nrow(rowData.spliced))
  start.splicedRef <- matrix(NA, dimnames=list(1:nrow(rowData.spliced),rownames(colData(readClassSeRef))), ncol=nrow(colData(readClassSeRef)), nrow = nrow(rowData.spliced))
  end.splicedRef <- start.splicedRef

  counts.splicedNew <- matrix(0, dimnames=list(1:nrow(rowData.spliced),rownames(colData(readClassSe))), ncol=nrow(colData(readClassSe)), nrow = nrow(rowData.spliced))
  start.splicedNew <- matrix(NA, dimnames=list(1:nrow(rowData.spliced),rownames(colData(readClassSe))), ncol=nrow(colData(readClassSe)), nrow = nrow(rowData.spliced))
  end.splicedNew <- start.splicedNew


  counts.splicedRef[!is.na(rowData.spliced$id.ref), ] <- as.matrix(assays(readClassSeRef)$counts[rowData.spliced$id.ref[!is.na(rowData.spliced$id.ref)],])

  start.splicedRef[!is.na(rowData.spliced$id.ref), ] <- as.matrix(assays(readClassSeRef)$start[rowData.spliced$id.ref[!is.na(rowData.spliced$id.ref)],])
  end.splicedRef[!is.na(rowData.spliced$id.ref), ] <- as.matrix(assays(readClassSeRef)$end[rowData.spliced$id.ref[!is.na(rowData.spliced$id.ref)],])


  counts.splicedNew[!is.na(rowData.spliced$id.new), ] <- as.matrix(assays(readClassSe)$counts[rowData.spliced$id.new[!is.na(rowData.spliced$id.new)],])
  start.splicedNew[!is.na(rowData.spliced$id.new), ] <- as.matrix(rowData.spliced[!is.na(rowData.spliced$id.new),'start.new'])
  end.splicedNew[!is.na(rowData.spliced$id.new), ] <- as.matrix(rowData.spliced[!is.na(rowData.spliced$id.new),'end.new'])



  counts.spliced <- cbind(counts.splicedRef, counts.splicedNew)
  start.spliced <- cbind(start.splicedRef, start.splicedNew)
  end.spliced <- cbind(end.splicedRef, end.splicedNew)
  rowData.spliced$start <- rowMins(start.spliced, na.rm=T)
  rowData.spliced$end <- rowMaxs(end.spliced, na.rm=T)
  rowData.spliced <- select(rowData.spliced, chr, start, end, strand, intronStarts, intronEnds) %>%
                      mutate(confidenceType='highConfidenceJunctionReads')

  se.spliced <- SummarizedExperiment::SummarizedExperiment(assays=SimpleList(counts=counts.spliced, start=start.spliced, end=end.spliced),
                                                               rowData = rowData.spliced,
                                                               colData = colDataCombined)



  ## create second SE object for unspliced Tx


  readClassSeRefTBL.unspliced <- filter(readClassSeRefTBL, confidenceType=='unsplicedNew')
  readClassSeTBL.unspliced <- filter(readClassSeTBL, confidenceType=='unsplicedNew')
  unsplicedRangesRef=GRanges(seqnames=readClassSeRefTBL.unspliced$chr,
                             ranges=IRanges(start=readClassSeRefTBL.unspliced$start,
                                            end=readClassSeRefTBL.unspliced$end),
                             strand=readClassSeRefTBL.unspliced$strand)
  unsplicedRangesNew=GRanges(seqnames=readClassSeTBL.unspliced$chr.rc,
                             ranges=IRanges(start=readClassSeTBL.unspliced$start,
                                            end=readClassSeTBL.unspliced$end),
                             strand=readClassSeTBL.unspliced$strand.rc)



  combinedSingleExonRanges <- reduce(c(unsplicedRangesRef,unsplicedRangesNew), ignore.strand=!stranded)

  rowData.unspliced <- as_tibble(as.data.frame(combinedSingleExonRanges, stringsAsFactors = FALSE)) %>%
    mutate_if(is.factor, as.character) %>%
    select(chr=seqnames, start, end, strand=strand) %>%
    mutate(intronStarts=NA, intronEnds=NA, confidenceType='unsplicedNew')

  overlapRefToCombined <-findOverlaps(unsplicedRangesRef,combinedSingleExonRanges, type='within', ignore.strand=!stranded, select='first')
  overlapNewToCombined <-findOverlaps(unsplicedRangesNew,combinedSingleExonRanges, type='within', ignore.strand=!stranded, select='first')

  countsRef.unspliced <- as_tibble(assays(readClassSeRef)[['counts']])[rowData(readClassSeRef)$confidenceType=='unsplicedNew',] %>%
                          mutate(index=overlapRefToCombined) %>%
                          group_by(index) %>%
                          summarise_all(sum)

  countsNew.unspliced <- as_tibble(assays(readClassSe)[['counts']])[rowData(readClassSe)$confidenceType=='unsplicedNew',] %>%
                          mutate(index=overlapNewToCombined) %>%
                          group_by(index) %>%
                          summarise_all(sum)

  startRef.unspliced <- as_tibble(assays(readClassSeRef)[['start']])[rowData(readClassSeRef)$confidenceType=='unsplicedNew',] %>%
                          mutate(index=overlapRefToCombined) %>%
                          group_by(index) %>%
                          summarise_all(min)

  startNew.unspliced <- readClassSeTBL %>% filter(confidenceType=='unsplicedNew') %>%
                                            select(start) %>%
                                            mutate(index=overlapNewToCombined) %>%
                                            group_by(index) %>%
                                            summarise_all(min)

  endRef.unspliced <- as_tibble(assays(readClassSeRef)[['end']])[rowData(readClassSeRef)$confidenceType=='unsplicedNew',] %>%
                        mutate(index=overlapRefToCombined) %>%
                        group_by(index) %>%
                        summarise_all(max)

  endNew.unspliced <- readClassSeTBL %>% filter(confidenceType=='unsplicedNew') %>%
                                          select(end) %>%
                                          mutate(index=overlapNewToCombined) %>%
                                          group_by(index) %>%
                                          summarise_all(max)

  counts.unsplicedRef <- matrix(0, dimnames=list(1:nrow(rowData.unspliced),rownames(colData(readClassSeRef))), ncol=nrow(colData(readClassSeRef)), nrow = nrow(rowData.unspliced))
  start.unsplicedRef <- matrix(NA, dimnames=list(1:nrow(rowData.unspliced),rownames(colData(readClassSeRef))), ncol=nrow(colData(readClassSeRef)), nrow = nrow(rowData.unspliced))
  end.unsplicedRef <- start.unsplicedRef

  counts.unsplicedNew <- matrix(0, dimnames=list(1:nrow(rowData.unspliced),rownames(colData(readClassSe))), ncol=nrow(colData(readClassSe)), nrow = nrow(rowData.unspliced))
  start.unsplicedNew <- matrix(NA, dimnames=list(1:nrow(rowData.unspliced),rownames(colData(readClassSe))), ncol=nrow(colData(readClassSe)), nrow = nrow(rowData.unspliced))
  end.unsplicedNew <- start.unsplicedNew


  counts.unsplicedRef[countsRef.unspliced$index, ] <- as.matrix(countsRef.unspliced[,colnames(counts.unsplicedRef)])
  start.unsplicedRef[countsRef.unspliced$index, ] <- as.matrix(startRef.unspliced[,colnames(start.unsplicedRef)])
  end.unsplicedRef[countsRef.unspliced$index, ] <- as.matrix(endRef.unspliced[,colnames(end.unsplicedRef)])
  counts.unsplicedNew[countsNew.unspliced$index, ] <- as.matrix(countsNew.unspliced[,colnames(counts.unsplicedNew)])
  start.unsplicedNew[countsNew.unspliced$index, ] <- as.matrix(startNew.unspliced[,'start'])
  end.unsplicedNew[countsNew.unspliced$index, ] <- as.matrix(endNew.unspliced[,'end'])



  counts.unspliced <- cbind(counts.unsplicedRef, counts.unsplicedNew)
  start.unspliced <- cbind(start.unsplicedRef, start.unsplicedNew)
  end.unspliced <- cbind(end.unsplicedRef, end.unsplicedNew)

  rowData.unspliced <- as_tibble(data.frame(combinedSingleExonRanges)) %>%
                        select(chr=seqnames, start, end, strand=strand) %>%
                        mutate(intronStarts=NA, intronEnds=NA, confidenceType='unsplicedNew')


  se.unspliced <- SummarizedExperiment::SummarizedExperiment(assays=SimpleList(counts=counts.unspliced, start=start.unspliced, end=end.unspliced),
                                                               rowData = rowData.unspliced,
                                                               colData = colDataCombined)

  se.combined <- SummarizedExperiment::rbind(se.spliced,se.unspliced)
  rownames(se.combined) <- 1:nrow(se.combined)
  rm(se.spliced, se.unspliced)
  return(se.combined)
  }
}


isore.extendAnnotations <- function(se,
                                    annotationGrangesList,
                                    filterSubsetTx = TRUE, # filter to remove read classes which are a subset of known transcripts. Also remove transcripts which are a subset of new transcripts (?)
                                    minAbsoluteReadCount = 2,  # minimun read count to consider a read class valid in a sample
                                    minRelativeReadCountByGene = 0.05,  ## minimum relative read count per gene, highly expressed genes will have many high read count low relative abundance transcripts that can be filtered
                                    minSampleNumber = 2,  # minimum sample number with minimum read count
                                    minimumExonDistance = 15)  # minum distance to known transcript to be considered valid as new
  {
  show('not yet implemented')

  filterSet1=(rowSums(assays(combinedTxCandidates)$counts>=minAbsoluteReadCount)>=minSampleNumber)

  ##compare with annotations, filter out complete matches and subset matches, calculate gene counts, filter by relative count
  cat('### [TODO] [optional]  classify readClasses ### \n')
  start.ptm <- proc.time()
  ## TODO: classify all read classes
  ## categories:
  ## compatible
  ## subset
  ## new transcript within annotation

  #readClassListFull <- classifyReadClasses(readClassListFull)

  end.ptm <- proc.time()
  cat(paste0('[TODO] [optional]  Finished  classifying readClasses in ', round((end.ptm-start.ptm)[3]/60,1), ' mins. \n'))

  return(annotationGrangesList)
}

isore.estimateDistanceToAnnotations <- function(seReadClass, annotationGrangesList, stranded=FALSE, prefix=''){
  cat('### calculate distance of read classes to annotations, basic filter for read-tx assignments ### \n')
  start.ptm <- proc.time()

  exonsByReadClass = rowRanges(seReadClass)
  readClassTable=as_tibble(rowData(seReadClass), rownames='readClassId')

  ## note/todo: here the stranded mode should always be used, need to check that in unstranded mode, readClasses without strand information from splice sites are '*'
  ## if stranded mode is turned off, then filtering needs to be adjusted to first select strandedMatches
  ## might not be a big issue (not clear)
  distTable <- calculateDistToAnnotation(exonsByReadClass,annotationGrangesList,maxDist = 35, primarySecondaryDist = 5, ignore.strand= !stranded)  # [readClassListFull$txTable$confidenceType=='highConfidenceJunctionReads' ]   ### change txId OK
  ## this line is removed, counts are in assays(se)

  distTable$readCount = assays(seReadClass)$counts[distTable$readClassId,]  # should actually be stored in counts, but is here to  assign genes based on high read counts
  distTable <- left_join(distTable, dplyr::select(readClassTable, readClassId, confidenceType)) %>% mutate(relativeReadCount=readCount/txNumberFiltered)
  distTable <- left_join(distTable,  as_tibble(mcols(annotationGrangesList)), by=c('annotationTxId'='TXNAME')) ## note: gene id still not unique, might need to assign after EM using empty read classes
  end.ptm <- proc.time()
  cat(paste0('Finished calculating distance of read classes to annotations in ', round((end.ptm-start.ptm)[3]/60,1), ' mins. \n'))




  cat('### assign unmatched readClasses to new geneIds ### \n')
  start.ptm <- proc.time()

  readClassToGeneIdTable <- dplyr::select(distTable, readClassId, GENEID, readCount) %>%group_by(GENEID) %>% mutate(geneCount = sum(readCount)) %>% distinct() %>%group_by(readClassId) %>% filter(geneCount==max(geneCount)) %>%filter(row_number()==1) %>% dplyr::select(readClassId, geneId=GENEID) %>% ungroup()


  newGeneCandidates <- (!readClassTable$readClassId %in% readClassToGeneIdTable$readClassId)
  readClassToGeneIdTableNew <- assignNewGeneIds(exonsByReadClass[newGeneCandidates], prefix=prefix, minoverlap=5, ignore.strand=F)
  readClassGeneTable <- rbind(readClassToGeneIdTable,readClassToGeneIdTableNew)
  readClassTable <- left_join(readClassTable, readClassGeneTable)
  end.ptm <- proc.time()
  cat(paste0('Finished assigning unmatched readClasses to new geneIds in ', round((end.ptm-start.ptm)[3]/60,1), ' mins. \n'))

  rm(list = c('newGeneCandidates','readClassGeneTable'))
  gc()

  ## implement next
  ############### FROM HERE ##############
  ## 2020-01-30:
  ## implement multi sample mode and read class annotation and filtering
  ########################################

  ## optional for multi sample quantification/ reconstruction
  cat('### [TODO] [optional] Combine read classes from multiple samples ### \n')
  start.ptm <- proc.time()
  end.ptm <- proc.time()
  cat(paste0('[TODO] [optional] Finished combining read classes from multiple samples in ', round((end.ptm-start.ptm)[3]/60,1), ' mins. \n'))

  cat('### [TODO] filter read classes/single sample or multi sample mode ### \n')
  start.ptm <- proc.time()
  end.ptm <- proc.time()
  cat(paste0('[TODO] Finished filtering read classes in ', round((end.ptm-start.ptm)[3]/60,1), ' mins. \n'))


  cat('### [TODO] [optional]  classify readClasses ### \n')
  start.ptm <- proc.time()
  ## TODO: classify all read classes
  ## categories:
  ## compatible
  ## subset
  ## new transcript within annotation

  #readClassListFull <- classifyReadClasses(readClassListFull)

  end.ptm <- proc.time()
  cat(paste0('[TODO] [optional]  Finished  classifying readClasses in ', round((end.ptm-start.ptm)[3]/60,1), ' mins. \n'))


  cat('### Create summarizedExperiment output ### \n')
  ## This chunk of code should be able to produce output required for quantification:
  ## assays: readClass count with empty read class(final read class that is based on transcript combination,i.e., equivalent class)
  ## rowData: can be empty
  ## metadata: eqClass to tx assignment; distTable for each rc and eqClass


  # bamFile.basename <- tools::file_path_sans_ext(basename(path(bamFile)))
  # counts <- matrix(readClassTable$readCount, dimnames = list(names(exonsByReadClass), bamFile.basename))
  # colDataDf <- DataFrame(name=bamFile.basename, row.names=bamFile.basename)
  #
  # # readTable is currently not returned
  # se <- SummarizedExperiment::SummarizedExperiment(assays=SimpleList(counts=counts),
  #                                                  rowRanges = exonsByReadClass,
  #                                                  colData = colDataDf,
  #                                                  metadata=list(distTable = distTable,
  #                                                                readClassTable = readClassTable))
   metadata(seReadClass)<-list(distTable=distTable)
  rowData(seReadClass) <- readClassTable

  rm(list=c('distTable','exonsByReadClass','readClassTable'))
  return(seReadClass)
}

classifyReadClasses <- function(readClassList) {

  exByTx_singleBpStartEnd <- cutStartEndFromGrangesList(readClassListFull$exonsByReadClass)
  spliceOverlaps=findSpliceOverlapsQuick(exByTx_singleBpStartEnd,exByTx_singleBpStartEnd)
  spliceOverlapsSelected =spliceOverlaps[mcols(spliceOverlaps)$compatible==TRUE,]
  txIsSubsetOf <- countQueryHits(spliceOverlapsSelected)-1
  txHasSubsetIn <- countSubjectHits(spliceOverlapsSelected)-1

  readClassListFull$readClassTable$txIsSubsetOf <- txIsSubsetOf


  setIncompatible <- (readClassListFull$readClassTable$readClassId %in% distTable$readClassId[distTable$compatible==FALSE])
  spliceOverlapsUnexplained=findSpliceOverlapsQuick(exByTx_singleBpStartEnd[setIncompatible],exByTx_singleBpStartEnd[setIncompatible])
  spliceOverlapsSelected =spliceOverlapsUnexplained[mcols(spliceOverlapsUnexplained)$compatible==TRUE,]

  #     txIsSubsetOfUnexplained <- rep(0, length(exByTx_singleBpStartEnd))
  #     txIsSubsetOfUnexplained[setIncompatible] <- countQueryHits(spliceOverlapsSelected)-1
  txHasSubsetInUnexplained <- rep(0, length(exByTx_singleBpStartEnd))
  txHasSubsetInUnexplained[setIncompatible] <- countSubjectHits(spliceOverlapsSelected)-1

  # readClassListFull$readClassTable$txIsSubsetOfUnexplained <- txIsSubsetOfUnexplained
  readClassListFull$readClassTable$txHasSubsetInUnexplained <- txHasSubsetInUnexplained


  distTableByGene <- left_join(distTable, dplyr::select(readClassListFull$readClassTable, readClassId,txIsSubsetOf, txHasSubsetInUnexplained)) %>% group_by(GENEID) %>% mutate(geneCount = sum(relativeReadCount), geneCountFractionCompatible = sum(relativeReadCount * compatible)/ geneCount)

  ## here ##
  ## here: add longest transcripts first, then validate with number of explained reads, predict tss/tes is separate task (?)
  ## return only new transcripts/genes
  ## add visualisation for new transcripts
  ##as.data.frame(filter(distTableByGene, GENEID == 'ENSG00000107104', txIsSubsetOf==0, confidenceType=='highConfidenceJunctionReads', compatible==FALSE))


}


assignNewGeneIds <- function(exByTx, prefix='', minoverlap=5, ignore.strand=F){

  exonSelfOverlaps <- findOverlaps(exByTx,exByTx,select = 'all',minoverlap = minoverlap, ignore.strand=ignore.strand)
  hitObject = tbl_df(exonSelfOverlaps) %>% arrange(queryHits, subjectHits)
  candidateList <- hitObject %>%
                    group_by(queryHits) %>%
                    filter( queryHits <= min(subjectHits), queryHits != subjectHits) %>%
                    ungroup()

  filteredOverlapList <- hitObject %>% filter(queryHits < subjectHits)

  #temp <- candidateList  ## %>% filter(queryHits<= subjectHits) ## remove identical hits will lead to wrong results
  #temp <- hitObject
  rm(list=c('exonSelfOverlaps','hitObject'))
  gc()
  length_tmp = 1
  while(nrow(candidateList)>length_tmp) {  # loop to include overlappng read classes which are not in order
    length_tmp = nrow(candidateList)
    temp <- left_join(candidateList,filteredOverlapList,by=c("subjectHits"="queryHits")) %>%
            group_by(queryHits) %>% filter(! subjectHits.y %in% subjectHits, !is.na(subjectHits.y)) %>%
            ungroup %>%
            dplyr::select(queryHits,subjectHits.y) %>%
            distinct() %>%
            dplyr::rename(subjectHits=subjectHits.y)
    candidateList <- rbind(temp, candidateList)
    while(nrow(temp)>0) {
      show('annotated transcripts from unknown genes by new gene id')

      show(nrow(candidateList))
      temp= left_join(candidateList,filteredOverlapList,by=c("subjectHits"="queryHits")) %>%
              group_by(queryHits) %>%
              filter(! subjectHits.y %in% subjectHits, !is.na(subjectHits.y)) %>%
              ungroup %>%
              dplyr::select(queryHits,subjectHits.y) %>%
              distinct() %>%
              dplyr::rename(subjectHits=subjectHits.y)

      candidateList <- rbind(temp, candidateList)
    }
    #length_tmp = nrow(candidateList)
    show('second loop')
    show(length_tmp)
    tst <- candidateList %>%
            group_by(subjectHits) %>%
            mutate(subjectCount = n()) %>%
            group_by(queryHits) %>%
            filter(max(subjectCount)>1) %>%
            ungroup()

    temp2= inner_join(tst,tst,by=c("subjectHits"="subjectHits")) %>%
            filter(queryHits.x!=queryHits.y)  %>%
            mutate(queryHits = if_else(queryHits.x > queryHits.y, queryHits.y, queryHits.x),
                   subjectHits = if_else(queryHits.x > queryHits.y, queryHits.x, queryHits.y)) %>%
            dplyr::select(queryHits,subjectHits) %>%
            distinct()
    candidateList <-  distinct(rbind(temp2, candidateList))
  }

  candidateList <- candidateList %>%
                    filter(! queryHits %in% subjectHits) %>%
                    arrange(queryHits, subjectHits)
  idToAdd <- (which(!(1:length(exByTx) %in% unique(candidateList$subjectHits))))

  candidateList <- rbind(candidateList, tibble(queryHits=idToAdd, subjectHits=idToAdd)) %>%
                    arrange(queryHits, subjectHits) %>%
                    mutate(geneId = paste('gene',prefix,'.',queryHits,sep='')) %>%
                    dplyr::select(subjectHits, geneId)
  candidateList$readClassId <- names(exByTx)[candidateList$subjectHits]
  candidateList <- dplyr::select(candidateList, readClassId, geneId)
  return(candidateList)
}

getMinimumEqClassByTx <- function(exonsByTranscripts) {

  exByTxAnnotated_singleBpStartEnd <- cutStartEndFromGrangesList(exonsByTranscripts)  # estimate overlap only based on junctions
  spliceOverlaps=findSpliceOverlapsQuick(exByTxAnnotated_singleBpStartEnd,exByTxAnnotated_singleBpStartEnd)  ## identify transcripts which are compatbile with other transcripts (subsets by splice sites)
  spliceOverlapsSelected =spliceOverlaps[mcols(spliceOverlaps)$compatible==TRUE,] ## select splicing compatible transcript matches

  minReadClassTable <- as_tibble(spliceOverlapsSelected) %>%
    dplyr::select(queryHits, subjectHits)
  minReadClassTable$queryTxId <- names(exByTxAnnotated_singleBpStartEnd)[minReadClassTable$queryHits]
  minReadClassTable$subjectTxId <- names(exByTxAnnotated_singleBpStartEnd)[minReadClassTable$subjectHits]
  minReadClassTable <- minReadClassTable %>%
    group_by(queryTxId) %>%
    arrange(queryTxId, subjectTxId) %>%
    mutate(eqClass = paste(subjectTxId, collapse='.'), minEqClassSize = n()) %>%
    dplyr::select(queryTxId, eqClass, minEqClassSize) %>%
    distinct()
  return(minReadClassTable)
}



calculateDistToAnnotation <- function(exByTx, exByTxRef, maxDist = 35, primarySecondaryDist = 5, ignore.strand=FALSE) {  #### HERE #### 2019-10-03: create table with read calss to transcript and read class to gene annotation. ANything iwhtout match is either new (annotated as new 5' exon, new 3' exon, new exon withmin 30bp threhsold(? etc) or discarded (?) (classificaiton has medium results... only in addition?)

  #(1)  find overlaps of read classes with annotated transcripts, allow for maxDist [b] distance for each exon; exons with size less than 35bp are dropped to find overlaps, but counted towards distance and compatibility
  spliceOverlaps <- findSpliceOverlapsByDist(exByTx,exByTxRef, maxDist = maxDist, firstLastSeparate = T, dropRangesByMinLength = T,cutStartEnd=TRUE, ignore.strand=ignore.strand)
  txToAnTable <- tbl_df(spliceOverlaps)%>%  group_by(queryHits)  %>%
    mutate(dist = uniqueLengthQuery + uniqueLengthSubject) %>%
    mutate(txNumber = n())

  # first round of filtering should only exclude obvious mismatches
  txToAnTableFiltered <- txToAnTable %>%
    group_by(queryHits)  %>%
    arrange(queryHits, dist) %>%
    filter(dist<=(min(dist)+primarySecondaryDist)) %>%
    filter(queryElementsOutsideMaxDist+subjectElementsOutsideMaxDist == min(queryElementsOutsideMaxDist+subjectElementsOutsideMaxDist )) %>%
    filter((uniqueStartLengthQuery<=primarySecondaryDist & uniqueEndLengthQuery<=primarySecondaryDist) == max(uniqueStartLengthQuery<=primarySecondaryDist & uniqueEndLengthQuery<=primarySecondaryDist)) %>%
    mutate(txNumberFiltered = n())


  # (2) calculate splice overlap for any not in the list (all hits have a unique new exon of at least 35bp length, might be new candidates)
  setTMP = unique(txToAnTableFiltered$queryHits)
  spliceOverlaps_rest <- findSpliceOverlapsByDist(exByTx[-setTMP],exByTxRef, maxDist = 0, type='any',  firstLastSeparate = T, dropRangesByMinLength=F, cutStartEnd=TRUE, ignore.strand=ignore.strand)
  txToAnTableRest <- tbl_df(spliceOverlaps_rest)%>%  group_by(queryHits)%>%
    mutate(dist = uniqueLengthQuery + uniqueLengthSubject) %>%
    mutate(txNumber = n())
  txToAnTableRest$queryHits <- (1:length(exByTx))[-setTMP][txToAnTableRest$queryHits]  # reassign IDs based on unfiltered list length

  # todo: check filters, what happens to reads with only start and end match?
  txToAnTableRest <- txToAnTableRest %>%
    group_by(queryHits)  %>%
    arrange(queryHits, dist) %>%
    filter(dist<=(min(dist)+primarySecondaryDist)) %>%
    filter(queryElementsOutsideMaxDist+subjectElementsOutsideMaxDist == min(queryElementsOutsideMaxDist+subjectElementsOutsideMaxDist )) %>%
    filter((uniqueStartLengthQuery<=primarySecondaryDist & uniqueEndLengthQuery<=primarySecondaryDist) == max(uniqueStartLengthQuery<=primarySecondaryDist & uniqueEndLengthQuery<=primarySecondaryDist)) %>%
    mutate(txNumberFiltered = n())

  # (3) find overlaps for remaining reads (reads which have start/end match, this time not cut and used to calculate distance)
  setTMPRest = unique(c(txToAnTableRest$queryHits,setTMP))
  txToAnTableRestStartEnd <- NULL
  if(length(exByTx[-setTMPRest])>0) {
    spliceOverlaps_restStartEnd <- findSpliceOverlapsByDist(exByTx[-setTMPRest],exByTxRef, maxDist = 0, type='any',  firstLastSeparate = T, dropRangesByMinLength=F, cutStartEnd=F, ignore.strand=ignore.strand)
    txToAnTableRestStartEnd <- tbl_df(spliceOverlaps_restStartEnd)%>%  group_by(queryHits) %>%
      mutate(dist = uniqueLengthQuery + uniqueLengthSubject + uniqueStartLengthQuery + uniqueEndLengthQuery) %>%
      mutate(txNumber = n())
    txToAnTableRestStartEnd$queryHits <- (1:length(exByTx))[-setTMPRest][txToAnTableRestStartEnd$queryHits]  # reassign IDs based on unfiltered list length

    # todo: check filters, what happens to reads with only start and end match?
    txToAnTableRestStartEnd <- txToAnTableRestStartEnd %>%
      group_by(queryHits)  %>%
      arrange(queryHits, dist) %>%
      filter(dist<=(min(dist)+primarySecondaryDist)) %>%
      mutate(txNumberFiltered = n())

  }

  txToAnTableFiltered <- rbind(txToAnTableFiltered, txToAnTableRest, txToAnTableRestStartEnd) %>% ungroup()

  txToAnTableFiltered$readClassId = names(exByTx)[txToAnTableFiltered$queryHits]
  txToAnTableFiltered$annotationTxId =names( exByTxRef)[txToAnTableFiltered$subjectHits]

  #  txToAnTableFiltered <- left_join(txToAnTableFiltered, select(txList$txTable, txId, readCount, confidenceType)) %>% mutate(relativeReadCount=readCount/txNumberFiltered)

  return(txToAnTableFiltered)


  ###

  ####       HERE WORK ON THIS TODO #################### 2019-10-07
  ## Todo: (1) add read classes that are not within the annotation (check done)
  ## (2) add single exon read classes
  ## (2.5) add more filters, i.e. low confidence reads vs high confidence, then add some to new annotations, don't assign others even if single hit (maybe distinct transcript, in that case only count to gene expression!)
  ## (3) clean up code
  ## write function to add exon rank and exon end rank for grangeslist !!!

  #                   tmpData=(test[which(test$subjectHits %in% which(txdbTablesList$txIdToGeneIdTable$GENEID=='ENSG00000132781')),])
  ## here: overlap with type='any', then implement further filter rules see below (2019-10-06) ###
  ##### HERE #####
  ## filter rules:
  ## (1) select minimum distance match (note: allow for a few base pairs error?)
  ## (2) select hits within (minimum unique query sequence, no start match)
  ## (3) select hits with minimum unqiue start sequence/end sequence

  ## then add matches with new exons (i.e. query not within subject, currently all removed, or not properly calculted?) add classificaiton into class (ie. new exon, new start, new end,...)

  ######################## HERE !############## create tx to annoation table
  #### somehow use start and end match information? because ic ut the last and first exon, these don't count to the distnace calculation. Maybe cut only within distance calucation to findoverlaps, bu calculate distance in addition for the first and last exon as extra information? ### HERE

  #     spliceOverlaps <- findSpliceOverlapsByDist(txExonsByCut,annotationsExonsByCut, maxDist = maxDist)
  #
  #     txToTranscriptTable <- tbl_df((spliceOverlaps[mcols(spliceOverlaps)$unique==TRUE]))
  #     so_2 <- spliceOverlaps_0bp[! (queryHits(spliceOverlaps_0bp) %in% queryHits(spliceOverlaps_0bp[mcols(spliceOverlaps_0bp)$unique==TRUE]))]
  #     spliceTableTMP <- tbl_df(so_2[mcols(so_2)$compatible==TRUE])
  #
  #
  #     readClassToTranscriptTable <- txList$txTable %>% select(txId)
  #     readClassToGeneTable <- txList$txTable %>% select(txId, geneId)

}










#=====================================================================
#                      development functions                         =
#=====================================================================
#return predictions for transcripts how likely they are true (annotated)
predictTrueTranscriptsDirectRna <- function(txDataTable, interactions = F, model = NULL, verbose = T, maxTestData=5000, cores=1) {

  require(doMC)
  registerDoMC(cores=cores)


  txDataTable <- txDataTable %>% dplyr::select(combinedTxId,exonOverlapGeneId,nAnnotatedTxByGene,confidenceType,firstJunctionCountPerRead, lastJunctionCountPerRead, firstJunctionInternalCountPerRead, lastJunctionInternalCountPerRead,readCount, relativeReadCount, uniqueReadCount, fullSpliceSupportReadCount, txSubsetNumber, txStartOverlapNumber, txEndOverlapNumber, uniqueReadCountByGene, fullSpliceSupportReadCountByGene, nTxByGene) %>%
    mutate( rel_firstJunctionCountPerRead = firstJunctionCountPerRead/readCount,
            rel_lastJunctionCountPerRead = lastJunctionCountPerRead/readCount,
            rel_firstJunctionInternalCountPerRead = firstJunctionInternalCountPerRead/readCount,
            rel_lastJunctionInternalCountPerRead = lastJunctionInternalCountPerRead/readCount,
            rel_readCount = readCount/uniqueReadCountByGene,
            txSubsetNumberBinary = (txSubsetNumber==0),
            rel_firstJunctionCountPerRead_log = pmax(log(rel_firstJunctionCountPerRead),0),
            rel_lastJunctionCountPerRead_log = pmax(log(rel_lastJunctionCountPerRead),0),
            rel_firstJunctionInternalCountPerRead_log = pmax(log(rel_firstJunctionInternalCountPerRead),-5),
            rel_lastJunctionInternalCountPerRead_log = pmax(log(rel_lastJunctionInternalCountPerRead),-5),
            rel_readCount_log=log(rel_readCount+0.001),
            fullSpliceSupportReadCount_log= log(fullSpliceSupportReadCount+1)
    )

  #                          txDataTableTest <- txDataTableTest %>% select(combinedTxId,exonOverlapGeneId,nAnnotatedTxByGene,confidenceType,firstJunctionCountPerRead, lastJunctionCountPerRead, firstJunctionInternalCountPerRead, lastJunctionInternalCountPerRead,readCount, relativeReadCount, uniqueReadCount, fullSpliceSupportReadCount, txSubsetNumber, txStartOverlapNumber, txEndOverlapNumber, uniqueReadCountByGene, fullSpliceSupportReadCountByGene, nTxByGene) %>%
  #     mutate( rel_firstJunctionCountPerRead = firstJunctionCountPerRead/readCount,
  #                                 rel_lastJunctionCountPerRead = lastJunctionCountPerRead/readCount,
  #                                 rel_firstJunctionInternalCountPerRead = firstJunctionInternalCountPerRead/readCount,
  #                                 rel_lastJunctionInternalCountPerRead = lastJunctionInternalCountPerRead/readCount,
  #                                 rel_readCount = readCount/uniqueReadCountByGene,
  #                                 txSubsetNumberBinary = (txSubsetNumber==0),
  #                        rel_firstJunctionCountPerRead_log = pmax(log(rel_firstJunctionCountPerRead),0),
  #                        rel_lastJunctionCountPerRead_log = pmax(log(rel_lastJunctionCountPerRead),0),
  #                        rel_firstJunctionInternalCountPerRead_log = pmax(log(rel_firstJunctionInternalCountPerRead),-5),
  #                        rel_lastJunctionInternalCountPerRead_log = pmax(log(rel_lastJunctionInternalCountPerRead),-5),
  #                         rel_readCount_log=log(rel_readCount+0.001),
  #                         fullSpliceSupportReadCount_log= log(fullSpliceSupportReadCount)
  #                         )
  #

  txDataTableTrain <-txDataTable %>% filter(nAnnotatedTxByGene>=0, readCount>1, uniqueReadCountByGene > 3 , fullSpliceSupportReadCountByGene>3, confidenceType == "highConfidenceJunctionReads", txSubsetNumber>=0)
  trainingLabels <- !grepl('^tx',txDataTableTrain$combinedTxId)

  txDataTableTrain <-txDataTable # %>% filter(nAnnotatedTxByGene>=0, readCount>1, uniqueReadCountByGene > 3 , fullSpliceSupportReadCountByGene>3, confidenceType == "highConfidenceJunctionReads", txSubsetNumber>=0)
  trainingLabels <- !grepl('^tx',txDataTableTrain$combinedTxId)


  modelData <- txDataTableTrain %>%
    dplyr::select(firstJunctionCountPerRead,
           lastJunctionCountPerRead,
           firstJunctionInternalCountPerRead,
           lastJunctionInternalCountPerRead,
           readCount,
           relativeReadCount,
           uniqueReadCount,
           fullSpliceSupportReadCount,
           txSubsetNumber,
           txStartOverlapNumber,
           txEndOverlapNumber,
           uniqueReadCountByGene,
           fullSpliceSupportReadCountByGene,
           nTxByGene,
           rel_readCount,
           rel_firstJunctionCountPerRead,
           rel_lastJunctionCountPerRead,
           rel_firstJunctionInternalCountPerRead,
           rel_lastJunctionInternalCountPerRead,
           rel_firstJunctionCountPerRead_log,
           rel_lastJunctionCountPerRead_log,
           rel_firstJunctionInternalCountPerRead_log,
           rel_lastJunctionInternalCountPerRead_log,
           rel_readCount_log,
           fullSpliceSupportReadCount_log,
           txSubsetNumberBinary)
  modelData[is.na(modelData)] <- 0

  testData <- txDataTable %>%
    dplyr::select(firstJunctionCountPerRead,
           lastJunctionCountPerRead,
           firstJunctionInternalCountPerRead,
           lastJunctionInternalCountPerRead,
           readCount,
           relativeReadCount,
           uniqueReadCount,
           fullSpliceSupportReadCount,
           txSubsetNumber,
           txStartOverlapNumber,
           txEndOverlapNumber,
           uniqueReadCountByGene,
           fullSpliceSupportReadCountByGene,
           nTxByGene,
           rel_readCount,
           rel_firstJunctionCountPerRead,
           rel_lastJunctionCountPerRead,
           rel_firstJunctionInternalCountPerRead,
           rel_lastJunctionInternalCountPerRead,
           rel_firstJunctionCountPerRead_log,
           rel_lastJunctionCountPerRead_log,
           rel_firstJunctionInternalCountPerRead_log,
           rel_lastJunctionInternalCountPerRead_log,
           rel_readCount_log,
           fullSpliceSupportReadCount_log,
           txSubsetNumberBinary)

  testData[is.na(testData)] <- 0

  mySample=sample(1:nrow(modelData),10000)

  f <- as.formula(y ~ .*.)

  y <- trainingLabels[mySample]
  x <- model.matrix(f, modelData[mySample,])[, -1]
  x_test <- model.matrix(~.*., testData)[, -1]
  # x <- x_test[mySample,]
  #  test=glmnet(x, y)
  # predictions= predict(test,newx=x_test)
  cv.fit=glmnet::cv.glmnet(x=x,y=y, parallel = T)
  ##  cv.fit=cv.glmnet(x=x,y=y,family='binomial', parallel = T)  ## slow but more correct


  f <- as.formula(y ~ .)

  y <- trainingLabels[mySample]
  x <- model.matrix(~., modelData)[, -1]
  x_test <- model.matrix(~., testData)[, -1]
  x <- x_test[mySample,]
  #  test=glmnet(x, y)
  # predictions= predict(test,newx=x_test)
  cv.fit=glmnet::cv.glmnet(x=x,y=y,family='binomial', parallel = T)


  predictions= predict(cv.fit,newx=x_test,s='lambda.min')

  j=1
  plot(c(0,1),c(0,1),xlim=c(0,1),ylim=c(0,1),ty='l')
  for(i in unique(txDataTable$confidenceType)){
    show(i)
    performance_tmp=myPerformance(grepl('^E',txDataTable$combinedTxId)[txDataTable$confidenceType==i],predictions[txDataTable$confidenceType==i])
    show( performance_tmp$AUC)
    lines(performance_tmp$FPR, performance_tmp$TPR,col=j)
    lines(performance_tmp$TPR, performance_tmp$precision,col=j)
    j=j+1
  }



}


#return predictions for transcripts how likely they are true (annotated)
evaluateTxdiscriminationModels <- function(txDataTable, interactions = F, model = NULL, verbose = T, maxTestData=5000, cores=1) {
  require(doMC)
  registerDoMC(cores=cores)


  txDataTable <- txDataTable %>% dplyr::select(combinedTxId,exonOverlapGeneId,nAnnotatedTxByGene,confidenceType,firstJunctionCountPerRead, lastJunctionCountPerRead, firstJunctionInternalCountPerRead, lastJunctionInternalCountPerRead,readCount, relativeReadCount, uniqueReadCount, fullSpliceSupportReadCount, txSubsetNumber, txStartOverlapNumber, txEndOverlapNumber, uniqueReadCountByGene, fullSpliceSupportReadCountByGene, nTxByGene) %>%
    mutate( rel_firstJunctionCountPerRead = firstJunctionCountPerRead/readCount,
            rel_lastJunctionCountPerRead = lastJunctionCountPerRead/readCount,
            rel_firstJunctionInternalCountPerRead = firstJunctionInternalCountPerRead/readCount,
            rel_lastJunctionInternalCountPerRead = lastJunctionInternalCountPerRead/readCount,
            rel_readCount = readCount/uniqueReadCountByGene,
            txSubsetNumberBinary = (txSubsetNumber==0),
            rel_firstJunctionCountPerRead_log = pmax(log(rel_firstJunctionCountPerRead),0),
            rel_lastJunctionCountPerRead_log = pmax(log(rel_lastJunctionCountPerRead),0),
            rel_firstJunctionInternalCountPerRead_log = pmax(log(rel_firstJunctionInternalCountPerRead),-5),
            rel_lastJunctionInternalCountPerRead_log = pmax(log(rel_lastJunctionInternalCountPerRead),-5),
            rel_readCount_log=log(rel_readCount+0.001),
            fullSpliceSupportReadCount_log= log(fullSpliceSupportReadCount)
    )

  resultMatrix <- matrix(0, ncol=11, nrow=10)
  colnames(resultMatrix) <- c('trainMode','testMode',paste('model', 1:(ncol(resultMatrix)-2), sep='_'))



  trainingTxData <-txDataTable %>% filter(nAnnotatedTxByGene>=0, readCount>1, uniqueReadCountByGene > 3 , fullSpliceSupportReadCountByGene>3, confidenceType == "highConfidenceJunctionReads", txSubsetNumber>=0)
  trainingLabels <- !grepl('^tx',trainingTxData$combinedTxId)


  modelData <- list()

  # model 1: baseline
  modelData[[1]] <- trainingTxData %>%
    dplyr::select(readCount,
           uniqueReadCountByGene)

  # # model2: count based
  #   modelData[[2]] <- trainingTxData %>%
  #                     select(firstJunctionCountPerRead,
  #                         lastJunctionCountPerRead,
  #                         firstJunctionInternalCountPerRead,
  #                         lastJunctionInternalCountPerRead,
  #                         readCount,
  #                         relativeReadCount,
  #                         uniqueReadCount,
  #                         fullSpliceSupportReadCount,
  #                         txSubsetNumber,
  #                         txStartOverlapNumber,
  #                         txEndOverlapNumber,
  #                         uniqueReadCountByGene,
  #                         fullSpliceSupportReadCountByGene,
  #                         nTxByGene,
  #                         rel_readCount)
  # model2: relative count based (normalised by gene count)
  modelData[[2]] <- trainingTxData %>%
    dplyr::select(firstJunctionCountPerRead,
           lastJunctionCountPerRead,
           firstJunctionInternalCountPerRead,
           lastJunctionInternalCountPerRead,
           readCount,
           relativeReadCount,
           uniqueReadCount,
           fullSpliceSupportReadCount,
           txSubsetNumber,
           txStartOverlapNumber,
           txEndOverlapNumber,
           uniqueReadCountByGene,
           fullSpliceSupportReadCountByGene,
           nTxByGene,
           rel_readCount,
           rel_firstJunctionCountPerRead,
           rel_lastJunctionCountPerRead,
           rel_firstJunctionInternalCountPerRead,
           rel_lastJunctionInternalCountPerRead,
           rel_readCount) %>%
    mutate(txSubsetNumber = (txSubsetNumber==0))

  # # model4: relative count based (normalised by gene count), minimal variables
  #    modelData[[4]] <- trainingTxData %>%
  #                     select(readCount,
  #                         uniqueReadCountByGene,
  #                         fullSpliceSupportReadCountByGene,
  #                         nTxByGene,
  #                         rel_readCount,
  #                         rel_firstJunctionCountPerRead,
  #                        rel_lastJunctionCountPerRead,
  #                        rel_firstJunctionInternalCountPerRead,
  #                        rel_lastJunctionInternalCountPerRead,
  #                        rel_readCount)

  #    #model5: relative read count + relative internal read counts, log transformd
  #
  #  modelData[[5]] <- trainingTxData %>%
  #                 select(rel_firstJunctionCountPerRead_log,
  #                        rel_lastJunctionCountPerRead_log,
  #                        rel_firstJunctionInternalCountPerRead_log,
  #                        rel_lastJunctionInternalCountPerRead_log,
  #                        rel_readCount_log,
  #                        txSubsetNumberBinary, fullSpliceSupportReadCount_log,nTxByGene, txStartOverlapNumber, txEndOverlapNumber)
  #
  #

  #model3 maximum model log based


  modelData[[3]] <- trainingTxData %>%
    dplyr::select(firstJunctionCountPerRead,
           lastJunctionCountPerRead,
           firstJunctionInternalCountPerRead,
           lastJunctionInternalCountPerRead,
           readCount,
           relativeReadCount,
           uniqueReadCount,
           fullSpliceSupportReadCount,
           txSubsetNumber,
           txStartOverlapNumber,
           txEndOverlapNumber,
           uniqueReadCountByGene,
           fullSpliceSupportReadCountByGene,
           nTxByGene,
           rel_readCount,
           rel_firstJunctionCountPerRead,
           rel_lastJunctionCountPerRead,
           rel_firstJunctionInternalCountPerRead,
           rel_lastJunctionInternalCountPerRead,
           rel_firstJunctionCountPerRead_log,
           rel_lastJunctionCountPerRead_log,
           rel_firstJunctionInternalCountPerRead_log,
           rel_lastJunctionInternalCountPerRead_log,
           rel_readCount_log,
           fullSpliceSupportReadCount_log,
           txSubsetNumberBinary)

  #model4 optimal model for 2 classes separate ?


  modelData[[4]] <- trainingTxData %>%
    dplyr::select(
      uniqueReadCount,
      fullSpliceSupportReadCount,
      uniqueReadCountByGene,
      rel_firstJunctionCountPerRead,
      rel_lastJunctionCountPerRead,
      rel_firstJunctionInternalCountPerRead,
      rel_lastJunctionInternalCountPerRead,
      rel_firstJunctionCountPerRead_log,
      rel_lastJunctionCountPerRead_log,
      rel_firstJunctionInternalCountPerRead_log,
      rel_lastJunctionInternalCountPerRead_log,
      rel_readCount_log,
      fullSpliceSupportReadCount_log)



  #################HERE @################

  n_train <- pmin(ceiling(length(unique(trainingTxData$exonOverlapGeneId)) * 0.15), maxTestData)
  n_test <- pmin(length(unique(trainingTxData$exonOverlapGeneId)) - n_train, maxTestData)




  myGeneSample <- sample(unique(trainingTxData$exonOverlapGeneId),n_train+n_test, replace=F)

  myIndex = 1
  for(myScenarioTrain in c('allTx','uniqueTx','subsetTx')) {

    if(myScenarioTrain=='allTx') {
      mySampleTrain <- trainingTxData$exonOverlapGeneId %in% myGeneSample[1:n_train]
    } else if(myScenarioTrain=='uniqueTx') {
      mySampleTrain <- trainingTxData$txSubsetNumber==0 & trainingTxData$exonOverlapGeneId %in% myGeneSample[1:n_train]
    } else if(myScenarioTrain=='subsetTx') {
      mySampleTrain <- trainingTxData$txSubsetNumber>0 & trainingTxData$exonOverlapGeneId %in% myGeneSample[1:n_train]
    }
    for(myScenarioTest in c('allTx','uniqueTx','subsetTx')) {
      if(myScenarioTest=='allTx') {
        mySampleTest <- trainingTxData$exonOverlapGeneId %in% myGeneSample[-(1:n_train)]
      } else if(myScenarioTest=='uniqueTx') {
        mySampleTest = trainingTxData$txSubsetNumber==0 & (trainingTxData$exonOverlapGeneId %in% myGeneSample[-(1:n_train)])
      } else if(myScenarioTest=='subsetTx') {
        mySampleTest = trainingTxData$txSubsetNumber>0 & (trainingTxData$exonOverlapGeneId %in% myGeneSample[-(1:n_train)])
      }
      resultMatrix[myIndex,'trainMode'] <- myScenarioTrain
      resultMatrix[myIndex,'testMode'] <- myScenarioTest
      show(myScenarioTrain)
      show(myScenarioTest)


      plot(c(0,1),c(0,1),xlim=c(0,1),ylim=c(0,1),ty='l')
      for(j in 1:(length(modelData))) {
        # test_model1 <- fitBinomialModel(trainingLabels[mySampleTrain], as.matrix(modelData[[j]][mySampleTrain,]), as.matrix(modelData[[j]][mySampleTest,]), show.cv=TRUE, maxSize.cv=10000, parallel = T)
        cv.fit=glmnet::cv.glmnet(x=as.matrix(modelData[[j]][mySampleTrain,]),y=trainingLabels[mySampleTrain],family='binomial',  parallel = T)
        predictions= predict(cv.fit,newx=as.matrix(modelData[[j]][mySampleTest,]),s='lambda.min')
        performance_tmp=myPerformance(trainingLabels[mySampleTest]==1,predictions)
        resultMatrix[myIndex,j+2] <- performance_tmp$AUC
        #
        # show(performance_tmp$AUC)
        # show( fisher.test(table(predictions>0,trainingLabels[mySampleTest])))
        # lines(performance_tmp$FPR, performance_tmp$TPR,col=j)
        # lines(performance_tmp$TPR, performance_tmp$precision,col=j)
      }
      j=j+1

      f <- as.formula(y ~ .*.)
      y <- trainingLabels[mySampleTrain]
      x <- model.matrix(f, modelData[[3]][mySampleTrain,])[, -1]
      x_test <- model.matrix(~.*., modelData[[3]][mySampleTest,])[, -1]
      #  test=glmnet(x, y)
      # predictions= predict(test,newx=x_test)
      cv.fit=glmnet::cv.glmnet(x=x,y=y, parallel = T)
      ##  cv.fit=cv.glmnet(x=x,y=y,family='binomial', parallel = T)  ## slow but more correct
      predictions= predict(cv.fit,newx=x_test,s='lambda.min')
      performance_tmp=myPerformance(trainingLabels[mySampleTest]==1,predictions)

      resultMatrix[myIndex,j+2] <- performance_tmp$AUC
      #
      # model.ksvm <- ksvm(x, y, kernel='vanilla', type="C-svc")  # build SVM model
      # show(model.ksvm)
      # predictions <- predict(model.ksvm,x_test, type='decision')  # predict labels for test data using the SVM model
      # performance_tmp=myPerformance(trainingLabels[mySampleTest]==1,predictions)

      #
      # show( fisher.test(table(predictions>0,trainingLabels[mySampleTest])))
      # show(performance_tmp$AUC)
      # lines(performance_tmp$FPR, performance_tmp$TPR,col=j)
      # lines(performance_tmp$TPR, performance_tmp$precision,col=j)



      myIndex=myIndex+1
      show(resultMatrix)
    }
  }

  return(resultMatrix)
}


#development code for a function to use barcode in cDNA data to identify full length transcripts
useBarcodeData <- function(combinedTxTable, readTxTable, barcodeData) {
  data=readRDS(file='/mnt/ont/pychopper/poreDecode-bam/GIS_Hct116_cDNA_Rep2-Run3_nanoporeDecode.rds')

  candidateGene='R2_54'
  candidateGene='R1_71'
  #candidateGene='ENSG00000126698'
  #candidateGene='ENSG00000167996'

  combinedTxTableSubset <- as.data.frame(combinedTxTable[which(combinedTxTable$exonOverlapGeneId==candidateGene & combinedTxTable$readCount>5),])
  combinedExonRangesExtendedSubset <- combinedExonRangesExtended[combinedTxTableSubset$combinedTxId]

  readTxTableSubset=readTxTable[(readTxTable$combinedTxId %in% combinedTxTable$combinedTxId[which(combinedTxTable$exonOverlapGeneId==candidateGene)]),]

  ##readTxTableSubset$originalStrand <- as.character(strand(reads[readTxTableSubset$readName]))

  # dataSubset <- as_tibble(data[data$names %in% readTxTableSubset$readName,]) %>% select(names, predictedStrand, strandScore, pa3Prime, pa3PrimeMinus)
  # test=left_join(readTxTableSubset, dataSubset, c('readName'='names'))
  # test$polyATail <- test$pa3Prime
  # test$polyATail[test$originalStrand=='-'] <-  test$pa3PrimeMinus[test$originalStrand=='-']
  # test$strandScoreTranscription <- test$strandScore
  # test$strandScoreTranscription[test$originalStrand=='-'] <-  (-test$strandScoreTranscription[test$originalStrand=='-'])



  # plotTracks(list(makeTrackFromGrangesList(combinedExonRangesExtendedSubset),makeTrackFromGrangesList(readGrglist[sample(readTxTableSubset$readName,50)])),showId=TRUE)
  # plotTracks(list(makeTrackFromGrangesList(combinedExonRangesExtendedSubset), makeTrackFromGrangesList(combinedExonRangesExtended[queryHits(overlapsCandidates)]),makeTrackFromGrangesList(readGrglist[readTxTableSubset$readName])),showId=TRUE)


  # ML model to predict annotated (true transcripts)


  ### todo: build multiple models (1 for uniquee TX, one for subset tx, use data from barcode search to identify true transcripts


  #dataSubset <- as_tibble(data[data$names %in% readTxTable$readName,]) %>% select(names, strandScore, pa3Prime, pa3PrimeMinus, barcode5Prime.plus.score, barcode5Prime.minus.score, barcode3Prime.plus.score, barcode3Prime.minus.score)

  # test=left_join(readTxTable, data, c('readName'='names'))
  # test$polyATail <- test$pa3Prime
  # test$polyATail[test$originalStrand=='-'] <-  test$pa3PrimeMinus[test$originalStrand=='-']
  # test$strandScoreTranscription <- test$strandScore
  # test$strandScoreTranscription[test$originalStrand=='-'] <-  (-test$strandScoreTranscription[test$originalStrand=='-'])

  ## how many training data points
  ## evaluation on unique Tx vs non-unique/subset Tx prediction
  ## evaluation on sequin data
  ## which features

  trainingTxData <- combinedTxTable %>% filter(nAnnotatedTxByGene>=0, readCount>1, uniqueReadCountByGene > 3 , fullSpliceSupportReadCountByGene>3, confidenceType == "highConfidenceJunctionReads", txSubsetNumber>=0)

  dataSubset <- as_tibble(data[data$names %in% readTxTable$readName[readTxTable$confidenceType=='highConfidenceJunctionReads' & readTxTable$combinedTxId %in% trainingTxData$combinedTxId],]) %>%
    dplyr::select(names, strandScore, pa3Prime, pa3PrimeMinus, barcode5Prime.plus.score, barcode5Prime.minus.score, barcode3Prime.plus.score, barcode3Prime.minus.score)

  readData=left_join(readTxTable[readTxTable$confidenceType=='highConfidenceJunctionReads' & readTxTable$combinedTxId %in% trainingTxData$combinedTxId,], dataSubset, c('readName'='names'))

  readData$polyATail <- readData$pa3Prime
  readData$polyATail[readData$originalStrand!=readData$strand] <-  readData$pa3PrimeMinus[readData$originalStrand!=readData$strand]
  readData$strandScoreTranscription <- readData$strandScore
  readData$strandScoreTranscription[readData$originalStrand!=readData$strand] <-  (-readData$strandScoreTranscription[readData$originalStrand!=readData$strand])

  readData$barcode5PrimeScore <- readData$barcode5Prime.plus.score
  readData$barcode5PrimeScore[readData$originalStrand!=readData$strand] <-  readData$barcode5Prime.minus.score[readData$originalStrand!=readData$strand]
  readData$barcode3PrimeScore <- readData$barcode3Prime.plus.score
  readData$barcode3PrimeScore[readData$originalStrand!=readData$strand] <-  readData$barcode3Prime.minus.score[readData$originalStrand!=readData$strand]

  readData2 <- readData %>% group_by(combinedTxId) %>%
    mutate(meanStrandScore = mean(strandScoreTranscription),
           meanPolyATail = mean(polyATail),
           meanBarcode5PrimeScore = mean(barcode5PrimeScore),
           meanBarcode3PrimeScore = mean(barcode3PrimeScore),
           relCountPolyA5 = sum(polyATail>5)/n(),
           relCountBarcode5Prime = sum(barcode5PrimeScore>0)/n(),
           relCountBarcode3Prime = sum(barcode3PrimeScore>0)/n()) %>%
    dplyr::select(combinedTxId, strand, meanStrandScore, meanPolyATail, meanBarcode5PrimeScore, meanBarcode3PrimeScore, relCountPolyA5, relCountBarcode5Prime, relCountBarcode3Prime) %>%
    distinct() %>% ungroup()

  trainingTxData <- left_join(trainingTxData, readData2, by='combinedTxId')

  trainingTxDataSelect <- trainingTxData %>% dplyr::select(firstJunctionCountPerRead, lastJunctionCountPerRead, firstJunctionInternalCountPerRead, lastJunctionInternalCountPerRead,readCount, relativeReadCount, uniqueReadCount, fullSpliceSupportReadCount, txSubsetNumber, txStartOverlapNumber, txEndOverlapNumber, uniqueReadCountByGene, fullSpliceSupportReadCountByGene, nTxByGene,
                                                    meanStrandScore,meanPolyATail,meanBarcode5PrimeScore,meanBarcode3PrimeScore,relCountPolyA5,relCountBarcode5Prime,relCountBarcode3Prime)

  # trainingTxDataSelect <- trainingTxData %>% select(firstJunctionCountPerRead, lastJunctionCountPerRead, firstJunctionInternalCountPerRead, lastJunctionInternalCountPerRead,readCount, relativeReadCount, uniqueReadCount, fullSpliceSupportReadCount, txSubsetNumber, txStartOverlapNumber, txEndOverlapNumber, uniqueReadCountByGene, fullSpliceSupportReadCountByGene, nTxByGene)

  trainingTxDataSelect <- trainingTxDataSelect %>%
    mutate( rel_firstJunctionCountPerRead = firstJunctionCountPerRead/readCount,
            rel_lastJunctionCountPerRead = lastJunctionCountPerRead/readCount,
            rel_firstJunctionInternalCountPerRead = firstJunctionInternalCountPerRead/readCount,
            rel_lastJunctionInternalCountPerRead = lastJunctionInternalCountPerRead/readCount,
            rel_readCount = readCount/uniqueReadCountByGene)

  trainingLabels <- !grepl('^tx',trainingTxData$combinedTxId)

  model1_data <- trainingTxDataSelect %>%
    mutate(txSubsetNumber = (txSubsetNumber==0))

  model2_data <- trainingTxDataSelect %>%
    dplyr::select(rel_firstJunctionCountPerRead,
           rel_lastJunctionCountPerRead,
           rel_firstJunctionInternalCountPerRead,
           rel_lastJunctionInternalCountPerRead,
           rel_readCount)
  #direct RNA
  #  model3_data <- trainingTxDataSelect %>%
  #                 select(rel_firstJunctionCountPerRead,
  #                        rel_lastJunctionCountPerRead,
  #                        rel_firstJunctionInternalCountPerRead,
  #                        rel_lastJunctionInternalCountPerRead,
  #                        rel_readCount,
  #                        txSubsetNumber, fullSpliceSupportReadCount,nTxByGene, fullSpliceSupportReadCount, txStartOverlapNumber, txEndOverlapNumber) %>%
  #                        mutate(txSubsetNumber = (txSubsetNumber==0))
  #
  model3_data <- trainingTxDataSelect %>%
    dplyr::select(rel_firstJunctionCountPerRead,
           rel_lastJunctionCountPerRead,
           rel_firstJunctionInternalCountPerRead,
           rel_lastJunctionInternalCountPerRead,
           rel_readCount,
           txSubsetNumber, fullSpliceSupportReadCount,nTxByGene, fullSpliceSupportReadCount, txStartOverlapNumber, txEndOverlapNumber) %>%
    mutate(txSubsetNumber = (txSubsetNumber==0),
           rel_firstJunctionCountPerRead = pmax(log(rel_firstJunctionCountPerRead),0),
           rel_lastJunctionCountPerRead = pmax(log(rel_lastJunctionCountPerRead),0),
           rel_firstJunctionInternalCountPerRead = pmax(log(rel_firstJunctionInternalCountPerRead),-5),
           rel_lastJunctionInternalCountPerRead = pmax(log(rel_lastJunctionInternalCountPerRead),-5),
           rel_readCount=log(rel_readCount+0.001),
           fullSpliceSupportReadCount= log(fullSpliceSupportReadCount))

  model3_data <- trainingTxDataSelect %>%
    dplyr::select(rel_firstJunctionCountPerRead,
           rel_lastJunctionCountPerRead,
           rel_firstJunctionInternalCountPerRead,
           rel_lastJunctionInternalCountPerRead,
           rel_readCount,
           txSubsetNumber, meanStrandScore,meanPolyATail,meanBarcode5PrimeScore,meanBarcode3PrimeScore,relCountPolyA5,relCountBarcode5Prime,relCountBarcode3Prime,
           fullSpliceSupportReadCount,nTxByGene, fullSpliceSupportReadCount, txStartOverlapNumber, txEndOverlapNumber) %>%
    mutate(txSubsetNumber = (txSubsetNumber==0))

  model3_data <- trainingTxDataSelect %>%
    dplyr::select(rel_firstJunctionCountPerRead,
           rel_lastJunctionCountPerRead,
           rel_firstJunctionInternalCountPerRead,
           rel_lastJunctionInternalCountPerRead,
           rel_readCount,
           txSubsetNumber, meanStrandScore,meanPolyATail,meanBarcode5PrimeScore,meanBarcode3PrimeScore,
           fullSpliceSupportReadCount,nTxByGene, fullSpliceSupportReadCount, txStartOverlapNumber, txEndOverlapNumber) %>%
    mutate(txSubsetNumber = (txSubsetNumber==0),
           rel_firstJunctionCountPerRead = pmax(log(rel_firstJunctionCountPerRead),0),
           rel_lastJunctionCountPerRead = pmax(log(rel_lastJunctionCountPerRead),0),
           rel_firstJunctionInternalCountPerRead = pmax(log(rel_firstJunctionInternalCountPerRead),-5),
           rel_lastJunctionInternalCountPerRead = pmax(log(rel_lastJunctionInternalCountPerRead),-5),
           rel_readCount=log(rel_readCount+0.001),
           fullSpliceSupportReadCount= log(fullSpliceSupportReadCount))


  #mySample <- sample(1:length(trainingLabels), 1500)
  #mySampleTest <- sample(1:length(trainingLabels), 1500)
  n_train <- 3000
  n_test <- 1000
  myGeneSample <- sample(unique(trainingTxData$exonOverlapGeneId[!grepl('^R_*',trainingTxData$exonOverlapGeneId)]),n_train+n_test, replace=F)

  mySampleTrain <- trainingTxData$exonOverlapGeneId %in% myGeneSample[1:n_train]
  mySampleTrain <- trainingTxData$txSubsetNumber==0 & trainingTxData$exonOverlapGeneId %in% myGeneSample[1:n_train]
  mySampleTrain <- trainingTxData$txSubsetNumber>0 & trainingTxData$exonOverlapGeneId %in% myGeneSample[1:n_train]


  mySampleTest <- trainingTxData$exonOverlapGeneId %in% myGeneSample[-(1:n_train)]
  mySampleTest = T ## test all
  mySampleTest = (!grepl('^R_*',trainingTxData$exonOverlapGeneId))
  mySampleTest = trainingTxData$txSubsetNumber==0 & (!trainingTxData$exonOverlapGeneId %in% myGeneSample[1:n_train]) & (!grepl('^R_*',trainingTxData$exonOverlapGeneId))
  mySampleTest = trainingTxData$txSubsetNumber>0 & (!trainingTxData$exonOverlapGeneId %in% myGeneSample[1:n_train]) & (!grepl('^R_*',trainingTxData$exonOverlapGeneId))
  ### evaluation on set of genes sequins
  testGenes <- unique(trainingTxData$exonOverlapGeneId[grepl('^R_*',trainingTxData$exonOverlapGeneId)])
  mySampleTest <- trainingTxData$exonOverlapGeneId %in% testGenes



  modelList <- list()


  require(doMC)
  registerDoMC(cores=7)

  test_model1 <- fitBinomialModel(trainingLabels[mySampleTrain], as.matrix(model1_data[mySampleTrain,]), as.matrix(model1_data[mySampleTest,]), show.cv=TRUE, maxSize.cv=10000, parallel = T)
  coef(test_model1[[2]])
  show( fisher.test(table(test_model1[[1]]>0,trainingLabels[mySampleTest])))
  tmp=myPerformance(trainingLabels[mySampleTest]==1,test_model1[[1]])
  show(tmp$AUC)
  plot(tmp$FPR, tmp$TPR, ty='l')
  abline(0,1)
  lines(tmp$TPR, tmp$precision)


  test_model2 <- fitBinomialModel(trainingLabels[mySampleTrain], as.matrix(model2_data[mySampleTrain,]), as.matrix(model2_data[mySampleTest,]), show.cv=TRUE, maxSize.cv=10000, parallel = T)
  coef(test_model2[[2]])
  show( fisher.test(table(test_model2[[1]]>0,trainingLabels[mySampleTest])))
  tmp=myPerformance(trainingLabels[mySampleTest]==1,test_model2[[1]])
  show(tmp$AUC)
  lines(tmp$FPR, tmp$TPR, col=2)
  lines(tmp$TPR, tmp$precision, col=2)


  test_model3 <- fitBinomialModel(trainingLabels[mySampleTrain], as.matrix(model3_data[mySampleTrain,]), as.matrix(model3_data[mySampleTest,]), show.cv=TRUE, maxSize.cv=10000, parallel = T)
  coef(test_model3[[2]])
  show( fisher.test(table(test_model3[[1]]>0,trainingLabels[mySampleTest])))
  tmp=myPerformance(trainingLabels[mySampleTest]==1,test_model3[[1]])
  show(tmp$AUC)
  lines(tmp$FPR, tmp$TPR, col=3)
  lines(tmp$TPR, tmp$precision, col=3)


  f <- as.formula(y ~ .*.)
  y <- trainingLabels[mySampleTrain]
  x <- model.matrix(f, model3_data[mySampleTrain,])[, -1]
  x_test <- model.matrix(~.*., model3_data[mySampleTest,])[, -1]
  #  test=glmnet(x, y)
  # predictions= predict(test,newx=x_test)
  #cv.fit=cv.glmnet(x=x,y=y, parallel = T)
  cv.fit=glmnet::cv.glmnet(x=x,y=y,family='binomial', parallel = T)  ## slow but more correct
  predictions= predict(cv.fit,newx=x_test,s='lambda.min')


  # coef(cv.fit)
  show( fisher.test(table(predictions>0,trainingLabels[mySampleTest])))
  tmp=myPerformance(trainingLabels[mySampleTest]==1,predictions)
  show(tmp$AUC)
  lines(tmp$FPR, tmp$TPR, col=4)
  lines(tmp$TPR, tmp$precision, col=4)


  ##correct but slow!
  #test=glmnet(x, y, family='binomial')
  #predictions= predict(test,newx=x_test, type='response')
  show( fisher.test(table(predictions[,40]>0.5,trainingLabels[mySampleTest])))

  tmp=myPerformance(trainingLabels[mySampleTest]==1,predictions[,40])
  show(tmp$AUC)
  lines(tmp$FPR, tmp$TPR, col=4, lwd=2)
  lines(tmp$TPR, tmp$precision, col=4, lwd=2)




  cv.fit=glmnet::cv.glmnet(x=x,y=y, parallel = T)
  predictions= predict(cv.fit,newx=x_test,s='lambda.min')


  cv.fit=glmnet::cv.glmnet(x=x,y=y,family='binomial', parallel = T)
  predictions= predict(cv.fit,newx=x_test,s='lambda.min')



  coef(cv.fit)
  show( fisher.test(table(predictions>0,trainingLabels[mySampleTest])))
  tmp=myPerformance(trainingLabels[mySampleTest]==1,predictions)
  show(tmp$AUC)
  lines(tmp$FPR, tmp$TPR, col=5)
  lines(tmp$TPR, tmp$precision, col=5)


}
findPolySites <- function(startRanges, endRanges, exonsByTxList, barcodeData=NULL, windowSize=3, minRelativeReadCount = 0.05) {
  ## notes: works only on specifically sorted list of exon ranges, probably not on reference ranges. only tested on direct RNA seq, no barcode information or polyA information used. provides rough tables, need to be evluated and improved. relies on thresholds, can be removed?
  ## todo: create granges table with colData that contains all the information, upload to UCSC for visualisation along trancripts.
  ## todo: create this for TSS, create score to filter true transcripts from degradation products
  #windowSize=3
  #minRelativeReadCount = 0.05

  unlistedExons <- unlist(exonsByTxList)
  #unlistedExonsReference <- unlist(exonsByTx)

  #good example: ENSG00000147669 tx.52
  #chr8:100,149,797-100,154,844
  #exonsByTxExample = exonsByTxListNew['tx.52']

  #implement based on exon, then do based on read-transcript assigmnent table
  # unlistedExons <- unlist(exonsByTxExample)

  # create TES table for output, which variables should we have here?
  #tesTable <- matrix(NA,ncol=12,nrow=0)
  #colnames(tesTable) <- c('geneId','txId','tesId', 'readCount', 'barcode3PrimeCount','meanPolyALength','polyACountMin3','readFraction10Bp','readFraction1Bp', 'expectedEndCount', 'distToTssDistribution','ditanceToAnnotated')


  lastExonsMinus <- unlistedExons[strand(unlistedExons)=='-'&!duplicated(names(unlistedExons))]
  lastExonsMinusReduced <- reduce(lastExonsMinus,with.revmap=T)

  lastExonsPlus <- unlistedExons[strand(unlistedExons)=='+'&!duplicated(names(unlistedExons), fromLast=TRUE)]
  lastExonsPlusReduced <- reduce(lastExonsPlus,with.revmap=T)

  endOverlaps <- findOverlaps(lastExonsPlusReduced, endRanges)
  startOverlaps <-findOverlaps(lastExonsMinusReduced, startRanges)

  distToEnd <-end(lastExonsPlusReduced[queryHits(endOverlaps)])-end(endRanges[subjectHits(endOverlaps)])
  distToEndMinus <-start(startRanges[subjectHits(startOverlaps)])- start(lastExonsMinusReduced[queryHits(startOverlaps)])

  #nullDistribution <- prop.table(table(distToEnd))

  windowDist = ceiling((distToEnd+1)/windowSize)
  windowDistMinus = ceiling((distToEndMinus+1)/windowSize)


  #might not need this, use binomial distribution
  #nullDistribution <- (prop.table(table(distToEnd[distToEnd<100])))  # valid for plus and minus
  #nullDistributionWindow <- (prop.table(table(windowDist[distToEnd<100])))
  #exonDistTable$p = as.double(nullDistributionWindow[pmin(30,exonDistTable$windowDist)])
  #exonDistTable$expected <- exonDistTable$p*exonDistTable$readCount
  #exonDistTable$likelihood <- log(exonDistTable$p*exonDistTable$count)
  #exonDistTable$obsExp <- exonDistTable$count/exonDistTable$expected


  tesTablePlus=createTesTable(queryHits(endOverlaps), windowDist, windowSize=windowSize, minRelativeReadCount = minRelativeReadCount)
  tesTableMinus=createTesTable(queryHits(startOverlaps), windowDistMinus, windowSize=windowSize, minRelativeReadCount = minRelativeReadCount)

  return(list(tesTablePlus,tesTableMinus))
  ## from here: combine tables, merge more aggressively, use barcodes and pllyA, and use number of TSS per annotated transcripts as evaluation
  ## until here, looks good. reduce to Tes table with information

  #
  # tesTable <- exonDistTable3 %>% group_by(newExonSplit, lastExon) %>% mutate(min_pbinom_vs_pos1 = min(pbinom_count*tesCandidate),windowSizeTES = sum(tesCandidate)) %>% filter(tesCandidateFinal == TRUE) %>% ungroup() %>% select()
  #
  # exonDistTable3 %>% group_by(newExonSplit, lastExon) %>% mutate(min_pbinom_vs_pos1 = min(pbinom_count*tesCandidate))
  # #####
  #
  #
  # as.data.frame(exonDistTable %>% filter(lastExon==61))
  # as.data.frame(exonDistTable %>% filter(lastExon==483))
  #
  #
  # exonDistTable$p = as.double(nullDistribution[pmin(999,exonDistTable$distToEnd)+1])
  # exonDistTable$expected <- exonDistTable$p*exonDistTable$readCount
  #
  #
  # exonDistTable <- exonDistTable %>% group_by(lastExon, distToEnd) %>% mutate(count=n()) %>% distinct() %>% arrange(lastExon, distToEnd)
  # exonDistTable <- exonDistTable %>% group_by(lastExon) %>% mutate(readCount = sum(count))
  #
  # exonDistTable$p = as.double(nullDistribution[pmin(999,exonDistTable$distToEnd)+1])
  # exonDistTable$expected <- exonDistTable$p*exonDistTable$readCount
  # exonDistTable <- exonDistTable %>% mutate(relCount=count/readCount)
  #
  #
  # nullDistribution
  #
  #
  # unlistedExonsExample <- unlist(exonsByTxExample)
  #
  # lastExonsMinusExample <- unlistedExonsExample[strand(unlistedExonsExample)=='-'&!duplicated(names(unlistedExonsExample))]
  # lastExonsPlusExample <- unlistedExonsExample[strand(unlistedExonsExample)=='+'&!duplicated(names(unlistedExonsExample), fromLast=TRUE)]
  # endOverlapsExample <- findOverlaps(lastExonsPlusExample, endRanges)
  #
  # distToEndExample <-end(lastExonsPlusExample[queryHits(endOverlapsExample)])- end(endRanges[subjectHits(endOverlapsExample)])
  #
  #
  #
  #
  # firstExonsPlus <- unlistedExons[strand(unlistedExons)=='+'&!duplicated(names(unlistedExons))]
  # firstExonsMinus <- unlistedExons[strand(unlistedExons)=='-'&!duplicated(names(unlistedExons), fromLast=TRUE)]
  # startOverlaps <- findOverlaps(firstExonsPlus, startRanges)
  #
  # estBetaParams <- function(mu, var) {
  #   alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  #   beta <- alpha * (1 / mu - 1)
  #   return(params = list(alpha = alpha, beta = beta))
  # }
  # tmp=((exonDistTable %>% filter(windowDist==1) %>% mutate(relCount=count/readCount) %>% select(relCount))$relCount)
  # #only test 10 windos for TES
  # tmp=((exonDistTable %>% filter(windowDist<10) %>% group_by(lastExon) %>% mutate(readCount = sum(count))  %>% filter(windowDist==2, first(count)>10) %>% mutate(relCount=count/readCount) %>% select(relCount))$relCount)
  # exonDistTable2 = exonDistTable %>% filter(windowDist<10) %>% group_by(lastExon) %>% mutate(readCount = sum(count))  %>% filter(first(count)>10) %>% mutate(relCount=count/readCount)
  #
  # %>% select(relCount))$relCount)
  #
  #
  # a=estBetaParams(mean(tmp),var(tmp))$alpha
  # b=estBetaParams(mean(tmp),var(tmp))$beta
  # hist(tmp,col=2,breaks=100,freq=F)
  # plot(x,dbeta(x,a,b),col=1,lwd=2)
  #
  # a=NULL
  # b=NULL
  # for(i in 1:5){
  # a[i]
  # }
  #
  # parameterTable <- exonDistTable %>% mutate(windowDist=pmin(30,windowDist)) %>% group_by(windowDist) %>% summarise (mu = mean(relCount), var=var(relCount), a=((1 - mu) / var - 1 / mu) * mu ^ 2, b=a * (1 / mu - 1))
  # parameterTable3 <- exonDistTable %>%  group_by(lastExon) %>% filter(first(count)>10) %>% mutate(windowDist=pmin(30,windowDist)) %>% group_by(windowDist) %>% summarise (mu = mean(relCount), var=var(relCount), a=((1 - mu) / var - 1 / mu) * mu ^ 2, b=a * (1 / mu - 1))
  #
  # test=pbeta(exonDistTable$relCount, parameterTable$a[pmin(30,exonDistTable$windowDist)], parameterTable$b[pmin(30,exonDistTable$windowDist)])
  # exonDistTable$p_count <- dbeta(exonDistTable$relCount, parameterTable$a[pmin(30,exonDistTable$windowDist)], parameterTable$b[pmin(30,exonDistTable$windowDist)])
  #  exonDistTable$p_count <- (dbeta(exonDistTable$relCount, parameterTable$a[pmin(30,exonDistTable$windowDist)], parameterTable$b[pmin(30,exonDistTable$windowDist)],log=T))
  #

}


##function to create TES tables from overlap object
createTesTable <- function(lastExon, windowDist, windowSize=3, minRelativeReadCount = 0.05)
{
  #exonDistTable=tbl_df(data.frame(lastExon=queryHits(endOverlaps), windowDist))
  exonDistTable=tbl_df(data.frame(lastExon=lastExon, windowDist=windowDist))

  exonDistTable <- exonDistTable %>% group_by(lastExon, windowDist) %>% mutate(count=n()) %>% distinct() %>% arrange(lastExon, windowDist)
  exonDistTable <- exonDistTable %>% group_by(lastExon) %>% mutate(readCount = sum(count), maxCountByExon = max(count), firstCountByExon = dplyr::first(count))




  exonDistTable <- exonDistTable %>% mutate(relCount=count/readCount)

  #only calculate parameters using TES with minimum read count at first position and only 10 positions
  exonDistTable_parameter = exonDistTable %>% filter(windowDist<=10) %>% group_by(lastExon) %>% mutate(readCount = sum(count))  %>% filter(dplyr::first(count)>10) %>% mutate(relCount=count/readCount)

  #model for 10 windows? currently not used
  #parameterTable2 <- exonDistTable_parameter %>% group_by(windowDist) %>% summarise (mu = mean(relCount), var=var(relCount), a=((1 - mu) / var - 1 / mu) * mu ^ 2, b=a * (1 / mu - 1))

  #including longer distant counts is more accurate when they are not excluded
  parameterTable3 <- exonDistTable %>%  group_by(lastExon) %>% filter(dplyr::first(count)>10) %>% mutate(windowDist=pmin(30,windowDist)) %>% group_by(windowDist) %>% summarise (mu = mean(relCount), var=var(relCount), a=((1 - mu) / var - 1 / mu) * mu ^ 2, b=a * (1 / mu - 1))

  #exonDistTable$dbeta_count <- (dbeta(exonDistTable$relCount, parameterTable3$a[pmin(30,exonDistTable$windowDist)], parameterTable3$b[pmin(30,exonDistTable$windowDist)],log=T))
  #exonDistTable$dbinom_count <- (dbinom(exonDistTable$count,exonDistTable$readCount, parameterTable3$mu[pmin(30,exonDistTable$windowDist)],log=T))


  #probabiliyt that observation is greater than expected based on binomial distr
  exonDistTable$pbinom_count <- round(pbinom(exonDistTable$count,exonDistTable$readCount, parameterTable3$mu[pmin(30,exonDistTable$windowDist)],lower.tail = F,log=T),2)



  # find new candidate TES within last exon based on binomial distribution (p<1e-05)
  tesCount = length(unique((exonDistTable$lastExon)))
  exonDistTable = exonDistTable %>% ungroup() %>% mutate(newExonSplit = cumsum((pbinom_count < (log(0.001)) & count>1 & relCount > minRelativeReadCount) | windowDist==1 | (count==maxCountByExon & count>firstCountByExon ))) %>% group_by(newExonSplit) %>% mutate(newDist = windowDist - dplyr::first(windowDist) + 1, newReadCount = sum(count))
  exonDistTable$pbinom_count_new <- round(pbinom(exonDistTable$count,exonDistTable$newReadCount, parameterTable3$mu[pmin(30,exonDistTable$newDist)],lower.tail = F,log=T),2)


  while(tesCount < length(unique((exonDistTable$newExonSplit))))
  {
    tesCount <- length(unique((exonDistTable$newExonSplit)))
    show('find new candidate TES, number:')
    show(tesCount)
    exonDistTable = exonDistTable %>% ungroup() %>% mutate(newExonSplit = cumsum((pbinom_count_new < (log(0.001)) & count>1 & relCount > minRelativeReadCount) | newDist==1)) %>% group_by(newExonSplit) %>% mutate(newDist = windowDist - dplyr::first(windowDist) + 1, newReadCount = sum(count))
    exonDistTable$pbinom_count_new <- round(pbinom(exonDistTable$count,exonDistTable$newReadCount, parameterTable3$mu[pmin(30,exonDistTable$newDist)],lower.tail = F,log=T),2)
  }

  #remove candidates which have read count (sum of all reads ending only in this end) within this exon of less than 5%
  exonDistTable2 = exonDistTable %>% ungroup() %>% mutate(tesCandidate = (newDist==1 & (relCount > minRelativeReadCount | newReadCount/readCount > minRelativeReadCount)), newExonSplit = cumsum(tesCandidate)) %>% group_by(newExonSplit, lastExon) %>% mutate(newDist = windowDist - dplyr::first(windowDist) + 1, newReadCount = sum(count))

  #merge TES based on nearby distance, always choose the one from the longer transcript to be more inclusive (no length normalisation required)
  exonDistTable3 = exonDistTable2 %>% ungroup() %>%
    mutate(tesCandidateFinal =
             tesCandidate==TRUE &
             (lag(tesCandidate,default=FALSE)==FALSE | lastExon!=lag(lastExon,default=0) | (windowDist - lag(windowDist,default=0) >2)) &
             (lag(tesCandidate, 2, default=FALSE)==FALSE | lastExon!=lag(lastExon, 2,default=0) | (windowDist - lag(windowDist, 2, default=0) >2) | ((lag(relCount, 2, default=1) > 0.15 & lag(count, 2, default=0) >2) |relCount<lag(relCount, 2, default=1) )),
           newExonSplit = cumsum(tesCandidateFinal)) %>%
    group_by(newExonSplit, lastExon) %>%
    mutate(newDist = windowDist - dplyr::first(windowDist) + 1, newReadCount = sum(count))



  # second round of merge: if remaining TES are within windowDistance of 2, and ( relative count of the smaller one is below 15% or read Count of smaller one is 1), then also merge
  tesTable <- exonDistTable3 %>%
    group_by(newExonSplit, lastExon) %>%
    mutate( relCountByLastExon = relCount,
            readCountByLastExon = readCount,
            windowDistanceToMaxTES = windowDist,
            readCountAtFirstTESPosition = count,
            readCountByTES = newReadCount,
            mergedTESCandidates = sum(tesCandidate),
            readCountByMergedTESCandidates = sum(count*tesCandidate),
            relCountByMergedTES = sum(relCount*tesCandidate),
            windowsCoveredByMergedTES = max(windowDist*tesCandidate)-1/max(1/windowDist*tesCandidate)+1,
            min_pbinom_vs_pos1 = min(pbinom_count*tesCandidate),
            windowSizeTES = sum(tesCandidate)) %>%
    filter(tesCandidateFinal == TRUE) %>%
    ungroup() %>%
    dplyr::select(lastExon,windowDistanceToMaxTES,readCountAtFirstTESPosition, readCountByLastExon, relCountByLastExon, readCountByTES, min_pbinom_vs_pos1,mergedTESCandidates, readCountByMergedTESCandidates,windowsCoveredByMergedTES, relCountByMergedTES)
  return(tesTable)
}

##########HERE FUSION DETECTION FROM R??? CAN USE ALIGNMENT LATER?
findFusionCandidates<-function(bamFile, barcodeFile) ### just some test code... candidate reads can be identified from minimap alignment
{

  #fusion rules:
  #(1) 2 GenomicAlignments
  #(1.5) minumum distance: 50KB(?)
  #(2) alignment length adds up to read length (+-range OK); i.e. 3' S/H  of 5'partner + 5' S/H of 3' partner == qwidth / read length
  #(3) start barcode present
  #(4) end barcode present on same strand
  # (5) no identical alignments (pseudogenes, repeats etc)



  #barcodeFile='/mnt/ont/pychopper/poreDecode-bam/GIS_K562_cDNAStranded_Rep2-Run3_nanoporeDecode.rds'
  #bamFile='/mnt/ont/s3.ontdata.store.genome.sg/Nanopore/03_Mapping/Grch38//minimap2-2.1-cDNAStranded/GIS_K562_cDNAStranded_Rep2-Run3/GIS_K562_cDNAStranded_Rep2-Run3.bam'
  #barcodeFile='/mnt/ont/pychopper/poreDecode-bam/GIS_K562_cDNA_Rep3-Run4_nanoporeDecode.rds'
  #bamFile='/mnt/ont/s3.ontdata.store.genome.sg/Nanopore/03_Mapping/Grch38//minimap2-2.1-cDNA/GIS_K562_cDNA_Rep3-Run4/GIS_K562_cDNA_Rep3-Run4.bam'
  #barcodeFile='/mnt/ont/pychopper/poreDecode-bam/GIS_K562_cDNA_Rep3-Run4_nanoporeDecode.rds'
  #bamFile='/mnt/ont/s3.ontdata.store.genome.sg/Nanopore/03_Mapping/Grch38//minimap2-2.1-cDNA/GIS_K562_cDNA_Rep3-Run4/GIS_K562_cDNA_Rep3-Run4.bam'

  #barcodeFile='/mnt/ont/pychopper/poreDecode-bam/GIS_K562_directcDNA_1_nanoporeDecode.rds'
  #bamFile='/mnt/ont/s3.ontdata.store.genome.sg/Nanopore/03_Mapping/Grch38//minimap2-2.1-cDNA/GIS_K562_directcDNA_1/GIS_K562_directcDNA_1.bam'

  reads=readGAlignments(bamFile, param=ScanBamParam(flag=scanBamFlag(isSecondaryAlignment=FALSE)), use.names=TRUE)
  readsDuplicated <- reads[duplicated(names(reads))|duplicated(names(reads),fromLast=TRUE)]

  mcols(readsDuplicated)$readName <- names(readsDuplicated)
  names(readsDuplicated)<-NULL
  cigOpLengths <- explodeCigarOpLengths(cigar(readsDuplicated), ops=CIGAR_OPS)
  cigOps <- explodeCigarOps(cigar(readsDuplicated), ops=CIGAR_OPS)
  mcols(readsDuplicated)$cigOp1 <- unlist(lapply(cigOps,'[[',1))
  mcols(readsDuplicated)$cigOpLength1 <- unlist(lapply(cigOpLengths,'[[',1))
  mcols(readsDuplicated)$cigOp2 <- unlist(lapply(cigOps,function(x) x[length(x)]))
  mcols(readsDuplicated)$cigOpLength2 <- unlist(lapply(cigOpLengths,function(x) x[length(x)]))


  readTable <- tbl_df(as.data.frame(readsDuplicated)) %>% group_by(readName) %>% filter(n()==2) %>% arrange(readName, cigOpLength1)
  readTable <- readTable %>% mutate(combinedLength = dplyr::first(cigOpLength2)+last(cigOpLength1), readLength = max(qwidth), dist=max(start)-min(end), equalChr = dplyr::first(seqnames)==last(seqnames) )

  if(!is.null(barcodeFile)){
    barcodeData <- readRDS(file=barcodeFile)
    readTable=inner_join(readTable,tbl_df(barcodeData),by=c('readName'='names'))
    readTableFiltered = readTable %>% filter(equalChr==FALSE, readLength==combinedLength, barcode5Prime.plus.score>10, barcode3Prime.plus.score>10 )

  } else {
    readTableFiltered = readTable %>% filter(equalChr==FALSE, abs(readLength-combinedLength) < 10)
  }
  readTableFiltered <- mutate(readTableFiltered,seqCombo=paste(dplyr::first(seqnames), last(seqnames), sep=':'))


  #
  # readTable2 = readTable %>% group_by(readName) %>% filter(n()==2)
  #
  # readTable3 <- readTable2 %>% arrange(readName, cigOpLength1) %>% mutate(combinedLength = dplyr::first(cigOpLength2)+last(cigOpLength1), readLength = max(qwidth), dist=max(start)-min(end), equalChr = first(seqnames)==last(seqnames) )
  #
  # readTable4=inner_join(readTable3,tbl_df(barcodeData),by=c('readName'='names'))
  #
  # readTable5 = readTable4 %>% filter(equalChr==FALSE, readLength==combinedLength, barcode5Prime.plus.score>10, barcode3Prime.plus.score>10 )
  #
  #
  # readSubset <- readGAlignments(bamFile, param=ScanBamParam(which=unlist(readGrglist[names(readGrglist) %in% readTable5$readName]),what = scanBamWhat()),use.names = T)
  # readSubsetTable5 <- readSubset[which(names(readSubset) %in% readTable5$readName)]
  # readSubsetTable5 <- readSubsetTable5[!duplicated(mcols(readSubsetTable5)$seq)]

}





#
# makeGRangesListFromFeatureFragments <- function(seqnames=Rle(factor()),
#                                                 fragmentStarts=list(),
#                                                 fragmentEnds=list(),
#                                                 fragmentWidths=list(),
#                                                 strand=character(0),
#                                                 sep=",")
# {
#     fragmentStarts <- normargListOfIntegers(fragmentStarts, sep,
#                                             "fragmentStarts")
#     nfrag_per_feature <- elementNROWS(fragmentStarts)
#     start <- unlist(fragmentStarts, recursive=FALSE, use.names=FALSE)
#
#     fragmentEnds <- normargListOfIntegers(fragmentEnds, sep,
#                                           "fragmentEnds")
#     nend_per_elt <- elementNROWS(fragmentEnds)
#     if (length(nend_per_elt) != 0L) {
#         if (length(nfrag_per_feature) == 0L)
#             nfrag_per_feature <- nend_per_elt
#         else if (!identical(nend_per_elt, nfrag_per_feature))
#             stop("'fragmentStarts' and 'fragmentEnds' have ",
#                  "incompatible \"shapes\"")
#     }
#     end <- unlist(fragmentEnds, recursive=FALSE, use.names=FALSE)
#
#     fragmentWidths <- normargListOfIntegers(fragmentWidths, sep,
#                                             "fragmentWidths")
#     nwidth_per_elt <- elementNROWS(fragmentWidths)
#     if (length(nwidth_per_elt) != 0L) {
#         if (length(nfrag_per_feature) == 0L)
#             nfrag_per_feature <- nwidth_per_elt
#         else if (!identical(nwidth_per_elt, nfrag_per_feature))
#             stop("\"shape\" of 'fragmentWidths' is incompatible ",
#                  "with \"shape\" of 'fragmentStarts' or 'fragmentEnds'")
#     }
#     width <- unlist(fragmentWidths, recursive=FALSE, use.names=FALSE)
#
#     ranges <- IRanges(start=start, end=end, width=width)
#     nfrag <- sum(nfrag_per_feature)
#     if (nfrag != length(ranges))
#         stop("GenomicRanges internal error in makeGRangesListFromFields(): ",
#              "nfrag != length(ranges). This should never happen. ",
#              "Please report.")
#     if (nfrag == 0L) {
#         ## Cannot blindly subset by FALSE because it doesn't work on a
#         ## zero-length Rle.
#         if (length(seqnames) != 0L)
#             seqnames <- seqnames[FALSE]
#         if (length(strand) != 0L)
#             strand <- strand[FALSE]
#     } else {
#         if (length(seqnames) != length(nfrag_per_feature) ||
#             length(strand) != length(nfrag_per_feature))
#             stop("length of 'seqnames' and/or 'strand' is incompatible ",
#                  "with fragmentStarts/Ends/Widths")
#         seqnames <- rep.int(seqnames, nfrag_per_feature)
#         strand <- rep.int(strand, nfrag_per_feature)
#     }
#     unlistData <- GRanges(seqnames=seqnames, ranges=ranges, strand=strand)
#     partitioning <- PartitioningByEnd(cumsum(nfrag_per_feature), names=NULL)
#     relist(unlistData, partitioning)
# }
#
showTopLines <- function(x, n=5) { show(as.data.frame(head(x,n)))}
