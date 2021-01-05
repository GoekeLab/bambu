
#' Isoform reconstruction using genomic alignments
#' @param readGrgList readGrgList
#' @param runName runName
#' @param stranded stranded
#' @param quickMode quickMode
#' @inheritParams bambu
#' @noRd
isore.constructReadClasses <- function(readGrgList,
                                       runName = "sample1", annotationGrangesList,
                                       genomeSequence = NULL, stranded = FALSE, verbose = FALSE) {
  unlisted_junctions <- unlistIntrons(readGrgList,
                                      use.ids = TRUE, use.names = FALSE)
  start.ptm <- proc.time()
  uniqueJunctions <- createJunctionTable(unlisted_junctions,
                                         genomeSequence = genomeSequence)
  # all seqlevels should be consistent, and drop those not in uniqueJunctions
  if (!all(GenomeInfoDb::seqlevels(unlisted_junctions) %in% 
           GenomeInfoDb::seqlevels(uniqueJunctions))) {
    unlisted_junctions <- GenomeInfoDb::keepSeqlevels(unlisted_junctions,
                                                      value = GenomeInfoDb::seqlevels(unlisted_junctions)[
                                                        GenomeInfoDb::seqlevels(unlisted_junctions) %in%
                                                          GenomeInfoDb::seqlevels(uniqueJunctions)], pruning.mode = "coarse")
    readGrgList <- GenomeInfoDb::keepSeqlevels(readGrgList,
                                               value = GenomeInfoDb::seqlevels(readGrgList)[ 
                                                 GenomeInfoDb::seqlevels(readGrgList) %in%
                                                   GenomeInfoDb::seqlevels(uniqueJunctions)], pruning.mode = "coarse")
  } # the seqleels will be made comparable for all ranges,
  # warning is shown if annotation is missing some
  if (!all(GenomeInfoDb::seqlevels(readGrgList) %in% 
           GenomeInfoDb::seqlevels(annotationGrangesList))) 
    message("not all chromosomes present in reference annotations,
            annotations might be incomplete. Please compare objects
            on the same reference")
  end.ptm <- proc.time()
  if (verbose) message("Finished creating junction list with splice motif
        in ", round((end.ptm - start.ptm)[3] / 60, 1), " mins.")
  exonsByReadClass <- generateExonsByReadClass(readGrgList,
                                               annotationGrangesList, unlisted_junctions, uniqueJunctions,
                                               stranded, verbose)
  counts <- matrix(mcols(exonsByReadClass)$readCount,
                   dimnames = list(names(exonsByReadClass), runName))
  colDataDf <- DataFrame(name = runName, row.names = runName)
  mcols(exonsByReadClass) <- mcols(exonsByReadClass)[, c("chr.rc", 
                                                         "strand.rc", "intronStarts", "intronEnds", "confidenceType")]
  se <- SummarizedExperiment(assays = SimpleList(counts = counts),
                             rowRanges = exonsByReadClass, colData = colDataDf)
  return(se)
}


#' Get unlisted intron ranges from exon ranges list
#' @noRd
unlistIntrons <- function(x, use.ids = TRUE, use.names = TRUE) {
  # License note: This function is adopted from the GenomicAlignments 
  # package (Author: Hervé Pagès, Valerie Obenchain, Martin Morgan)
  # License Artistic-2.0
  # https://doi.org/doi:10.18129/B9.bioc.GenomicAlignments
  
  flat <- unlist(x, use.names = FALSE)
  gaps <- gaps(ranges(x))
  
  firstseg <- start(PartitioningByWidth(x))
  seqnms <- rep(seqnames(flat)[firstseg], elementNROWS(gaps))
  strand <- rep(strand(flat)[firstseg], elementNROWS(gaps))
  
  gr <- GenomicRanges::GRanges(seqnms, unlist(gaps,
                                              use.names = use.names), strand)
  if (use.ids & !is.null(mcols(x, use.names = FALSE)$id)) 
    mcols(gr)$id <- rep(mcols(x)$id, elementNROWS(gaps))
  return(gr)
}


#' generate exonByReadClass
#' @noRd
generateExonsByReadClass <- function(readGrgList, annotationGrangesList, 
                                     unlisted_junctions, uniqueJunctions, stranded, verbose){
  GenomeInfoDb::seqlevels(readGrgList) <-
    unique(c(GenomeInfoDb::seqlevels(readGrgList),
             GenomeInfoDb::seqlevels(annotationGrangesList)))
  GenomeInfoDb::seqlevels(annotationGrangesList) <- 
    GenomeInfoDb::seqlevels(readGrgList)
  readClassListSpliced <- createModelforJunctionReads(
    readGrgList, annotationGrangesList, unlisted_junctions,
    uniqueJunctions, stranded, verbose)
  # seqlevels are made equal (added for chromosomes missing in any of them)
  start.ptm <- proc.time()
  singleExonReads <- unlist(readGrgList[elementNROWS(readGrgList) == 1],
                            use.names = FALSE)
  mcols(singleExonReads)$id <- mcols(readGrgList[
    elementNROWS(readGrgList) == 1])$id
  referenceExons <- unique(c(GenomicRanges::granges(unlist(
    readClassListSpliced[mcols(readClassListSpliced)$confidenceType ==
                           "highConfidenceJunctionReads" &
                           mcols(readClassListSpliced)$strand.rc != "*"], use.names = FALSE)), 
    GenomicRanges::granges(unlist(annotationGrangesList,
                                  use.names = FALSE))))
  readClassListUnsplicedWithAnnotation <- constructUnsplicedReadClasses(
    granges = singleExonReads, grangesReference = referenceExons,
    confidenceType = "unsplicedWithin", stranded = stranded)
  singleExonReads <- singleExonReads[!mcols(singleExonReads)$id %in%
                                       readClassListUnsplicedWithAnnotation$readIds]
  referenceExons <- reduce(singleExonReads, ignore.strand = !stranded)
  readClassListUnsplicedReduced <- constructUnsplicedReadClasses(
    granges = singleExonReads, grangesReference = referenceExons,
    confidenceType = "unsplicedNew", stranded = stranded)
  end.ptm <- proc.time()
  if (verbose) message("Finished create single exon transcript models
        (read classes) in ", round((end.ptm - start.ptm)[3] / 60, 1), " mins.")
  exonsByReadClass <- c(readClassListSpliced,
                        readClassListUnsplicedWithAnnotation$exonsByReadClass,
                        readClassListUnsplicedReduced$exonsByReadClass)
  return(exonsByReadClass)
}

#' create transcript model for splice junction reads
#' @param readGrgList reads GRangesList
#' @param annotationGrangesList annotation GRangesList
#' @param unlisted_junctions unlisted_junctions
#' @param uniqueJunctions uniqueJunctions
#' @param stranded stranded
#' @param verbose verbose
#' @noRd
createModelforJunctionReads <- function(readGrgList, annotationGrangesList,
                                        unlisted_junctions, uniqueJunctions, stranded, verbose) {
  GenomeInfoDb::seqlevels(unlisted_junctions) <-
    GenomeInfoDb::seqlevels(readGrgList)
  GenomeInfoDb::seqlevels(uniqueJunctions) <- 
    GenomeInfoDb::seqlevels(readGrgList)
  uniqueAnnotatedIntrons <- unique(unlistIntrons(annotationGrangesList,
                                                 use.names = FALSE, use.ids = FALSE))
  junctionTables <- junctionStrandCorrection(uniqueJunctions,
                                             unlisted_junctions, uniqueAnnotatedIntrons,
                                             stranded = stranded, verbose = verbose)
  uniqueJunctions <- junctionTables[[1]][, c("score", "spliceMotif",
                                             "spliceStrand", "junctionStartName", "junctionEndName",
                                             "startScore", "endScore", "id")]
  unlisted_junctions <- junctionTables[[2]]
  uniqueJunctions$annotatedJunction <- (!is.na(GenomicRanges::match(
    uniqueJunctions,uniqueAnnotatedIntrons)))
  uniqueJunctions$annotatedStart <- uniqueJunctions$junctionStartName %in%
    uniqueJunctions$junctionStartName[uniqueJunctions$annotatedJunction]
  uniqueJunctions$annotatedEnd <- uniqueJunctions$junctionEndName %in%
    uniqueJunctions$junctionEndName[uniqueJunctions$annotatedJunction]
  uniqueJunctions <- correctJunctionFromPrediction(uniqueJunctions, verbose)
  start.ptm <- proc.time()
  readClassListSpliced <- constructSplicedReadClassTables(
    uniqueJunctions = uniqueJunctions,
    unlisted_junctions = unlisted_junctions,
    readGrgList = readGrgList,
    stranded = stranded)
  end.ptm <- proc.time()
  if (verbose)
    message("Finished create transcript models (read classes) for reads with
    spliced junctions in ", round((end.ptm - start.ptm)[3] / 60, 1)," mins.")
  return(readClassListSpliced)
}


#' correct junction from prediction
#' @param uniqueJunctions uniqueJunctions
#' @param verbose verbose
#' @noRd
correctJunctionFromPrediction <- function(uniqueJunctions, verbose) {
  start.ptm <- proc.time()
  if (sum(uniqueJunctions$annotatedJunction) > 5000 &
      sum(!uniqueJunctions$annotatedJunction) > 4000) {
    uniqJunctionsNmodels <-
      findUniqueJunctions(uniqueJunctions, NULL, verbose)
    uniqueJunctions <- uniqJunctionsNmodels$uniqueJunctions
    junctionModel <- uniqJunctionsNmodels$predictSpliceSites[[2]]
  } else {
    junctionModel <- standardJunctionModels_temp
    uniqJunctionsNmodels <-
      findUniqueJunctions(uniqueJunctions, junctionModel, verbose)
    uniqueJunctions <- uniqJunctionsNmodels$uniqueJunctions
    message("Junction correction with not enough data,
            precalculated model is used")
  }
  end.ptm <- proc.time()
  if (verbose) 
    message("Model to predict true splice sites built in ",
            round((end.ptm - start.ptm)[3] / 60, 1), " mins.")
  start.ptm <- proc.time()
  uniqueJunctions <- findHighConfidenceJunctions( junctions = uniqueJunctions,
                                                  junctionModel = junctionModel, verbose = verbose)
  uniqueJunctions$mergedHighConfJunctionIdAll_noNA <- 
    uniqueJunctions$mergedHighConfJunctionId
  uniqueJunctions$mergedHighConfJunctionIdAll_noNA[
    is.na(uniqueJunctions$mergedHighConfJunctionId)] <- 
    names(uniqueJunctions[is.na(uniqueJunctions$mergedHighConfJunctionId)])
  uniqueJunctions$strand.mergedHighConfJunction <- 
    as.character(strand(
      uniqueJunctions[uniqueJunctions$mergedHighConfJunctionIdAll_noNA]))
  end.ptm <- proc.time()
  if (verbose) 
    message("Finished correcting junction based on set of high confidence
            junctions in ", round((end.ptm - start.ptm)[3] / 60, 1), " mins.")
  return(uniqueJunctions)
}


#' find unique junctions
#' @noRd
findUniqueJunctions <- function(uniqueJunctions, junctionModel, verbose){
  predictSpliceSites <- predictSpliceJunctions(
    annotatedJunctions = uniqueJunctions,
    junctionModel = junctionModel,
    verbose = verbose)
  uniqueJunctions <- predictSpliceSites[[1]][, c(
    "score", "spliceMotif",
    "spliceStrand", "junctionStartName", "junctionEndName",
    "startScore", "endScore", "annotatedJunction",
    "annotatedStart", "annotatedEnd")]
  return(list("uniqueJunctions" = uniqueJunctions, 
              "predictSpliceSites" = predictSpliceSites))
}
