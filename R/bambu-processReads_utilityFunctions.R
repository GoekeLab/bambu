
#' Isoform reconstruction using genomic alignments
#' @param readGrgList readGrgList
#' @param runName runName
#' @param stranded stranded
#' @param quickMode quickMode
#' @inheritParams bambu
#' @noRd
isore.constructReadClasses <- function(readGrgList,
                                       runName = "sample1", annotationGrangesList,
                                       genomeSequence = NULL, stranded = FALSE, ncore = 1, verbose = FALSE) {
  unlisted_junctions <- unlistIntrons(readGrgList,
                                      use.ids = TRUE, use.names = FALSE)
  start.ptm <- proc.time()
  uniqueJunctions <- createJunctionTable(unlisted_junctions,
                                         ncore = ncore, genomeSequence = genomeSequence)
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
