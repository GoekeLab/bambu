#' Isoform reconstruction using genomic alignments
#' @param unlisted_junctions unlisted_junctions
#' @param annotations annotations
#' @param genomeSequence genomeSequence
#' @param stranded stranded
#' @param verbose verbose
#' @inheritParams bambu
#' @noRd
isore.constructJunctionTables <- function(unlisted_junctions, annotations,
                                          genomeSequence, stranded = FALSE, 
                                          verbose = FALSE) {
  start.ptm <- proc.time()
  #summarise junction counts and strand for all reads
  uniqueJunctions <- createJunctionTable(unlisted_junctions,
                                         genomeSequence = genomeSequence)
  end.ptm <- proc.time()
  if (verbose) message("Finished creating junction list with splice motif
        in ", round((end.ptm - start.ptm)[3] / 60, 1), " mins.")
  
  uniqueAnnotatedIntrons <- unique(unlistIntrons(annotations, 
                                                 use.ids = FALSE))
  # correct strand of junctions based on (inferred) strand of reads
  strand(uniqueJunctions) <- junctionStrandCorrection(uniqueJunctions,
                                                      unlisted_junctions, 
                                                      uniqueAnnotatedIntrons,
                                                      stranded = stranded, 
                                                      verbose = verbose)
  # add annotation labels to junctions
  mcols(uniqueJunctions) <- tibble(as.data.frame(uniqueJunctions)) %>% 
    mutate(annotatedJunction = (!is.na(GenomicRanges::match(uniqueJunctions,
                                                            uniqueAnnotatedIntrons)))) %>% 
    group_by(seqnames) %>% 
    mutate(annotatedStart = start %in% start[annotatedJunction],
           annotatedEnd = end %in% end[annotatedJunction]) %>% ungroup() %>%
    select(score, spliceMotif, spliceStrand, junctionStartName, junctionEndName,
           startScore, endScore, id, annotatedJunction, annotatedStart, annotatedEnd)
  # correct junction coordinates using logistic regression classifier
  uniqueJunctions <- correctJunctionFromPrediction(uniqueJunctions, verbose)
  return(uniqueJunctions)
}



#' Isoform reconstruction using genomic alignments
#' @param readGrgList readGrgList
#' @param runName runName
#' @param stranded stranded
#' @inheritParams bambu
#' @noRd
  isore.constructReadClasses <- function(readGrgList, runName = "sample1", 
                                         annotationGrangesList, 
                                         genomeSequence = NULL, stranded = FALSE, 
                                         verbose = FALSE) {
  unlisted_junctions <- unlistIntrons(readGrgList, use.ids = TRUE)
  start.ptm <- proc.time()
  uniqueJunctions <- createJunctionTable(unlisted_junctions,
                                         genomeSequence = genomeSequence)
  end.ptm <- proc.time()
  if (verbose) message("Finished creating junction list with splice motif
        in ", round((end.ptm - start.ptm)[3] / 60, 1), " mins.")
  
 #  # all seqlevels should be consistent, and drop those not in uniqueJunctions
 #  if (!all(GenomeInfoDb::seqlevels(unlisted_junctions) %in% 
 #           GenomeInfoDb::seqlevels(uniqueJunctions))) {
 #    unlisted_junctions <- keepSeqLevelsByReference(unlisted_junctions, uniqueJunctions)
 #    readGrgList <- keepSeqLevelsByReference(readGrgList, uniqueJunctions)
 #  } # the seqlevels will be made comparable for all ranges,
 #  if (!all(GenomeInfoDb::seqlevels(readGrgList) %in% 
 #           GenomeInfoDb::seqlevels(annotationGrangesList))) 
 #    message("not all chromosomes present in reference annotations,
 #            annotations might be incomplete. Please compare objects
 #            on the same reference")
 # combinedSeqLevels <- unique(c(GenomeInfoDb::seqlevels(readGrgList),
 #                               GenomeInfoDb::seqlevels(annotationGrangesList)))
 #  GenomeInfoDb::seqlevels(readGrgList) <- combinedSeqLevels
 #  GenomeInfoDb::seqlevels(annotationGrangesList) <- combinedSeqLevels
 #  GenomeInfoDb::seqlevels(unlisted_junctions) <- combinedSeqLevels
 #  GenomeInfoDb::seqlevels(uniqueJunctions) <- combinedSeqLevels
  uniqueAnnotatedIntrons <- unique(unlistIntrons(annotationGrangesList, 
                                                 use.ids = FALSE))
  strand(uniqueJunctions) <- junctionStrandCorrection(uniqueJunctions,
                                                      unlisted_junctions, 
                                                      uniqueAnnotatedIntrons,
                                                      stranded = stranded, 
                                                      verbose = verbose)
  
  mcols(uniqueJunctions) <- tibble(as.data.frame(uniqueJunctions)) %>% 
    mutate(annotatedJunction = (!is.na(GenomicRanges::match(uniqueJunctions, uniqueAnnotatedIntrons)))) %>% group_by(seqnames) %>% 
    mutate(annotatedStart = start %in% start[annotatedJunction],
           annotatedEnd = end %in% end[annotatedJunction]) %>% ungroup() %>%
    select(score, spliceMotif, spliceStrand, junctionStartName, junctionEndName,
           startScore, endScore, id, annotatedJunction, annotatedStart, annotatedEnd)
  
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

  exonsByReadClass <- generateExonsByReadClass(readGrgList, 
                                               annotationGrangesList, 
                                               readClassListSpliced, 
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
unlistIntrons <- function(x, use.ids = TRUE, use.names = FALSE) {
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


#' Create Junction tables from unlisted junction granges
#' @importFrom BiocParallel bppram bpvec
#' @noRd
createJunctionTable <- function(unlisted_junctions,
                                genomeSequence = NULL) {
  # License note: This function is adopted from the GenomicAlignments package 
  #genomeSequence <- checkInputSequence(genomeSequence)
  if (!all(GenomeInfoDb::seqlevels(unlisted_junctions) %in%
           GenomeInfoDb::seqlevels(genomeSequence))) {
    message("not all chromosomes present in reference genome sequence,
            ranges are dropped (createJunctionTable())")
    unlisted_junctions <- GenomeInfoDb::keepSeqlevels(unlisted_junctions,
                                                      value = GenomeInfoDb::seqlevels(unlisted_junctions)[
                                                        GenomeInfoDb::seqlevels(unlisted_junctions) %in% 
                                                          GenomeInfoDb::seqlevels(genomeSequence)],
                                                      pruning.mode = "coarse")
  }
  
  uniqueJunctions <- sort(unique(BiocGenerics::unstrand(unlisted_junctions)))
  names(uniqueJunctions) <- paste("junc", seq_along(uniqueJunctions),
                                  sep = ".")
  plus_score <- countMatches(uniqueJunctions,
                             unlisted_junctions[strand(unlisted_junctions)=='+'], 
                             ignore.strand=T)
  minus_score <- countMatches(uniqueJunctions,
                              unlisted_junctions[strand(unlisted_junctions)=='-'], 
                              ignore.strand=T)

  junctionSeqStart <- BSgenome::getSeq(genomeSequence,
                                       IRanges::shift(flank(uniqueJunctions,width = 2), 2))#shift from IRanges
  junctionSeqEnd <- BSgenome::getSeq(genomeSequence,
                                     IRanges::shift(flank(uniqueJunctions,width = 2, start = FALSE), -2))
 
  mcols(uniqueJunctions) <- DataFrame(tibble(
    chr=as.factor(seqnames(uniqueJunctions)), 
    start=start(uniqueJunctions),
    end=end(uniqueJunctions),
    score=plus_score+minus_score,
    plus_score = plus_score,
    minus_score = minus_score,
    spliceMotif=paste(junctionSeqStart, junctionSeqEnd, sep = "-"),
    spliceStrand = spliceStrand(spliceMotif),
    junctionStartName = paste(chr, start, sep = ":"),
    junctionEndName = paste(chr, end, sep = ":"),
    id = seq_along(uniqueJunctions)) %>%
      group_by(chr, start) %>% 
      mutate(startScore=sum(score)) %>% 
      group_by(chr, end) %>%  
      mutate(endScore=sum(score)) %>%
      ungroup() %>%
      select(score, plus_score, minus_score, spliceMotif, spliceStrand,
             junctionStartName, junctionEndName, startScore, endScore, id))
  strand(uniqueJunctions) <- uniqueJunctions$spliceStrand
  return(uniqueJunctions)
}

# can be removed since not used anymore
#' #' Keep only the seqlevels that also occur in the reference granges
#' #' @noRd
#' keepSeqLevelsByReference <- function(gr, ref) {
#'   gr <- GenomeInfoDb::keepSeqlevels(gr,
#'                                     value = GenomeInfoDb::seqlevels(gr)[
#'                                       GenomeInfoDb::seqlevels(gr) %in% 
#'                                         GenomeInfoDb::seqlevels(ref)],
#'                                     pruning.mode = "coarse")
#'   return(gr)
#' }

#' Function to create a object that can be queried by getSeq
#' Either from fa file, or BSGenome object
#' @importFrom BiocParallel bppram bpvec
#' @noRd
checkInputSequence <- function(genomeSequence) {
  if (is.null(genomeSequence)) stop("Reference genome sequence is missing,
        please provide fasta file or BSgenome name, see available.genomes()")
  if (methods::is(genomeSequence, "character")) {
    if (grepl(".fa", genomeSequence)) {
      if (.Platform$OS.type == "windows") {
        genomeSequence <- Biostrings::readDNAStringSet(genomeSequence)
        newlevels <- unlist(lapply(strsplit(names(genomeSequence)," "),
                                   "[[", 1))
        names(genomeSequence) <- newlevels
      } else {
        genomeSequence <- Rsamtools::FaFile(genomeSequence)
      }
    } else {
      genomeSequence <- BSgenome::getBSgenome(genomeSequence)
    }
  }
  return(genomeSequence)
}

#' @param motif motif
#' @noRd
spliceStrand <- function(motif) {
  NATURAL_INTRON_MOTIFS_RC <- as.character(Biostrings::reverseComplement(
    Biostrings::DNAStringSet(GenomicAlignments::NATURAL_INTRON_MOTIFS)))
  
  motifStrand <- ifelse(motif %in% GenomicAlignments::NATURAL_INTRON_MOTIFS,
                        "+", "*")
  motifStrand[motif %in% NATURAL_INTRON_MOTIFS_RC] <- "-"
  return(motifStrand)
}

#' generate exonByReadClass
#' @noRd
generateExonsByReadClass <- function(readGrgList, annotationGrangesList, readClassListSpliced, stranded, verbose){

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


#' JUNCTIONSTRANDCORRECTION
#' @noRd
junctionStrandCorrection <- function(uniqueJunctions, unlisted_junctions,
                                     uniqueAnnotatedIntrons, stranded, verbose = FALSE) {
  # note: strand sometimes incorrectly infered based on motifs, might 
  # introduce systematic errors due to alignment (biased to splice motifs)
  
  uniqueJunctionsUpdate <- uniqueJunctions
  # make a copy to revert to if strand correction does not improve results
  annotatedIntronNumber <- evalAnnotationOverlap(uniqueJunctions,
                                                 uniqueAnnotatedIntrons,ignore.strand = FALSE)["TRUE"]
  if (verbose) {
    message("before strand correction, annotated introns:")
    message(annotatedIntronNumber)
    message(annotatedIntronNumber / length(uniqueJunctions))
  }
  # infer strand for each read based on strand of junctions
  strandStep <- TRUE
  while (strandStep) { # iterate twice to improve strand prediction w.t.
    # mean junction counts, annotate junction strand with read strand
    if (stranded==FALSE) { # update junction strand score
      
      
      strandScoreByRead <- 
        updateStrandScoreByRead(unlisted_junctions,
                                uniqueJunctionsUpdate)
    } else{ # just use strand from reads for stranded data
      strandScoreByRead <- uniqueJunctionsUpdate$minus_score -
        uniqueJunctionsUpdate$plus_score
    }
    # overwrite info from motif which increases overlap with known junc
    strand(uniqueJunctionsUpdate[strandScoreByRead < 0]) <- "+"
    strand(uniqueJunctionsUpdate[strandScoreByRead > 0]) <- "-"
    updatedList <- updateJunctionwimprove(annotatedIntronNumber,
                                          uniqueJunctions, uniqueJunctionsUpdate, uniqueAnnotatedIntrons,
                                          strandStep, verbose)
    annotatedIntronNumber <- updatedList$annotatedIntronNumber
    uniqueJunctions <- updatedList$uniqueJunctions
    strandStep <- updatedList$strandStep
  }
  return(strand(uniqueJunctions))
}


#' Evaluate annoation overlap
#' @noRd
evalAnnotationOverlap <- function(intronRanges, uniqueAnnotatedIntrons,
                                  ignore.strand = FALSE) {
  return(table(!is.na(GenomicRanges::match(intronRanges,
                                           uniqueAnnotatedIntrons, ignore.strand = ignore.strand ))))
}


#' This function assigns a strand to each read based on the majority of junctions
#' The strand of the junctions is infered by the sequence in createJunctionTables
#' @noRd
updateStrandScoreByRead <- function(unlisted_junctions,
                                    uniqueJunctions){

  allJunctionToUniqueJunctionMatch <- match(unlisted_junctions, uniqueJunctions, ignore.strand=T)
  
  unlisted_junction_granges_strandList <-
    splitAsList(strand(uniqueJunctions)[
      allJunctionToUniqueJunctionMatch],
      mcols(unlisted_junctions)$id)
  
  strandJunctionSum <-
    as.integer(sum(unlisted_junction_granges_strandList == "-") -
    sum(unlisted_junction_granges_strandList == "+"))
  
   readStrand <- factor(rep("*", length(unlisted_junction_granges_strandList)), levels=c('+','-','*'))
   readStrand[strandJunctionSum < 0] <- "+"
   readStrand[strandJunctionSum > 0] <- "-"
   
  
   strand_unlisted_junctions <-
     readStrand[match(mcols(unlisted_junctions)$id,
                      as.integer(names(unlisted_junction_granges_strandList)))]
  plus_score <- countMatches(uniqueJunctions,
                             unlisted_junctions[strand_unlisted_junctions =='+'], 
                             ignore.strand=T)
  minus_score <- countMatches(uniqueJunctions,
                              unlisted_junctions[strand_unlisted_junctions =='-'], 
                              ignore.strand=T)
    return(minus_score-plus_score)
}



#' update junctions object if strand prediction improves overlap 
#' with annotations
#' @param annotatedIntronNumber annotatedIntronNumber
#' @param uniqueJunctions uniqueJunctions
#' @param uniqueJunctionsUpdate uniqueJunctionsUpdate
#' @param uniqueAnnotatedIntrons uniqueAnnotatedIntrons
#' @param strandStep strandStep
#' @param verbose verbose
#' @noRd
updateJunctionwimprove <- function(annotatedIntronNumber, uniqueJunctions,
                                   uniqueJunctionsUpdate, uniqueAnnotatedIntrons, strandStep, verbose) {
  annotatedIntronNumberNew <- evalAnnotationOverlap(uniqueJunctionsUpdate,
                                                    uniqueAnnotatedIntrons, ignore.strand = FALSE)["TRUE"]
  if (annotatedIntronNumberNew > annotatedIntronNumber & !is.na(
    annotatedIntronNumber)) {
    # update junctions object if strand prediction improves overlap
    # with annotations
    if (verbose) {
      message("after strand correction, annotated introns:")
      message(annotatedIntronNumberNew)
      message(annotatedIntronNumberNew / length(uniqueJunctionsUpdate))
    }
    annotatedIntronNumber <- annotatedIntronNumberNew
    uniqueJunctions <- uniqueJunctionsUpdate
  } else {
    strandStep <- FALSE
  }
  outputList <- list(
    "strandStep" = strandStep,
    "annotatedIntronNumber" = annotatedIntronNumber,
    "uniqueJunctions" = uniqueJunctions
  )
  return(outputList)    
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

