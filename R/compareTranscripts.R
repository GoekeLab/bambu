### examples for test purposes
## Expected annotations of transcripts used in test query
# 'ENST00000344579', # exon skipping, alternative TSS (-48), +, ENSG00000158109
# 'ENST00000270792', # intron retention subject 1(last exon),alt.TSS,alt.TES, +,
# 'ENST00000410032', # alternative first exon, exon skipping query: 2, 
#                    # exon skipping subject: 0, alternative TSS (2bp only), 
#                    # internalFirstExon.subject +
# 'ENST00000468178', # alternative last exon +
# 'ENST00000485956', # alternative first exon, alternative last exon,
# #exon skipping subject = 1, internal first exon query, +
# 'ENST00000530807', # exon skipping query 1, alternative TSS (-17),  -
# 'ENST00000409894', # alternative 3' exon splice site, exon skipping query 2,
# #alternative TSS, alterantive TES, +, ENSG00000125630
# 'ENST00000524447',  # alternative TSS, alternative last exon (internal), 
# #alternative exon 3' end,-, ENSG00000165916
# 'ENST00000591696' # alternative TSS, alternative 3' exon (2), 
# #alternative 5' exon (1) alternative TES, ,+,ENSG00000141349
############################################################
# query <- readRDS(system.file("extdata", 
#     "annotateSpliceOverlapByDist_testQuery.rds",
#     package = "bambu"))
# subject <- readRDS(system.file("extdata", 
#        "annotateSpliceOverlapByDist_testSubject.rds",
#        package = "bambu"))
############################################################

#' annotate splice overlap by distance
#' @description This function takes in a GRangesList (query)
#' and a target GRangesList (subject). The function creates
#' an annotation table in tibble by comparing ranges entries
#' from transcripts between the query and subject GRangesLists.
#' @usage compareTranscripts(query, subject)
#' @params query a GRangesList
#' @params subject a GRangesList
#' @return a tibble with the following annotations:
#' \itemize{
#'    \item alternativeFirstExon
#'    \item alternativeTSS
#'    \item internalFirstExon.query
#'    \item internalFirstExon.subject
#'    \item alternativeLastExon
#'    \item alternativeTES
#'    \item internalLastExon.query
#'    \item internalLastExon.subject
#'    \item intronRetention.subject
#'    \item intronRetention.query
#'    \item exonSkipping.query
#'    \item exonSkipping.subject
#'    \item exon5prime (splicing)
#'    \item exon3prime (splicing)
#' }
#' @examples
#' query <- readRDS(system.file("extdata", 
#'     "annotateSpliceOverlapByDist_testQuery.rds",
#'     package = "bambu"))
#' subject <- readRDS(system.file("extdata", 
#'     "annotateSpliceOverlapByDist_testSubject.rds",
#'     package = "bambu"))
#' annotationTable <- compareTranscripts(query, subject)
#' @noRd
compareTranscripts <- function(query, subject) {
    subjectFullRng <- ranges(subject)
    queryFullRng <- ranges(query)
    strand <- as.character(getStrandFromGrList(query))
    queryStartRng <- selectStartEndExonFromRangesList(queryFullRng, strand,
        "start")
    subjectStartRng <- selectStartEndExonFromRangesList(subjectFullRng, strand,
        "start")
    queryEndRng <- selectStartEndExonFromRangesList(queryFullRng, strand, 
        "end")
    subjectEndRng <- selectStartEndExonFromRangesList(subjectFullRng, strand, 
        "end")
    querySpliceRng <- ranges(myGaps(query))
    querySpliceRng[elementNROWS(querySpliceRng) == 0] <-
        IRanges(start = 1,end = 1) # add mock intron
    subjectSpliceRng <- ranges(myGaps(subject))
    subjectSpliceRng[elementNROWS(subjectSpliceRng) == 0] <- 
        IRanges(start = 1,end = 1)# add mock intron
    annotatedTable <- tibble(queryId = names(query), 
        subjectId = names(subject), strand = strand)
    annotatedTable$alternativeFirstExon <- 
        alternativeStartEndExon(queryStartRng, subjectStartRng)
    annotatedTable$alternativeTSS <- calculateTerminalDistance(queryStartRng, 
        subjectStartRng, annotatedTable$alternativeFirstExon, strand, "start")
    annotatedTable$internalFirstExon.query <- annotateInternalStartEnd(
        queryStartRng, subjectFullRng, annotatedTable$alternativeFirstExon)
    annotatedTable$internalFirstExon.subject <- annotateInternalStartEnd(
        subjectStartRng, queryFullRng, annotatedTable$alternativeFirstExon)
    annotatedTable$alternativeLastExon <-
        alternativeStartEndExon(queryEndRng, subjectEndRng)
    annotatedTable$alternativeTES <- calculateTerminalDistance(queryEndRng, 
        subjectEndRng, annotatedTable$alternativeLastExon, strand, "end")
    annotatedTable$internalLastExon.query <- annotateInternalStartEnd(
        queryEndRng, subjectFullRng, annotatedTable$alternativeLastExon)
    annotatedTable$internalLastExon.subject <- annotateInternalStartEnd(
        subjectEndRng, queryFullRng, annotatedTable$alternativeLastExon)
    annotatedTable$intronRetention.subject <- 
        annotateIntronRetent(querySpliceRng, subjectFullRng)
    annotatedTable$intronRetention.query <- 
        annotateIntronRetent(subjectSpliceRng, queryFullRng)
    annotatedTable$exonSkipping.query <- annotateExonSkip(querySpliceRng, 
        subjectFullRng, subjectStartRng, subjectEndRng)
    annotatedTable$exonSkipping.subject <- annotateExonSkip(subjectSpliceRng, 
        queryFullRng, queryStartRng, queryEndRng) 
    exonSpliceTable <- annotateExonSplice(querySpliceRng, subjectFullRng,
        subjectStartRng, subjectEndRng, annotatedTable$strand)
    annotatedTable <- cbind(annotatedTable, exonSpliceTable)
    return(annotatedTable)
}
