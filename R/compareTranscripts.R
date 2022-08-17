### examples for test purposes
## Expected annotations of transcripts used in test query
# 'ENST00000344579', # exon skipping, alternative TSS (-48), +, ENSG00000158109
# 'ENST00000270792', # intron retention subject 1(last exon),alt.TSS,alt.TES, +,
# 'ENST00000410032', # alternative first exon, exon skipping query: 2, 
#                    # exon skipping subject: 0, alternative TES (2bp only), 
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

#' @title compare alternatively-spliced transcripts
#' @param query a GRangesList of transcripts
#' @param subject a GRangesList of transcripts
#' @details This function compares two alternatively-spliced transcripts and returns 
#' a \code{tibble} object that determines the type of the alternative splicing between 
#' the query and the subject transcript. Alternative splicing includes exons skipping, 
#' intron retention, alternative 5' exon start site, alternative 3' exon end site, alternative
#' first exon, alternative last exon, alternative transcription start site (TSS), 
#' alternative transcription end site (TES), internal first exon, internal last exon, or a 
#' mixed combination of them. 
#' @seealso \code{Value} for more details about each of the alternative splicing events. 
#' @return a \code{tibble} object with the following columns:
#' 
#' Remark: two exons, one each from query and subject transcript are defined to be equivalent if one of 
#' the exons fully covers the another exon. 
#' 
#' \itemize{
#'     \item alternativeFirstExon: FALSE if query transcript and subject transcript have equivalent
#'     first exon. TRUE otherwise. 
#'     \item alternativeTSS: +k if the query initiates the transcription k sites earlier 
#'     than the subject transcript, -k if the query initiates the transcription k sites later 
#'     than the subject transcript. 
#'     \item internalFirstExon.query: TRUE if the first exon of the query transcript is equivalent to
#'     one of the exon of the subject transcript (except the first exon). FALSE otherwise. 
#'     \item internalFirstExon.subject: TRUE if the first exon of the subject transcript is equivalent to
#'     one of the exon of the query transcript (except the first exon). FALSE otherwise. 
#'     \item alternativeLastExon: FALSE if query transcript and subject transcript have equivalent last
#'     exon. TRUE otherwise. 
#'     \item alternativeTES: +k if the query transcript ends the transcription k sites later than the 
#'     subject transcript, -k if the query transcript ends the transcription k sites earlier than the 
#'     subject transcript. 
#'     \item internalLastExon.query: TRUE if the last exon of the query transcript is equivalent to one 
#'     of the exon in the subject transcript (except the last exon). FALSE otherwise. 
#'     \item internalLastExon.subject: TRUE if the last exon of the subject transcript is equivalent to 
#'     one of the exon in the query transcript (except the last exon). FALSE otherwise. 
#'     \item intronRetention.subject: k if there is/are k intron(s) in the subject transcript relative to the 
#'     query transcript. 
#'     \item intronRetention.query: k if there is/are k intron(s) in the query transcript relative to the 
#'     subject transcript. 
#'     \item exonSkipping.query: k if there is/are k exon(s) in the subject transcript (except the first and last) not in the
#'     query transcript. 
#'     \item exonSkipping.subject: k if there is/are k exon(s) in the query transcript (except the first and last) not in the 
#'     subject transcript. 
#'     \item exon5prime (splicing): k if there is/are k equivalent exon(s) having different 5' start site (except the 5' start site of first exon).
#'     \item exon3prime (splicing): k if there is/are k equivalent exon(s) having different 3' end site (except the 3' end site of the last exon).
#' }
#' @importFrom dplyr tibble 
#' @examples
#' query <- readRDS(system.file("extdata", 
#'     "annotateSpliceOverlapByDist_testQuery.rds",
#'     package = "bambu"))
#' subject <- readRDS(system.file("extdata", 
#'     "annotateSpliceOverlapByDist_testSubject.rds",
#'     package = "bambu"))
#' compareTranscriptsTable <- compareTranscripts(query, subject)
#' 
#' compareTranscriptsTable
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
