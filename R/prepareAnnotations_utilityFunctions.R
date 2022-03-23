#' Prepare annotation granges object from GTF file
#' @title Prepare annotation granges object from GTF file into a 
#' GRangesList object
#' @param file a GTF file
#' @return A \code{\link{GRangesList}} object
#' @details Unlike \code{\link{readFromGTF}}, this function finds out the
#' equivalence classes between the transcripts,
#' with \code{\link{mcols}} data having three columns:
#' \itemize{
#'     \item TXNAME specifying prefix for new gene Ids (genePrefix.number),
#'     defaults to empty
#'     \item GENEID indicating whether filter to remove read classes
#'     which are a subset of known transcripts(), defaults to TRUE
#'     \item eqClass specifying minimun read count to consider a read class
#'     valid in a sample, defaults to 2
#'   }
#' @importFrom GenomicRanges makeGRangesListFromDataFrame 
#' @noRd
prepareAnnotationsFromGTF <- function(file) {
    if (missing(file)) {
        stop("A GTF file is required.")
    } else {
        data <- utils::read.delim(file, header = FALSE, comment.char = "#")
        colnames(data) <- c("seqname", "source", "type", "start", "end",
            "score", "strand", "frame", "attribute")
        data <- data[data$type == "exon", ]
        data$strand[data$strand == "."] <- "*"
        data$GENEID <- gsub("gene_id (.*?);.*", "\\1", data$attribute)
        data$TXNAME <- gsub(".*transcript_id (.*?);.*", "\\1", data$attribute)
        multiTxCheck <- as_tibble(data) %>% select(seqname, GENEID) %>% distinct() %>% group_by(GENEID) %>% 
            mutate(n=n(), id=paste0('-',row_number()))
        if(any(multiTxCheck$n>1)) { # identical TXNAMES
            warning('Annotations contain duplicated transcript names
                        Transcript names will be made unique')
            uniqueNamesTbl <- as_tibble(data) %>% 
                select(seqname, TXNAME, GENEID) %>% 
                mutate(order=row_number()) %>% left_join(multiTxCheck) %>%
               mutate(gene_unique = paste0(GENEID, ifelse(n==1,'', id)),
                      tx_unique = paste0(TXNAME, ifelse(n==1,'', id))) %>%
                arrange(order)
            
            data$TXNAME_Original <- data$TXNAME
            data$GENEID_Original <- data$GENEID
            data$TXNAME <- uniqueNamesTbl$tx_unique
            data$GENEID <- uniqueNamesTbl$gene_unique
            }
        geneData <- unique(data[, c("TXNAME", "GENEID")])
        grlist <- makeGRangesListFromDataFrame(
        data[, c("seqname", "start", "end", "strand", "TXNAME")],
            split.field = "TXNAME", keep.extra.columns = TRUE)
        grlist <- grlist[IRanges::order(start(grlist))]
        unlistedExons <- unlist(grlist, use.names = FALSE)
        partitioning <- PartitioningByEnd(cumsum(elementNROWS(grlist)),
            names = NULL)
        txIdForReorder <- togroup(PartitioningByWidth(grlist))
        exon_rank <- lapply(elementNROWS(grlist), seq, from = 1)
        exon_rank[which(unlist(unique(strand(grlist))) == "-")] <- lapply(
            exon_rank[which(unlist(unique(strand(grlist))) == "-")], rev
            ) # * assumes positive for exon ranking
        names(exon_rank) <- NULL
        unlistedExons$exon_rank <- unlist(exon_rank)
        unlistedExons <- unlistedExons[order(txIdForReorder,
            unlistedExons$exon_rank)]
        # exonsByTx is always sorted by exon rank, not by strand,
        # make sure that this is the case here
        unlistedExons$exon_endRank <- unlist(lapply(elementNROWS(grlist),
            seq, to = 1), use.names = FALSE)
        unlistedExons <- unlistedExons[order(txIdForReorder,
            start(unlistedExons))]
        mcols(unlistedExons) <- mcols(unlistedExons)[, c("exon_rank",
            "exon_endRank")]
        grlist <- relist(unlistedExons, partitioning)
        # sort the grlist by start position, ranked by exon number
        mcols(grlist) <- DataFrame(geneData[(match(names(grlist),
            geneData$TXNAME)), ])
        mcols(grlist)$txid <- seq_along(grlist)
        minEqClasses <- getMinimumEqClassByTx(grlist)
        if(!identical(names(grlist),minEqClasses$queryTxId)) warning('eq classes might be incorrect')
        mcols(grlist)$eqClass <- minEqClasses$eqClass
        mcols(grlist)$eqClassById <- minEqClasses$eqClassById
    }
    return(grlist)
}


#' Get minimum equivalent class by Transcript
#' equivalent classes are identical if transcripts are ordered alphabetically,
#' annotations=prepareAnnotations('annotations.gtf')
#' eqList <- unstrsplit(splitAsList(mcols(annotations)$TXNAME[unlist(mcols(annotations)$eqClassById)],
#' rep(mcols(annotations)$TXNAME, elementNROWS(mcols(annotations)$eqClassById))), sep='.')
#' all(mcols(annotations)eqClass == eqList)
#' 
#' @param exonsByTranscripts exonsByTranscripts
#' @importFrom dplyr tibble
#' @noRd
getMinimumEqClassByTx <- function(exonsByTranscripts) {
    exByTxAnnotated_singleBpStartEnd <-
        cutStartEndFromGrangesList(exonsByTranscripts)
    # estimate overlap only based on junctions
    spliceOverlaps <- findSpliceOverlapsQuick(
        exByTxAnnotated_singleBpStartEnd,
        exByTxAnnotated_singleBpStartEnd
        )
    ## identify transcripts compatible with other (subsets by splice sites)
    spliceOverlaps <- spliceOverlaps[mcols(spliceOverlaps)$compatible == TRUE, ]
    ## select splicing compatible transcript matches
    if(identical(mcols(exonsByTranscripts)$txid, seq_along(exonsByTranscripts))) {
        queryTxId <- queryHits(spliceOverlaps)
        subjectTxId <- subjectHits(spliceOverlaps)
        if(! identical(unique(queryTxId), mcols(exonsByTranscripts)$txid)) warning('eq classes might be incorrect')
        eqClassById <- unname(splitAsList(subjectTxId[order(queryTxId, subjectTxId)], queryTxId))
    } else {
        eqClassById <- rep(NA, length(exonsByTranscripts)) 
        }
    queryTxId <-
        names(exByTxAnnotated_singleBpStartEnd)[queryHits(spliceOverlaps)]
    subjectTxId <-
        names(exByTxAnnotated_singleBpStartEnd)[subjectHits(spliceOverlaps)]
    subjectTxId <- subjectTxId[order(queryTxId, subjectTxId)]
    queryTxId <- sort(queryTxId)
    eqClass <- unstrsplit(splitAsList(subjectTxId, queryTxId), sep = ".")
    eqClass <- eqClass[match(names(exonsByTranscripts), names(eqClass))]

    return(DataFrame(queryTxId = names(eqClass), eqClass = unname(eqClass), eqClassById=eqClassById))
}

