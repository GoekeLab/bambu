#' Reduce transcript expression to gene expression
#' @title transcript to gene expression
#' @param se a summarizedExperiment object from \code{\link{bambu}}
#' @return A SummarizedExperiment object
#' @import data.table 
#' @export
#' @examples
#' se <- readRDS(system.file("extdata",
#'     "seOutput_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.rds",
#'     package = "bambu"
#' ))
#' transcriptToGeneExpression(se)
transcriptToGeneExpression <- function(se) {
    counts <- assays(se)$counts
    runnames <- colnames(counts)[-1]
    rowDataSe <- as.data.table(rowData(se))
    
    counts  = fac2sparse(factor(rowData(se)$GENEID, levels = unique(rowData(se)$GENEID))) %*% counts
    if(!is.null(metadata(se)$incompatibleCounts)){
        incompatibleCounts <- metadata(se)$incompatibleCounts
        if("nonuniqueCounts" %in% names(metadata(se))){
            incompatibleCounts = incompatibleCounts + metadata(se)$nonuniqueCounts
        }
        incompatibleCounts = Matrix(incompatibleCounts[match(rownames(counts), rownames(incompatibleCounts)),], sparse = TRUE)
        counts = counts + incompatibleCounts
    }
    counts.total = colSums(counts)
    counts.total[counts.total==0] = 1
    counts.CPM = counts/counts.total * 10^6

    ## geneRanges
    exByGene <- reducedRangesByGenes(rowRanges(se))
    if ("txClassDescription" %in% colnames(rowDataSe)) {
        rowDataSe <- rowDataSe[, .(TXNAME, GENEID, txClassDescription)]
        rowDataSe[, newGeneClass := ifelse(grepl("ENSG", GENEID),
            "annotation", unique(txClassDescription)), by = GENEID]
        mcols(exByGene) <- unique(rowDataSe[, .(GENEID,
            newGeneClass)])[match(names(exByGene), GENEID)]
    }
    ## SE
    RowNames <- rownames(counts)
    ColNames <- colnames(counts)
    ColData <- colData(se)
    ColData@rownames <- ColNames
    ColData@listData$name <- ColNames
    seOutput <- SummarizedExperiment(
    assays = SimpleList(counts = counts,
            CPM = counts.CPM),
        rowRanges = exByGene[RowNames],
        colData = ColData)
    
    return(seOutput)
}

addTssId <- function(se){
  seTssTable <- tibble(start = unlist(endoapply(start(rowRanges(se)), function(x) x[1])), 
                       end = unlist(endoapply(end(rowRanges(se)), function(x) x[length(x)])),
                       strand = as.character(getStrandFromGrList(rowRanges(se))),
                       chr = as.character(getChrFromGrList(rowRanges(se)))) %>%
    mutate(tssRanges = ifelse(strand != "-", start, end)) %>%
    group_by(chr, strand, tssRanges) %>% 
    mutate(TSSID = paste0("BambuTss", cur_group_id())) %>%
    ungroup()
  rowData(se)$TSSID <- seTssTable$TSSID
  rowData(se)$tssRanges <- seTssTable$tssRanges
  return(se)
}

transcriptToTssTesExpression <- function(se, feature = "tss") {
  counts <- assays(se)$counts
  if(feature == "tss"){
      counts  <- fac2sparse(rowData(se)$globleTssId) %*% counts
      ranges <- metadata(rowRanges(se))$tss_clusters
  }
  if(feature == "tes"){
      counts  <- fac2sparse(rowData(se)$globleTesId) %*% counts
      ranges <- metadata(rowRanges(se))$tes_clusters
  }
  RowNames <- rownames(counts)
  ColNames <- colnames(counts)
  ColData <- colData(se)
  ColData@rownames <- ColNames
  ColData@listData$name <- ColNames
  seOutput <- SummarizedExperiment(
    assays = SimpleList(counts = counts),
    rowRanges = ranges[RowNames],
    colData = ColData)
  rowRanges(seOutput) <- split(rowRanges(seOutput), names(rowRanges(seOutput)))
  #only output the TES with read support!
  seOutput <- seOutput[rowSums(assays(seOutput)$counts) > 0, ]
  return(seOutput)
}
