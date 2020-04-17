#' @noRd
transcriptToGeneExpression<- function(se, annotationGrangesList){
  counts <- as.data.table(assays(se)$estimates,keep.rownames = TRUE)
  counts <- melt(counts, id.vars = "rn", measure.vars = colnames(counts)[-1])
  setnames(counts, "rn","TXNAME")

  counts <- as.data.table(rowData(se))[,.(TXNAME,GENEID)][counts, on = "TXNAME"]

  counts[, valueGene:=sum(value), by = list(variable, GENEID)]
  counts[, valueGeneCPM:=valueGene/sum(value)*10^6, by = list(variable)]

  ## counts
  counts_gene <- dcast(unique(counts[,.(GENEID, variable, valueGene)]), GENEID ~ variable, value.var = "valueGene")
  counts_gene_CPM <- dcast(unique(counts[,.(GENEID, variable, valueGeneCPM)]), GENEID ~ variable, value.var = "valueGeneCPM")

  ## geneRanges
  gene_tx_map <- unlist(annotationGrangesList)
  elementData <- as.data.table(elementMetadata(annotationGrangesList))[,.(TXNAME,GENEID,newTxClass)]
  elementData[, newGeneClass := ifelse(grepl("ENSG",GENEID),"annotation",unique(newTxClass)), by = GENEID]

  tmp <- data.table(TXNAME = names(gene_tx_map),
                    start = start(gene_tx_map),
                    end = end(gene_tx_map),
                    strand = as.character(strand(gene_tx_map)),
                    seqnames = as.character(seqnames(gene_tx_map)))
  tmp <- elementData[,.(TXNAME,GENEID,newGeneClass)][tmp, on = "TXNAME"]

  gene_tmp <- tmp[, list(seqnames = unique(seqnames),
                         start = min(start),
                         end = max(end),
                         strand = unique(strand)), by = list(GENEID, newGeneClass)]

  gene_ranges <- GRanges(seqnames = Rle(gene_tmp$seqnames),
                         ranges = IRanges(gene_tmp$start, end = gene_tmp$end, names = gene_tmp$GENEID),
                         strand = Rle(strand(gene_tmp$strand)))
  mcols(gene_ranges)$newGeneClass <- gene_tmp$newGeneClass



  ## SE
  counts_gene <- setDF(counts_gene)
  RowNames <- counts_gene$GENEID
  rownames(counts_gene) <- RowNames
  counts_gene_CPM <- setDF(counts_gene_CPM)
  rownames(counts_gene_CPM) <- RowNames
  ColNames <- colnames(counts_gene)[-1]

  seOutput <- SummarizedExperiment(assays = SimpleList(counts = as.matrix(counts_gene[,-1],
                                                                ncol = length(ColNames),
                                                                dimnames = list(RowNames, ColNames)),
                                                      CPM = as.matrix(counts_gene_CPM[match(RowNames, counts_gene_CPM$GENEID),-1],
                                                            ncol =  length(ColNames),
                                                            dimnames = list(RowNames, ColNames))),
                                                     rowRanges = gene_ranges[RowNames],
                                                     colData = colData(se))
  return(seOutput)
}
