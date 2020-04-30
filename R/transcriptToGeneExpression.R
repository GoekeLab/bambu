#' Reduce transcript expression to gene expression
#' @noRd
transcriptToGeneExpression<- function(se){
  counts <- as.data.table(assays(se)$counts,keep.rownames = TRUE)
  counts <- melt(counts, id.vars = "rn", measure.vars = colnames(counts)[-1])
  setnames(counts, "rn","TXNAME")

  rowDataSe <- as.data.table(rowData(se))

  counts <- rowDataSe[,.(TXNAME,GENEID)][counts, on = "TXNAME"]

  counts[, valueGene:=sum(value), by = list(variable, GENEID)]
  counts[, valueGeneCPM:=valueGene/sum(value)*10^6, by = list(variable)]

  ## counts
  counts_gene <- dcast(unique(counts[,.(GENEID, variable, valueGene)]), GENEID ~ variable, value.var = "valueGene")
  counts_gene_CPM <- dcast(unique(counts[,.(GENEID, variable, valueGeneCPM)]), GENEID ~ variable, value.var = "valueGeneCPM")

  ## geneRanges
  exByGene <- txRangesToGeneRanges(rowRanges(se),TXNAMEGENEID_Map = rowDataSe[,.(TXNAME,GENEID)])


  if("newTxClass" %in% colnames(rowDataSe)){
    rowDataSe <- rowDataSe[,.(TXNAME,GENEID,newTxClass)]
    rowDataSe[, newGeneClass := ifelse(grepl("ENSG",GENEID),"annotation",unique(newTxClass)), by = GENEID]
    mcols(exByGene) <- unique(rowDataSe[,.(GENEID,newGeneClass)])[match(names(exByGene),GENEID)]
  }




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
                                                     rowRanges = exByGene[RowNames],
                                                     colData = colData(se))
  return(seOutput)
}
