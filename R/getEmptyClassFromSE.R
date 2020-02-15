getEmptyClassFromSE <- function(se = se, txdbTablesList = NULL){
  distTable <- data.table(metadata(se)$distTable)[,.(readClassId, annotationTxId, readCount, GENEID)]
  distTable[, eqClass:=paste(sort(unique(annotationTxId)),collapse='.'), by = list(readClassId,GENEID)]

  eqClassCountTable <- distTable[,list(readCount=sum(readCount)), by = list(eqClass, GENEID)]
  eqClassTable <- data.table(txdbTablesList$txIdToGeneIdTable)
  eqClassCountTable <- eqClassCountTable[, list(TXNAME=unlist(strsplit(eqClass,'\\.'))), by = list(eqClass,GENEID, readCount)]
  #check that if observed classes is more than reference eqClasses

  eqClassCountTable <- unique(merge(eqClassCountTable,eqClassTable,all = TRUE, on = c('GENEID','eqClass','TXNAME'))) # merge should be performed on both sides
  # if there exists new isoforms from eqClassCountTable, it would not be found in eqClassTable, keep them
  eqClassCountTable[is.na(readCount), readCount:=0]

  ## remove empty read class where there is no shared class found
  eqClassCountTable[,sum_nobs:=sum(readCount), by = list(GENEID, TXNAME)]

  eqClassCountTable <- unique(eqClassCountTable[sum_nobs>0,.(GENEID, eqClass, readCount,TXNAME)])

  setnames(eqClassCountTable,old = c("TXNAME","GENEID","eqClass","readCount") , new = c("tx_id","gene_id","read_class_id","nobs"))


  return(eqClassCountTable)

}
