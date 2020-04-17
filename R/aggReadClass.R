#' CheckReadClassTxAssignmentUniqueness
#' @param dt A data.table object
#' @noRd
aggReadClass <- function(dt){
  dt[, eqClass:=paste(sort(unique(tx_sid)), collapse = '.'), by = list(read_class_sid,gene_sid)]
  dt[, read_class_sid_stored:=read_class_sid]
  eqClassVec <- unique(dt$eqClass)
  dt[, read_class_sid:=match(eqClass,eqClassVec)]
  dt[, nobs_stored:=nobs]
  dt[, nobs:=NULL]
  tmp <- unique(dt[,.(read_class_sid, read_class_sid_stored, nobs_stored)])
  tmp[, nobs:=sum(nobs_stored), by = read_class_sid]
  tmp <- unique(tmp[,.(read_class_sid, nobs)])
  dt <- unique(dt[,.(read_class_sid, tx_sid,gene_sid)])[tmp, on = 'read_class_sid']
  return(list(dt,eqClassVec))
}
