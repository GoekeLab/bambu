#' CheckReadClassTxAssignmentUniqueness
#'@title
#'@param mapping
#'@param nobs
aggReadClass <- function(dt){
  dt[, tx_p:=paste(sort(unique(tx_sid)), collapse = '.'), by = list(read_class_sid,gene_sid)]
  dt[, read_class_sid_stored:=read_class_sid]
  equiClassVec <- unique(dt$equiClass_p)
  dt[, read_class_sid:=match(equiClass_p,equiClassVec)]
  dt[, nobs_stored:=nobs]
  dt[, nobs:=NULL]
  tmp <- unique(dt[,.(read_class_sid, read_class_sid_stored, nobs_stored)])
  tmp[, nobs:=sum(nobs_stored), by = read_class_sid]
  tmp <- unique(tmp[,.(read_class_id, nobs)])
  dt <- unique(dt[,.(read_class_sid, tx_sid,gene_sid)])[tmp, on = 'read_class_sid']
  return(list(dt,equiClassVec))
}
