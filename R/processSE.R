process_se <- function(se){
  obj <- assay(se)
  annotation <- metadata(se)[,c("TXNAME","GENEID","eqClass")]
  annotation <- data.table(annotation)
  setnames(annotation,old = c("TXNAME","GENEID","eqClass") , new = c("tx_id","gene_id","read_class_id"))
  obj <- data.table(read_class_id = colnames(obj), nobs = as.numeric(obj))

  if(!all(unique(obj$read_class_id) %in% unique(annotation$read_class_id))){
    stop("The names of obj is not compatible with that provided in annotation.")
  }
  dt <- obj[annotation, on = 'read_class_id']
  dt[is.na(nobs), nobs:=0]
  return(dt)
}






