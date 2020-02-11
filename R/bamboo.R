#' Transcript abundance quantification and isoform recontruction
#'@title Transcript abundance quantification and isoform recontruction
#'@param obj  A named vector of observed counts for each read class
#'@param annotation A \code{.data.frame} object that maps read class to transcript, with first column being read class name or ids, second column being transcript id or names
#'@param col_indexes A vector that indicates the columns for read class, transcript, and gene.
#'@export
#'@examples
#' \dontrun{
#'  test.bam <- system.file("extdata", "GIS_HepG2_cDNAStranded_Rep5-Run4_chr9_108865774_109014097.bam", package = "bamboo")
#'  fa.file <- system.file("extdata", "Homo_sapiens.GRCh38.dna_sm.primary_assembly_chr9.fa.gz", package = "bamboo")
#'  bamboo(obj = test.bam, fa.file = fa.file)
#'  }
bamboo <- function(obj,  algo.control,...){
  input.list <- list(obj,  algo.control)
  ##----step1: check input arguments and preprocess input arguments
  input.list <- checkInput(input.list)

  dt <- input.list[[1]]
  algo.control <- input.list[[2]]
  ##----step2: match to simple numbers to increase claculation efficiency
  geneVec <- unique(dt$gene_id)
  txVec <- unique(dt$tx_id)
  readclassVec <- unique(dt$read_class_id)
  dt[, gene_sid:=match(gene_id, geneVecTotal)]
  dt[, tx_sid:=match(tx_id, txVecTotal)]
  dt[, read_class_sid:=match(tx_id, txVecTotal)]


  ##----step3: aggregate read class
  temp <- aggReadClass(dt)
  dt <- temp[[1]]
  equiClassVec <- temp[[2]]

  ##----step4: quantification
  start.time <- proc.time()
  outList <- abundance_quantification(dt,
                                      mc.cores = algo.control[["ncore"]],
                                      method = algo.control[["method"]],
                                      conv.control = algo.control[["convcontrol"]])
  end.time <- proc.time()
  cat(paste0('Finished EM estimation in ', round((end.time-start.time)[3]/60,1), ' mins', ' \n'))


  theta_est <- outList[[1]]
  theta_est[, `:=`(tx_name = txVecTotal[as.numeric(tx_sid)],
                   gene_name = geneVecTotal[gene_sid] )]
  theta_est[,`:=`(tx_sid=NULL, gene_sid=NULL)]


  b_est <- outList[[2]]
  b_est[, `:=`(gene_name = geneVecTotal[gene_sid], equiClass = equiClassVec[read_class_sid])]
  b_est[, `:=`(gene_sid = NULL,read_class_sid=NULL)]




  theta_est <- theta_est[,.(tx_name, estimates)]
  theta_est[,`:=`(CPM = estimates/sum(estimates)*(10^6))]
  theta_est <- setDF(theta_est)
  rownames(theta_est) <- theta_est$tx_name
  theta_est$tx_name <- NULL


  b_est[, nobs:=NULL]
  b_est <- setDF(b_est)

  est.list <- list(counts = theta_est, metadata = b_est)
  return(est.list)
}

setMethod("bamboo", signature(obj = "data.table"),
          definition = function(obj){
            bamboo(obj,...)
          })
## method when output from buildTranscriptModel is provided
setMethod("bamboo", signature(obj = "summarizedExperiment"),
          definiton = function(obj){
            dt <- process_se(obj)
            rowData <- metadata(dt)[,c("TXNAME","GeneID","eqClass")]
            rownames(rowData) <- rowData$TXNAME
            rowData$TXNAME <- NULL
            ## To do:
            ## task1: optional: to implement filtering function
            ## task2: to implement for multiple samples, when multiple samples are provided, run txdbtableslist for one time
            est <- bamboo(dt)
            counts <- est[["counts"]]
            seOutput <- summarizedExperiment(assays = list(counts = counts[match(rownames(rowData),rownames(counts)),]),
                                             rowData = rowData,
                                             metadata = est[["metadata"]]) # transcript annotation with read class information
            return(seOutput)
          })



setMethod("bamboo", signature(obj = "character"),
          function(obj){
            if(is.null(txdbTablesList)){
              if(is.null(txdb)){
                stop("txdb object is missing!")
              }else{
                txdbTablesList <- prepareAnnotations(txdb)
              }
            }
            if(is.null(fa.file)){
              stop("Genome fa file is missing!")
            }else{
              if(class(fa.file)!='FaFile') {
                fa.file <- FaFile(fa.file)
              }else{
                stop('Genome Fa is not a fa file')
              }
            }
            if(is.null(ir.control)){
              ir.control <- list(stranded = FALSE,
                                 protocol = NULL,
                                 prefix = '',
                                 minimumReadSupport = 2,
                                 minimumTxFraction = 0.02)
            }else{
              if(is.null(ir.control[['stranded']])){
                ir.control[['stranded']] <- FALSE
              }
              if(is.null(ir.control[['prefix']])){
                ir.control[['prefix']] <- ''
              }
              if(is.null(ir.control[['minimumReadSupport']])){
                ir.control[['minimumReadSupport']] <- 2
              }
              if(ir.control[['minimumTxFraction']]){
                ir.control[['minimumTxFraction']] <- 0.02
              }
            }
            start.time <- proc.time()
            se  <- isore(obj,txdbTablesList = txdbTablesList,
                                         genomeFA = fa.file,
                                         stranded = ir.control[['stranded']],
                                         protocol = ir.control[['protocol']],
                                         prefix = ir.control[['prefix']],
                                         minimumReadSupport= ir.control[['minimumReadSupport']],
                                         minimumTxFraction = ir.control[['minimumTxFraction']])
            end.time <- proc.time()
            cat(paste0('Finished build transcript models in ', round((end.time-start.time)[3]/60,1), ' mins', ' \n'))

            seOutput <- bamboo(se)
            return(seOutput)
          }
)
