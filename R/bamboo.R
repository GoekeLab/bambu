#' @rdname bamboo
#' @export
setGeneric("bamboo", function(object,...) standardGeneric("bamboo"))

#' Transcript abundance quantification and isoform recontruction
#' @title Transcript abundance quantification and isoform recontruction
#' @param obj  A named vector of observed counts for each read class
#' @param annotation A \code{.data.frame} object that maps read class to transcript, with first column being read class name or ids, second column being transcript id or names
#' @param col_indexes A vector that indicates the columns for read class, transcript, and gene.
#' @useDynLib bamboo
#' @importFrom Rcpp sourceCpp
#' @export
#' @examples
#' \dontrun{
#'  test.bam <- system.file("extdata", "GIS_HepG2_cDNAStranded_Rep5-Run4_chr9_108865774_109014097.bam", package = "bamboo")
#'  fa.file <- system.file("extdata", "Homo_sapiens.GRCh38.dna_sm.primary_assembly_chr9.fa.gz", package = "bamboo")
#'  bamboo(obj = test.bam, fa.file = fa.file)
#'  }
bamboo <- function(object,algo.control = NULL,...){
  if(is.null(object)){
    stop("Input object is missing.")
  }else if(any(!(c('gene_id','tx_id','read_class_id','nobs') %in% colnames(object)))){
    stop("Columns gene_id, tx_id, read_class_id, nobs, are missing from object.")
  }
  dt <- object

  ## check quantification parameters
   default.algo.control <- list(ncore = parallel::detectCores(),
                               method = "two-step",
                               convcontrol = 10^(-3))

  if(is.null(algo.control)){
    algo.control <- default.algo.control
  }else{
    if(is.null(algo.control[["ncore"]])){
      algo.control[["ncore"]] <- default.algo.control[["ncore"]]
    }else if(as.numeric(algo.control[["ncore"]])>default.algo.control[["ncore"]]){
      algo.control[["ncore"]] <- default.algo.control[["ncore"]]
    }
    if(is.null(algo.control[["method"]])){
      algo.control[["method"]] <- default.algo.control[["method"]]
    }else if(!(algo.control[["method"]] %in% c("one-step","two-step"))){
      algo.control[["method"]] <- default.algo.control[["method"]]
    }
    if(is.null(algo.control[["convcontrol"]])){
      algo.control[["convcontrol"]] <- default.algo.control[["convcontrol"]]
    }
  }

  ##----step2: match to simple numbers to increase claculation efficiency
  geneVec <- unique(dt$gene_id)
  txVec <- unique(dt$tx_id)
  readclassVec <- unique(dt$read_class_id)
  dt[, gene_sid:=match(gene_id, geneVec)]
  dt[, tx_sid:=match(tx_id, txVec)]
  dt[, read_class_sid:=match(read_class_id, readclassVec)]

  dt[,`:=`(tx_id = NULL, gene_id = NULL, read_class_id = NULL)]

  ##----step3: aggregate read class
  temp <- aggReadClass(dt)
  dt <- temp[[1]]
  eqClassVec <- temp[[2]]

  ##----step4: quantification
  start.time <- proc.time()
  outList <- abundance_quantification(dt,
                                      mc.cores = algo.control[["ncore"]],
                                      method = algo.control[["method"]],
                                      conv.control = algo.control[["convcontrol"]])
  end.time <- proc.time()
  cat(paste0('Finished EM estimation in ', round((end.time-start.time)[3]/60,1), ' mins', ' \n'))


  theta_est <- outList[[1]]
  theta_est[, `:=`(tx_name = txVec[as.numeric(tx_sid)],
                   gene_name = geneVec[gene_sid] )]
  theta_est[,`:=`(tx_sid=NULL, gene_sid=NULL)]


  b_est <- outList[[2]]
  b_est[, `:=`(gene_name = geneVec[gene_sid], eqClass = eqClassVec[as.numeric(read_class_sid)])]
  b_est[, `:=`(gene_sid = NULL,read_class_sid=NULL)]




  theta_est <- theta_est[,.(tx_name, estimates)]
  theta_est[,`:=`(CPM = estimates/sum(estimates)*(10^6))]
  theta_est <- setDF(theta_est)
  rownames(theta_est) <- theta_est$tx_name
  theta_est$tx_name <- NULL

  b_est <- setDF(b_est)

  est.list <- list(counts = theta_est, metadata = b_est)
  return(est.list)
}

# setMethod("bamboo", signature("data.table"),
#           function(object){
#             bamboo(object,...)
#           })
## method when output from buildTranscriptModel is provided

se.bamboo <- function(object,algo.control = NULL){
  dt <- process_se(object)
  rowData <- metadata(dt)[,c("tx_id","gene_id","read_class_id")]
  rownames(rowData) <- rowData$tx_id
  rowData$TXNAME <- NULL
  ## To do:
  ## task1: optional: to implement filtering function
  ## task2: to implement for multiple samples, when multiple samples are provided, run txdbtableslist for one time
  est <- bamboo(dt,algo.control)
  counts <- est[["counts"]]
  seOutput <- summarizedExperiment(assays = list(counts = counts[match(rownames(rowData),rownames(counts)),]),
                                   rowData = rowData,
                                   metadata = est[["metadata"]]) # transcript annotation with read class information
  return(seOutput)
}


setOldClass("summarizedExperiment")

#' Accessors for the 'counts' slot of a DESeqDataSet object.
#'
#' The counts slot holds the count data as a matrix of non-negative integer
#' count values, one row for each observational unit (gene or the like), and one
#' column for each sample.
#'
#' @docType methods
#' @name counts
#' @rdname counts
#' @aliases counts counts,DESeqDataSet-method counts<-,DESeqDataSet,matrix-method
#'
#' @param object a \code{DESeqDataSet} object.
#' @param normalized logical indicating whether or not to divide the counts by
#' the size factors or normalization factors before returning
#' (normalization factors always preempt size factors)
#' @param replaced after a \code{DESeq} call, this argument will return the counts
#' with outliers replaced instead of the original counts, and optionally \code{normalized}.
#' The replaced counts are stored by \code{DESeq} in \code{assays(object)[['replaceCounts']]}.
#' @param value an integer matrix
#' @author Simon Anders
#' @seealso \code{\link{sizeFactors}}, \code{\link{normalizationFactors}}
#'
#' @examples
#'
#' dds <- makeExampleDESeqDataSet(m=4)
#' head(counts(dds))
#'
#' dds <- estimateSizeFactors(dds) # run this or DESeq() first
#' head(counts(dds, normalized=TRUE))
#'
#' @export
setMethod("bamboo", signature(object = "summarizedExperiment"),se.bamboo)




bam.bamboo <- function(object,algo.control = NULL, fa.file=NULL, txdb=NULL, txdbTablesList=NULL){
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
      fa.file <- Rsamtools::FaFile(fa.file)
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

  seOutput <- bamboo(se, algo.control)
  return(seOutput)
}


#' Accessors for the 'counts' slot of a DESeqDataSet object.
#'
#' The counts slot holds the count data as a matrix of non-negative integer
#' count values, one row for each observational unit (gene or the like), and one
#' column for each sample.
#'
#' @docType methods
#' @name counts
#' @rdname counts
#' @aliases counts counts,DESeqDataSet-method counts<-,DESeqDataSet,matrix-method
#'
#' @param object a \code{DESeqDataSet} object.
#' @param normalized logical indicating whether or not to divide the counts by
#' the size factors or normalization factors before returning
#' (normalization factors always preempt size factors)
#' @param replaced after a \code{DESeq} call, this argument will return the counts
#' with outliers replaced instead of the original counts, and optionally \code{normalized}.
#' The replaced counts are stored by \code{DESeq} in \code{assays(object)[['replaceCounts']]}.
#' @param value an integer matrix
#' @author Simon Anders
#' @seealso \code{\link{sizeFactors}}, \code{\link{normalizationFactors}}
#'
#' @examples
#'
#' dds <- makeExampleDESeqDataSet(m=4)
#' head(counts(dds))
#'
#' dds <- estimateSizeFactors(dds) # run this or DESeq() first
#' head(counts(dds, normalized=TRUE))
#'
#' @export
setMethod("bamboo", signature("character"), bam.bamboo)

