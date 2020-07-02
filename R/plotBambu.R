
#' @noRd
plot <- function(se, ..., group.variable = NULL, type = c("annotation","pca","heatmap"), gene_id = NULL, transcript_id = NULL){
    UseMethod("plot")
   }
#' plotSEOuptut
#' @title plot.bambu
#' @param se An summarized experiment object obtained from \code{\link{bambu}} or \code{\link{transcriptToGeneExpression}}.
#' @param group.variable Variable for grouping in plot, has be to provided if choosing to plot PCA.
#' @param type plot type variable, a values of annotation for a single gene with heatmap for isoform expressions,  pca,  or heatmap, see details.
#' @param gene_id specifying the gene_id for plotting gene annotation, either gene_id or transcript_id has to be provided when type = "annotation".
#' @param transcript_id specifying the transcript_id for plotting transcript annotation, either gene_id or transcript_id has to be provided when type = "annotation"
#' @details \code{\link{type}} indicates the type of plots to be plotted. There are two types of plots can be chosen, PCA or heatmap.
#' @return A heatmap plot for all samples
#' @importFrom stats prcomp
#' @importFrom ggplot2 ggplot
#' @importFrom ComplexHeatmap Heatmap
#' @importFrom circlize colorRamp2
#' @importFrom RColorBrewer brewer.pal
#' @importFrom ggbio autoplot
#' @importFrom gridExtra grid.arrange
#' @export
#' @examples 
#' se <- readRDS(system.file("extdata", "seOutputCombined_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.rds", package = "bambu"))
#' colnames(se) <- colData(se)$name <- c("sample1","sample2")
#' assays(se)$CPM[,2]  <- pmax(0, rnorm(length(assays(se)$CPM[,2]),assays(se)$CPM[,2],10))
#' plot.bambu(se, type = "heatmap")
plot.bambu <- function(se,group.variable = NULL, type = c("annotation","pca","heatmap"), gene_id = NULL, transcript_id = NULL){

  if(type == "annotation"){
    if(is.null(gene_id)&(is.null(transcript_id))){
      stop("Please provide the gene_id(s) of the gene of interest or transcript_id(s) for the transcripts of interest!")
    }
    if(ncol(rowData(se))==0){
      if(is.null(gene_id)){
        stop("Please provide the gene_id(s) of the gene of interest when gene expression are provided!")
      }
      if(!all(gene_id %in% rownames(se))){
        stop("all(gene_id %in% rownames(se)) condition is not satisfied!")
      }
      geneRanges <- rowRanges(se)[gene_id]
      names(geneRanges) <- paste0(gene_id,":", unlist(lapply(strand(geneRanges),function(x) unique(as.character(x)))))
      
      p_annotation <- ggbio::autoplot(geneRanges, group.selfish = TRUE)
      p_expression <- ggbio::autoplot(as.matrix(log2(assays(se)$CPM[gene_id,]+1)))
      
      p <- gridExtra::grid.arrange(p_annotation@ggplot, p_expression)
      return(p)
    }else{
      if(!is.null(transcript_id)){
        if(!all(transcript_id %in% rownames(se))){
          stop("all(transcript_id %in% rownames(se)) condition is not satisfied!")
        }
          txRanges <- rowRanges(se)[transcript_id]
          names(txRanges) <- paste0(transcript_id,":", unlist(lapply(strand(txRanges),function(x) unique(as.character(x)))))
          p_annotation <- ggbio::autoplot(txRanges, group.selfish = TRUE)
          p_expression <- ggbio::autoplot(as.matrix(log2(assays(se)$CPM[transcript_id,]+1)),axis.text.angle = 45)
          
          
          p <- gridExtra::grid.arrange(p_annotation@ggplot, p_expression, heights = c(1,1))
         
        return(p)
      }else{
        if(!all(gene_id %in% rowData(se)$GENEID)){
          stop("all(gene_id %in% rowData(se)$GENEID) condition is not satisfied!")
        }
        p <- lapply(gene_id, function(g){
          txVec <- rowData(se)[rowData(se)$GENEID == g,]$TXNAME
          txRanges <- rowRanges(se)[txVec]
          names(txRanges) <- paste0(txVec,":", unlist(lapply(strand(txRanges),function(x) unique(as.character(x)))))
          p_annotation <- ggbio::autoplot(txRanges, group.selfish = TRUE)
          p_expression <- ggbio::autoplot(as.matrix(log2(assays(se)$CPM[txVec,]+1)),axis.text.angle = 45, hjust = 1)
          
          
          p <- gridExtra::grid.arrange(p_annotation@ggplot, p_expression, top = g, heights = c(1,1))
          return(p)
        })
        
        return(p)
      }
      
    }
  }
  
  #= 
  count.data <- assays(se)$CPM
  if(length(apply(count.data,1,sum)>10)>100){
    count.data <- count.data[apply(count.data,1,sum)>10,]
  }

  count.data <- log2(count.data+1)




    if(type == "pca"){
      if(!is.null(group.variable)){
        sample.info <- as.data.table(as.data.frame(colData(se)[,c("name",group.variable)]))
        setnames(sample.info, 1:2, c("runname","groupVar"))

        pca_result <- prcomp(t(as.matrix(count.data))) ## can't really cluster them nicely
        plotData <- data.table(pca_result$x[,1:2],keep.rownames = TRUE)
        setnames(plotData, 'rn','runname')
        if(!all(plotData$runname %in% sample.info$runname)){
          stop("all(plotData$runname %in% sample.info$runname) is not satisfied!")
        }
        plotData <- sample.info[plotData, on = 'runname']

        p <- ggplot(plotData, aes(x = PC1, y = PC2))+
          geom_point(aes(col = groupVar))+
          ylab(paste0('PC2 (',round(pca_result$sdev[2]/sum(pca_result$sdev)*100,1),'%)'))+
          xlab(paste0('PC1 (',round(pca_result$sdev[1]/sum(pca_result$sdev)*100,1),'%)'))+
          theme_minimal()

      }else{
        pca_result <- prcomp(t(as.matrix(count.data))) ## can't really cluster them nicely
        plotData <- data.table(pca_result$x[,1:2],keep.rownames = TRUE)
        setnames(plotData, 'rn','runname')

        p <- ggplot(plotData, aes(x = PC1, y = PC2))+
          geom_point(aes(col = runname))+
          ylab(paste0('PC2 (',round(pca_result$sdev[2]/sum(pca_result$sdev)*100,1),'%)'))+
          xlab(paste0('PC1 (',round(pca_result$sdev[1]/sum(pca_result$sdev)*100,1),'%)'))+
          theme_minimal()

      }

      return(p)
    }

  if(type == "heatmap"){

    corData <- cor(count.data,method = "spearman")
    col_fun = circlize::colorRamp2(seq(floor(range(corData)[1]*10)/10,ceiling(range(corData)[2]*10)/10,length.out = 8), RColorBrewer::brewer.pal(8,"Blues"))

    if(!is.null(group.variable)){
      sample.info <- as.data.table(as.data.frame(colData(se)[,c("name",group.variable)]))
      setnames(sample.info, 1:2, c("runname","groupVar"))
      
      if(!all(colnames(count.data) %in% sample.info$runname)){
        stop("all(colnames(count.data) %in% sample.info$runname) is not satisfied!")
      }
      topAnnotation <- ComplexHeatmap::HeatmapAnnotation(group = sample.info[match(colnames(count.data), runname)]$groupVar)

      p <- ComplexHeatmap::Heatmap(corData, name = "Sp.R", col = col_fun,
                   top_annotation = topAnnotation, show_row_names = FALSE,column_names_gp = grid::gpar(fontsize = 9))
    }else{

      p <- ComplexHeatmap::Heatmap(corData, name = "Sp.R", col = col_fun, show_row_names = FALSE,column_names_gp = grid::gpar(fontsize = 9))
      

    }
    return(p)
  }

}
