#' plotSEOuptut
#' @title plot.bambu
#' @param se An summarized experiment object obtained from \code\link{bambu} or \code{\link{transcriptToGene}}.
#' @param group.variable Variable for grouping in plot, has be to provided if choosing to plot PCA.
#' @param type plot type variable, a values of pca,  or heatmap, see \code{\link{details}}.
#' @details \code{\link{type}} indicates the type of plots to be plotted. There are two types of plots can be chosen, PCA or heatmap.
#' @return A heatmap plot for all samples
#' @importFrom stats prcomp
#' @importFrom ggplot2 ggplot
#' @importFrom ComplexHeatmap Heatmap
#' @importFrom circlize colorRamp2
#' @importFrom RColorBrewer brewer.pal
#' @export
plot.bambu <- function(se, group.variable = NULL, type ){

  count.data <- assays(se)$CPM
  if(length(apply(count.data,1,sum)>10)>100){
    count.data <- count.data[apply(count.data,1,sum)>10,]
  }

  count.data <- log2(count.data+1)




    if(type == "pca"){
      if(!is.null(group.variable)){
        sample.info <- as.data.table(as.data.frame(colData(se)[,c("name",group.variable)]))
        setnames(sample.info, 1:2, c("runname","groupVar"))

        pca_result <- prcomp(t(as.matrix(count.data)), scale = TRUE) ## can't really cluster them nicely
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
        pca_result <- prcomp(t(as.matrix(count.data)), scale = TRUE) ## can't really cluster them nicely
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
    col_fun = colorRamp2(seq(floor(range(corData)[1]*10)/10,ceiling(range(corData)[2]*10)/10,length.out = 8), brewer.pal(8,"Blues"))

    if(!is.null(group.variable)){
      if(!all(colnames(count.data) %in% sample.info$runname)){
        stop("all(colnames(count.data) %in% sample.info$runname) is not satisfied!")
      }
      topAnnotation <- HeatmapAnnotation(group = sample.info[match(colnames(count.data), runname)]$groupVar)

      p <- Heatmap(corData, name = "Sp.R", col = col_fun,
                   top_annotation = topAnnotation)

    }else{

      p <- Heatmap(corData, name = "Sp.R", col = col_fun)

    }
    return(p)
  }

}
