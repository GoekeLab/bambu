#' plotSEOuptut
#' @title plot.bambu
#' @param se An summarized experiment object obtained from \code{\link{bambu}}
#' or \code{\link{transcriptToGeneExpression}}.
#' @param group.variable Variable for grouping in plot, has be to provided if
#' choosing to plot PCA.
#' @param type plot type variable, a values of annotation for a single gene with
#' heatmap for isoform expressions,  pca,  or heatmap, see details.
#' @param gene_id specifying the gene_id for plotting gene annotation, either
#' gene_id or transcript_id has to be provided when type = "annotation".
#' @param transcript_id specifying the transcript_id for plotting transcript
#' annotation, either gene_id or transcript_id has to be provided when 
#' type = "annotation"
#' @details \code{\link{type}} indicates the type of plots to be plotted. There
#' are two types of plots can be chosen, PCA or heatmap.
#' @return A heatmap plot for all samples
#' @export
#' @examples
#' se <- readRDS(system.file("extdata",
#' "seOutputCombined_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.rds",
#' package = "bambu"))
#' plotBambu(se, type = "PCA")
plotBambu <- function(se, group.variable = NULL,
    type = c("annotation", "pca", "heatmap"),
    gene_id = NULL, transcript_id = NULL) {
    if (type == "annotation") {
        p <- plotAnnotation(se, gene_id, transcript_id)
        return(p)
    }
    # =
    count.data <- assays(se)$CPM
    count.data <- count.data[apply(count.data, 1, sd) >
        quantile(apply(count.data, 1, sd), 0.50), ]
    count.data <- log2(count.data + 1)
    if (type == "pca") {
        p <- plotPCA(se, count.data, group.variable)
        return(p)
    }
    if (type == "heatmap") {
        p <- plotHeatmap(se, count.data, group.variable)
        return(p)
    }
}
