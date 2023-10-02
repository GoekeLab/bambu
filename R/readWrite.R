#' @title Write Bambu results to GTF and transcript/gene-count files
#' @param se a \code{\link{SummarizedExperiment}} object 
#' from \code{\link{bambu}}.
#' @param path the destination of the output files 
#' (gtf, transcript counts, and gene counts)
#' @param prefix the prefix of the output files
#' @details The function will write the output from Bambu to files. The 
#' annotations will be written to a .gtf file, transcript counts (total counts, 
#' CPM, full-length counts, partial-length counts, and unique counts) and gene counts
#' will be written to .txt files.  
#' @export
#' @examples
#' se <- readRDS(system.file("extdata",
#'     "seOutput_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.rds",
#'     package = "bambu"
#' ))
#' path <- tempdir()
#' writeBambuOutput(se, path)
writeBambuOutput <- function(se, path, prefix = "") {
    if (missing(se) | missing(path)) {
        stop("Both summarizedExperiment object from bambu and
            the path for the output files are required.")
    } else {
        outdir <- paste0(path, "/")
        if (!dir.exists(outdir))
            dir.create(outdir, recursive = TRUE)

        transcript_grList <- rowRanges(se)
        transcript_gtffn <- paste(outdir, prefix,
            "extended_annotations.gtf", sep = "")
        gtf <- writeToGTF(annotation = transcript_grList,
            file = transcript_gtffn)
        
        for(d in names(assays(se))){
            writeCountsOutput(se, varname=d,
                             feature='transcript',outdir, prefix)
        }
        
        seGene <- transcriptToGeneExpression(se[,1])
        
        for (i in seq(ceiling(length(colnames(se)) / 100))){
          if (i == ceiling(length(colnames(se)) / 100)){
              seGene <- cbind(seGene, transcriptToGeneExpression(se[,(100*i-99):length(colnames(se))]))
          } else{
              seGene <- cbind(seGene, transcriptToGeneExpression(se[,(100*i-99):(100*i)]))
          }
        }
        
        seGene <- seGene[,-1]
        
        writeCountsOutput(seGene, varname='counts', feature='gene',outdir, prefix)
    }
}


#' helper function to write counts 
#' @noRd
writeCountsOutput <- function(se, varname = "counts",
                              feature = "transcript", outdir, prefix){
  
    estimatesfn <- paste(outdir, prefix, varname,"_",feature,".txt", sep = "")
    if(!is(assays(se)[[varname]], "sparseMatrix")){
      estimates <- data.table(as.data.frame(assays(se)[[varname]]),
                              keep.rownames = TRUE) 
      if(feature == "transcript"){
        setnames(estimates, "rn", "TXNAME")
        geneIDs <- data.table(as.data.frame(rowData(se))[,c("TXNAME","GENEID")])
        estimates <- geneIDs[estimates, on = "TXNAME"]
      }else{
        setnames(estimates, "rn","GENEID")
      } 
      
      utils::write.table(estimates, file = estimatesfn, sep = "\t", quote = FALSE, row.names = FALSE)
      
    } else{
        estimates <- assays(se)[[varname]]
        if (feature == "transcript"){
          Matrix::writeMM(estimates, estimatesfn)
          geneIDs <- data.table(as.data.frame(rowData(se))[,c("TXNAME","GENEID")])
          utils::write.table(geneIDs, file = paste0(outdir, "txANDgenes.txt"), 
                             sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
          utils::write.table(colnames(se), file = paste0(outdir, "barcodes.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
          
        } else{
          Matrix::writeMM(estimates, estimatesfn)
          utils::write.table(names(se), file = paste0(outdir, "genes.txt"), 
                             sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
        }
    }
}

#' Write annotation GRangesList into a GTF file
#' @title write GRangeslist into GTF file
#' @param annotation a \code{GRangesList} object
#' @param file the output gtf file name
#' @param geneIDs an optional dataframe of geneIDs (column 2) with
#'     the corresponding transcriptIDs (column 1)
#' @return gtf a GTF dataframe
#' @importFrom dplyr select as_tibble mutate %>% left_join arrange group_by
#'     ungroup recode_factor
#' @importFrom methods is
#' @export
#' @examples
#' outputGtfFile <- tempfile()
#' gr <- readRDS(system.file("extdata",
#'     "annotationGranges_txdbGrch38_91_chr9_1_1000000.rds",
#'     package = "bambu"
#' ))
#' writeToGTF(gr, outputGtfFile)
writeToGTF <- function(annotation, file, geneIDs = NULL) {
    if (missing(annotation) | missing(file)) {
        stop("Both GRangesList and the name of the output file are required.")
    } else if (!is(annotation, "CompressedGRangesList")) {
        stop("The inputted GRangesList is of the wrong class.")
    }
    df <- as_tibble(annotation)
    df$exon_rank <- paste('exon_number "', df$exon_rank, '";', sep = "")
    if (missing(geneIDs)) {
        if (!is.null(mcols(annotation, use.names = FALSE)$GENEID)) {
            geneIDs <- as_tibble(mcols(annotation, use.names = FALSE)[,
                c("TXNAME", "GENEID")])
            geneIDs$seqnames <- unlist(unique(seqnames(annotation)))
            geneIDs$strand <- unlist(unique(strand(annotation)))
        }
    }
    df <- left_join(df, geneIDs[, c("TXNAME", "GENEID")], 
        by = c("group_name" = "TXNAME"))
    df$group_name <- paste('transcript_id "', df$group_name, '";', sep = "")
    df$GENEID <- paste('gene_id "', df$GENEID, '";', sep = "")
    dfExon <- mutate(df, source = "Bambu", feature = "exon", score = ".",
        frame = ".", attributes = paste(GENEID, group_name, exon_rank)) %>%
        select(seqnames, source, feature, start, end, score,
        strand, frame, attributes, group_name)
    dfTx <- as.data.frame(range(ranges(annotation)))
    dfTx <-
        left_join(dfTx, geneIDs, by = c("group_name" = "TXNAME"))
    dfTx$group_name <-
        paste('transcript_id "', dfTx$group_name, '";', sep = "")
    dfTx$GENEID <- paste('gene_id "', dfTx$GENEID, '";', sep = "")

    dfTx <- mutate(dfTx,source = "Bambu", feature = "transcript", score = ".",
        frame = ".", attributes = paste(GENEID, group_name)) %>%
        select(seqnames, source, feature, start, end, score,
        strand, frame, attributes, group_name)

    gtf <- rbind(dfTx, dfExon) %>% group_by(group_name) %>%
        arrange(as.character(seqnames), start) %>% ungroup() %>%
        select(seqnames, source, feature, start, end, score,
        strand, frame, attributes)
    gtf <- mutate(gtf, strand = recode_factor(strand, `*` = "."))
    utils::write.table(gtf, file = file, quote = FALSE, row.names = FALSE,
        col.names = FALSE, sep = "\t")
}


#' Outputs GRangesList object from reading a GTF file
#' @title convert a GTF file into a GRangesList
#' @param file a .gtf file
#' @param keep.extra.columns a vector with names of columns to keep from 
#' the the attributes in the gtf file. For ensembl, this could be 
#' keep.extra.columns=c('gene_name','gene_biotype',
#' 'transcript_biotype', 'transcript_name')
#' @return grlist a \code{GRangesList} object, with two columns
#' \describe{
#'     \item{TXNAME}{specifying prefix for new gene Ids (genePrefix.number),
#'         defaults to empty}
#'     \item{GENEID}{indicating whether filter to remove read classes which are
#'         a subset of known transcripts(), defaults to TRUE}
#'   }
#' @importFrom GenomicRanges makeGRangesListFromDataFrame
#' @export
#' @examples
#' gtf.file <- system.file("extdata",
#'     "Homo_sapiens.GRCh38.91_chr9_1_1000000.gtf",
#'     package = "bambu"
#' )
#' readFromGTF(gtf.file)
readFromGTF <- function(file, keep.extra.columns = NULL){
    if (missing(file)) {
        stop('A GTF file is required.')
    }else{
        data <- utils::read.delim(file,header = FALSE, comment.char = '#')
        colnames(data) <- c("seqname","source","type","start",
            "end","score","strand","frame","attribute")
    data <- data[data$type == 'exon',]
    data$strand[data$strand == '.'] <- '*'
    data$GENEID = gsub('gene_id (.*?);.*','\\1',data$attribute)
    data$TXNAME = gsub('.*transcript_id (.*?);.*', '\\1',data$attribute)
    data$exon_rank = gsub('.*exon_number (.*?);.*', '\\1',data$attribute)
    if (!is.null(keep.extra.columns)) {
        for (extraColumn in seq_along(keep.extra.columns)) {
            data[,keep.extra.columns[extraColumn]] <-
                gsub(paste0('.*',keep.extra.columns[extraColumn],
                ' (.*?);.*'), '\\1',data$attribute)
            data[grepl(';', data[,keep.extra.columns[extraColumn]]),
                keep.extra.columns[extraColumn]] <- ''
        }
    }
    grlist <- makeGRangesListFromDataFrame(
        data[,c('seqname', 'start','end','strand','TXNAME','exon_rank')],
        split.field = 'TXNAME',keep.extra.columns = TRUE)
    grlist <- grlist[IRanges::order(start(grlist))]
    geneData = (unique(data[,c('TXNAME', 'GENEID',keep.extra.columns)]))
    mcols(grlist) <-
        DataFrame(geneData[(match(names(grlist), geneData$TXNAME)),])
    }
    return(grlist)
}
