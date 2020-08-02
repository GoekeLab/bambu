#' Outputs a GTF file, transcript-count file, and gene-count file from bambu
#' @title Write bambu results to GTF and transcript/gene-count files
#' @param se a \code{\link{SummarizedExperiment}} object from \code{\link{bambu}}
#' @param path the destination of the output files (gtf, transcript counts, and gene counts)
#' @return The function will generate three files, a .gtf file for the annotations, 
#' two .txt files for transcript and gene counts respectively. 
#' @export
#' @examples 
#' se <- readRDS(system.file("extdata", 
#' "seOutput_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.rds", 
#' package = "bambu"))
#' path <- tempdir()
#' writeBambuOutput(se,path)
writeBambuOutput <- function(se, path, prefix=''){
  if (missing(se) | missing(path)){
    stop('Both summarizedExperiment object from bambu and the path for the output files are required.')
  }else{
    outdir <- paste0(path,"/")
    if (!dir.exists(outdir)){
      dir.create(outdir, recursive = TRUE)
    }
    transcript_grList <- rowRanges(se)
    transcript_gtffn <- paste(outdir, prefix, "extended_annotations.gtf",sep="")
    gtf <- writeToGTF(annotation=transcript_grList,file=transcript_gtffn)
    transcript_counts <- as.data.frame(assays(se)$counts)
    geneIDs <- data.frame(mcols(transcript_grList, use.names=FALSE)$TXNAME,mcols(transcript_grList, use.names=FALSE)$GENEID)
    colnames(geneIDs) <- c("TXNAME","GENEID")
    transcript_counts <- cbind(geneIDs,transcript_counts)  
    transcript_countsfn <- paste(outdir,prefix, "counts_transcript.txt",sep="")
    write.table(transcript_counts, file= transcript_countsfn, sep="\t",quote=FALSE,row.names= FALSE)
    gene_se <- transcriptToGeneExpression(se)
    gene_counts <- as.data.frame(assays(gene_se)$counts)
    gene_countsfn <- paste(outdir, prefix, "counts_gene.txt",sep="")
    write.table(gene_counts, file= gene_countsfn, sep="\t", quote=FALSE)
  }
}
#' Write annotation GRangesList into a GTF file
#' @title write GRangeslist into GTF file
#' @param annotation a \code{\link{GRangesList}} object
#' @param file the output gtf file name
#' @param geneIDs an optional dataframe of geneIDs (column 2) with the corresponding transcriptIDs (column 1)
#' @return gtf a GTF dataframe
#' @export
#' @examples
#' outputGtfFile <- tempfile()
#' gr <- readRDS(system.file("extdata", 
#' "annotationGranges_txdbGrch38_91_chr9_1_1000000.rds", 
#' package = "bambu"))
#' writeToGTF(gr, outputGtfFile)
writeToGTF <- function (annotation,file,geneIDs=NULL) {
  if (missing(annotation) | missing(file)){
    stop('Both GRangesList and the name of the output file are required.')
  }else if (!is(annotation,"CompressedGRangesList")){
    stop('The inputted GRangesList is of the wrong class.')
  }
  df <- as_tibble(annotation)
  df$exon_rank <- paste('exon_number "',df$exon_rank,'";',sep= '')
  if (missing(geneIDs)){
    if (!is.null(mcols(annotation, use.names=FALSE)$GENEID)){
      geneIDs <- as_tibble(mcols(annotation, use.names=FALSE)[, c('TXNAME','GENEID')])
      geneIDs$seqnames = unlist(unique(seqnames(annotation)))
      geneIDs$strand = unlist(unique(strand(annotation)))
    }
  }
  df <- left_join(df,geneIDs[,c('TXNAME','GENEID')],by=c("group_name"="TXNAME"),all=TRUE)
  df$group_name <- paste('transcript_id "',df$group_name,'";',sep= '')
  df$GENEID <- paste('gene_id "',df$GENEID,'";',sep= '')
  dfExon <- mutate(df, source='Bambu',
                   feature='exon',
                   score='.',
                   frame='.',
                   attributes=paste(GENEID, group_name, exon_rank)) %>%
    dplyr::select( seqnames, source, feature, start, end, score, strand, frame, attributes, group_name)
  dfTx <- as.data.frame(range(ranges(annotation)))
  dfTx <- left_join(dfTx,geneIDs,by=c("group_name"="TXNAME"),all=TRUE)
  dfTx$group_name <- paste('transcript_id "',dfTx$group_name,'";',sep= '')
  dfTx$GENEID <- paste('gene_id "',dfTx$GENEID,'";',sep= '')
  
  dfTx <- mutate(dfTx, source='Bambu',
                 feature='transcript',
                 score='.',
                 frame='.',
                 attributes=paste(GENEID, group_name)) %>%
    dplyr::select( seqnames, source, feature, start, end, score, strand,frame, attributes, group_name)
  gtf <- rbind(dfTx, dfExon) %>%
    group_by(group_name) %>% 
    arrange(as.character(seqnames), start) %>%
    ungroup() %>%
    dplyr::select(seqnames, source, feature, start, end, score, strand, frame, attributes)
  gtf <- mutate(gtf, strand=recode_factor(strand, `*`="."))
  write.table(gtf, file= file, quote=FALSE, row.names=FALSE, col.names=FALSE, sep = "\t")
} 


#' Outputs GRangesList object from reading a GTF file
#' @title convert a GTF file into a GRangesList
#' @param file a .gtf file
#' @param keep.extra.columns a vector with names of columns to keep from the the attributes in the gtf file. For ensembl, this could be keep.extra.columns=c('gene_name','gene_biotype','transcript_biotype', 'transcript_name')
#' @return grlist a \code{\link{GRangesList}} object, with two columns
#' \itemize{
#'   \item TXNAME specifying prefix for new gene Ids (genePrefix.number), defaults to empty
#'   \item GENEID indicating whether filter to remove read classes which are a subset of known transcripts(), defaults to TRUE
#'   }
#' @export
#' @examples
#' gtf.file <- system.file("extdata", 
#' "Homo_sapiens.GRCh38.91_chr9_1_1000000.gtf", 
#' package = "bambu")
#' readFromGTF(gtf.file)
readFromGTF <- function(file, keep.extra.columns = NULL){
  if (missing(file)){
    stop('A GTF file is required.')
  }else{
    data=read.delim(file,header=FALSE, comment.char='#')

    colnames(data) <- c("seqname","source","type","start","end","score","strand","frame","attribute")
    data <- data[data$type=='exon',]
        data$strand[data$strand=='.'] <- '*'
    data$GENEID = gsub('gene_id (.*?);.*','\\1',data$attribute)
    data$TXNAME=gsub('.*transcript_id (.*?);.*', '\\1',data$attribute)
    if(!is.null(keep.extra.columns)){
      for(extraColumn in seq_along(keep.extra.columns)){
      data[,keep.extra.columns[extraColumn]] <- gsub(paste0('.*',keep.extra.columns[extraColumn],' (.*?);.*'), '\\1',data$attribute)
      data[grepl(';', data[,keep.extra.columns[extraColumn]]),keep.extra.columns[extraColumn]] <- ''
      }
      }
    grlist <- makeGRangesListFromDataFrame(
      data[,c('seqname', 'start','end','strand','TXNAME')],split.field='TXNAME',keep.extra.columns = TRUE)
    grlist <- grlist[IRanges::order(start(grlist))]
    
    geneData=(unique(data[,c('TXNAME', 'GENEID',keep.extra.columns)]))
    mcols(grlist) <- DataFrame(geneData[(match(names(grlist), geneData$TXNAME)),])
    
  }
  return (grlist)
}

