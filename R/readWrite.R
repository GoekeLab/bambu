#' Outputs a GTF file, transcript-count file, and gene-count file from bambu
#' @title Write bambu results to GTF and transcript/gene-count files
#' @param se a \code{\link{SummarizedExperiment}} object from \code{\link{bambu}}
#' @param path the destination of the output files (gtf, transcript counts, and gene counts)
#' @return The function will generate three files, a \code{\link{.gtf}} file for the annotations, 
#' two \code{\link{.txt}} files for transcript and gene counts respectively. 
#' @export
writeBambuOutput <- function(se,path){
  if (missing(se) | missing(path)){
    stop('Both summarizedExperiment object from bambu and the path for the output files are required.')
  }else{
    outdir <- paste0(path,"/")
    if (!dir.exists(outdir)){
      dir.create(outdir, recursive = TRUE)
    }
    transcript_grList <- rowRanges(se)
    transcript_gtffn <- paste(path,"transcript_exon.gtf",sep="")
    gtf <- writeToGTF(annotation=transcript_grList,file=transcript_gtffn)
    transcript_counts <- as.data.frame(assays(se)$counts)
    geneIDs <- data.frame(mcols(transcript_grList, use.names=FALSE)$TXNAME,mcols(transcript_grList, use.names=FALSE)$GENEID)
    colnames(geneIDs) <- c("TXNAME","GENEID")
    transcript_counts <- cbind(geneIDs,transcript_counts)  
    transcript_countsfn <- paste(path,"counts_transcript.txt",sep="")
    write.table(transcript_counts, file= transcript_countsfn, sep="\t",quote=FALSE,row.names= FALSE)
    gene_se <- transcriptToGeneExpression(se)
    gene_counts <- as.data.frame(assays(gene_se)$counts)
    gene_countsfn <- paste(path,"counts_gene.txt",sep="")
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
#' @param file a \code{\link{.gtf}} file
#' @return grlist a \code{\link{GRangesList}} object, with two columns
#' \itemize{
#'   \item TXNAME specifying prefix for new gene Ids (genePrefix.number), defaults to empty
#'   \item GENEID indicating whether filter to remove read classes which are a subset of known transcripts(), defaults to TRUE
#'   }
#' @export
readFromGTF <- function(file){
  if (missing(file)){
    stop('A GTF file is required.')
  }else{
    data=read.delim(file,header=FALSE, comment.char='#')

    colnames(data) <- c("seqname","source","type","start","end","score","strand","frame","attribute")
    data <- data[data$type=='exon',]
        data$strand[data$strand=='.'] <- '*'
    data$GENEID = gsub('gene_id (.*?);.*','\\1',data$attribute)
    data$TXNAME=gsub('.*transcript_id (.*?);.*', '\\1',data$attribute)
    grlist <- makeGRangesListFromDataFrame(
      data[,c('seqname', 'start','end','strand','TXNAME')],split.field='TXNAME',keep.extra.columns = TRUE)
    grlist <- grlist[IRanges::order(start(grlist))]
    
    geneData=(unique(data[,c('TXNAME', 'GENEID')]))
    mcols(grlist) <- DataFrame(geneData[(match(names(grlist), geneData$TXNAME)),])
    
  }
  return (grlist)
}

