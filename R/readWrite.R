#' Outputs a GTF file, transcript-count file, and gene-count file from bambu
#' @title write bambu results to GTF and transcript/gene-count files
#' @param se a summarizedExperiment object from \code{\link{bambu}}
#' @param path the destination of the output files (gtf, transcript counts, and gene counts)
#' @return gtf a GTF dataframe
#' @export
write.bambu <- function(se,path){
  if (missing(se) | missing(path)){
    stop('Both summarizedExperiment object from bambu and the path for the output files are required.')
  }else{
    outdir <- strsplit(path,"/")[[1]][1]
    if (dir.exists(outdir) == FALSE){
      dir.create(outdir)
    }
    transcript_grList <- rowRanges(se)
    transcript_gtffn <- paste(path,"transcript_exon.gtf",sep="")
    gtf <- write.gtf(annotation=transcript_grList,file=transcript_gtffn)
    transcript_counts <- as.data.frame(assays(se)$counts)
    geneIDs <- data.frame(mcols(transcript_grList, use.names=FALSE)$TXNAME,mcols(transcript_grList, use.names=FALSE)$GENEID)
    colnames(geneIDs) <- c("TXNAME","GENEID")
    transcript_counts <- cbind(geneIDs,transcript_counts)  
    transcript_countsfn <- paste(path,"counts_transcript.txt",sep="")
    write.table(transcript_counts, file= transcript_countsfn, sep="\t",quote=FALSE)
    gene_se <- transcriptToGeneExpression(se)
    gene_counts <- as.data.frame(assays(gene_se)$counts)
    gene_countsfn <- paste(path,"counts_gene.txt",sep="")
    write.table(gene_counts, file= gene_countsfn, sep="\t", quote=FALSE)
  }
}
#' Outputs a GTF file for the nanorna-bam nextflow pipeline 
#' @title write GRangeslist into GTF file
#' @param annotation a GRangesList object
#' @param file the output gtf file name
#' @param geneIDs an optional dataframe of geneIDs (column 2) with the corresponding transcriptIDs (column 1)
#' @return gtf a GTF dataframe
#' @export
write.gtf <- function (annotation,file,geneIDs=NULL) {
  if (missing(annotation) | missing(file)){
    stop('Both GRangesList and the name of file are required.')
  }else if (class(annotation) != "CompressedGRangesList"){
    stop('The inputted GRangesList is of the wrong class.')
  }else {
    df <- as.data.frame(annotation)
    df$exon_rank <- paste('exon_number "',df$exon_rank,'";',sep= '')
    if (missing(geneIDs)){
      if (!is.null(mcols(annotation, use.names=FALSE)$GENEID)){
        geneIDs <- data.frame(mcols(annotation, use.names=FALSE)$TXNAME,mcols(annotation, use.names=FALSE)$GENEID)
        colnames(geneIDs) <- c("TXNAME","GENEID")
      }
    }
    df <- merge(df,geneIDs,by.x="group_name",by.y="TXNAME",all=TRUE)
    df$group_name <- paste('transcript_id "',df$group_name,'";',sep= '')
    df$GENEID <- paste('gene_id "',df$GENEID,'";',sep= '')
    gtf_exon <- data.frame(seqname=df$seqnames, source= "Bambu",feature= "exon",
                          start=df$start,end=df$end,score=".",strand=df$strand,frame=".",
                          attributes= paste(df$GENEID,df$group_name,df$exon_rank))
    df_end <- df[order(df$end, decreasing = TRUE),]
    df_end <- data.frame(group_name=df_end$group_name,uend=df_end$end)
    df_end <- df_end[!duplicated(df_end$group_name),]
    df <- df[!duplicated(df$group_name),]
    df <- merge(df,df_end,by="group_name",all=TRUE)
    gtf_trns <- data.frame(seqname=df$seqnames, source= "Bambu",feature= "transcript",
                             start=df$start,end=df$uend,score=".",strand=df$strand,frame=".",
                             attributes= paste(df$GENEID,df$group_name))
    }
    gtf <- rbind(gtf_trns,gtf_exon)
    gtf <- gtf[order(gtf$attributes),]
    write.table(gtf, file= file, quote=FALSE, row.names=FALSE, col.names=FALSE, sep = "\t")
  } 
#' Outputs GRangesList object from reading a GTF file
#' @title convert a GTF file into a GRangesList
#' @param file a GTF file
#' @return grlist a GRangesList object
#' @export
read.gtf <- function(file){
  if (missing(file)){
    stop('A GTF file is required.')
  }else{
    data=read.delim(file,header=FALSE)
    colnames(data) <- c("seqname","source","type","start","end","score","strand","frame","attribute")
    data$strand[data$strand=='.'] <- '*'
    data$GENEID = gsub('gene_id (.*); tra.*','\\1',data$attribute)
    data$TXNAME=gsub('.*transcript_id (.*?);.*', '\\1',data$attribute)
    #geneData=unique(data[,c('txid','gene_id','refMatch')]))
    geneData=unique(data[,c('TXNAME', 'GENEID')])
    grlist <- makeGRangesListFromDataFrame(
      data[data$type=='exon',],split.field='TXNAME', names.field='TXNAME',keep.extra.columns = TRUE)
    geneData=(unique(data[,c('TXNAME', 'GENEID')]))
    mcols(grlist) <- DataFrame(geneData[(match(names(grlist), geneData$txid)),])
  }
  return (grlist)
}

#' Outputs GRangesList object from reading a GTF file for running bambu
#' @title convert a GTF file into a GRangesList
#' @param file a GTF file
#' @return grlist a GRangesList object
#' @export
prepareAnnotationsFromGTF_draft <- function(file){
  if (missing(file)){
    stop('A GTF file is required.')
  }else{
    data <- read.delim(file,header=FALSE,comment.char='#')
    colnames(data) <- c("seqname","source","type","start","end","score","strand","frame","attribute")
    data <- data[data$type=='exon',]
    data$strand[data$strand=='.'] <- '*'
    data$GENEID = gsub('gene_id (.*?);.*','\\1',data$attribute)
    data$TXNAME=gsub('.*transcript_id (.*?);.*', '\\1',data$attribute)
    data$exon_rank=as.integer(gsub('.*exon_number (.*?);.*', '\\1',data$attribute))
    geneData=unique(data[,c('TXNAME', 'GENEID')])
    grlist <- makeGRangesListFromDataFrame(
      data[,c('seqname', 'start','end','strand','exon_rank','TXNAME')],split.field='TXNAME',keep.extra.columns = TRUE)
    
    unlistedExons <- unlist(grlist, use.names = FALSE)
    partitioning <- PartitioningByEnd(cumsum(elementNROWS(grlist)), names=NULL)
    txIdForReorder <- togroup(PartitioningByWidth(grlist))
    unlistedExons <- unlistedExons[order(txIdForReorder, unlistedExons$exon_rank)] #'exonsByTx' is always sorted by exon rank, not by strand, make sure that this is the case here
    unlistedExons$exon_endRank <- unlist(sapply(elementNROWS(grlist),seq,to=1), use.names=FALSE)
    unlistedExons <- unlistedExons[order(txIdForReorder, start(unlistedExons))]
#    mcols(unlistedExons) <- mcols(unlistedExons)[,c('exon_rank','exon_endRank')]
    grlist <- relist(unlistedExons, partitioning)
    minEqClasses <- getMinimumEqClassByTx(grlist)
    mcols(grlist) <- DataFrame(geneData[(match(names(grlist), geneData$TXNAME)),])
    mcols(grlist)$eqClass <- minEqClasses$eqClass[match(names(grlist),minEqClasses$queryTxId)]
  }
  return (grlist)
}
