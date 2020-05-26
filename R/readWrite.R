#' Outputs a GTF file, transcript-count file, and gene-count file from bambu
#' @title transcript to gene expression
#' @param se a summarizedExperiment object from \code{\link{bambu}}
#' @param path the destination of the output files (gtf, transcript counts, and gene counts)
#' @export
write.bambu <- function(transcript_se,path){
  if (missing(transcript_se) | missing(path)){
    stop('Both summarizedExperiment object from bambu and the path for the output files are required.')
  }else{
    dir.create(path)
    transcript_grList <- rowRanges(transcript_se)
    transcript_gtffn <- paste(path,"bambu_transcript_exon.gtf",sep="/")
    transcript_geneIDs <-as.data.frame(rowData(transcript_se))[,c(1,2)]
    write.gtf(grList=transcript_grList,file=transcript_gtffn,geneIDs=transcript_geneIDs)
    transcript_counts <- as.data.frame(assays(transcript_se)$counts)
    transcript_countsfn <- paste(path,"counts_transcript.txt",sep="/")
    write.table(transcript_counts, file= transcript_countsfn, sep="\t",quote=FALSE)
    gene_se <- transcriptToGeneExpression(se)
    gene_counts <- as.data.frame(assays(gene_se)$counts)
    gene_countsfn <- paste(path,"counts_gene.txt",sep="/")
    write.table(gene_counts, file= gene_countsfn, sep="\t", quote=FALSE)
    }
}
#' Outputs a GTF file for the nanorna-bam nextflow pipeline 
#' @title transcript to gene expression
#' @param grList a GRangesList object
#' @param file the output gtf file name
#' @param geneIDs an optional dataframe of geneIDs (column 2) with the corresponding transcriptIDs (column 1)
#' @export
write.gtf <- function (grList,file,geneIDs) {
  if (missing(grList) | missing(file)){
    stop('Both GRangesList and the name of file are required.')
  }else if (class(grList) != "CompressedGRangesList"){
    stop('The inputted GRangesList is of the wrong class.')
  }else {
    df <- as.data.frame(grList)
    df$exon_rank <- paste('exon_number "',df$exon_rank,'";',sep= '')
    if (missing(geneIDs)){
      df$group_name <- paste('transcript_id "',df$group_name,'";',sep= '')
      gtf_exon <- data.frame(seqname=df$seqnames, source= "Bambu",feature= "exon",
                             start=df$start,end=df$end,score=".",strand=df$strand,frame=".",
                             attributes= paste(df$group_name,df$exon_rank))
    }else{
      df <- merge(df,geneIDs,by.x="group_name",by.y="TXNAME",all=TRUE)
      df$group_name <- paste('transcript_id "',df$group_name,'";',sep= '')
      df$GENEID <- paste('gene_id "',df$GENEID,'";',sep= '')
      gtf_exon <- data.frame(seqname=df$seqnames, source= "Bambu",feature= "exon",
                             start=df$start,end=df$end,score=".",strand=df$strand,frame=".",
                             attributes= paste(df$GENEID,df$group_name,df$exon_rank))
    }
    df_end <- df[order(df$end, decreasing = TRUE),]
    df_end <- data.frame(group_name=df_end$group_name,uend=df_end$end)
    df_end <- df_end[!duplicated(df_end$group_name),]
    df <- df[!duplicated(df$group_name),]
    df <- merge(df,df_end,by="group_name",all=TRUE)
    if (missing(geneIDs)){
      gtf_trns <- data.frame(seqname=df$seqnames, source= "Bambu",feature= "transcript",
                             start=df$start,end=df$uend,score=".",strand=df$strand,frame=".",
                             attributes= paste(df$group_name))
    }else{
      gtf_trns <- data.frame(seqname=df$seqnames, source= "Bambu",feature= "transcript",
                             start=df$start,end=df$uend,score=".",strand=df$strand,frame=".",
                             attributes= paste(df$GENEID,df$group_name))
    }
    gtf <- rbind(gtf_trns,gtf_exon)
    gtf <- gtf[order(gtf$attributes),]
  } 
  write.table(gtf, file= file, quote=FALSE, row.names=FALSE, col.names=FALSE, sep = "\t")
}
#' Outputs GRangesList object from reading a GTF file
#' @title transcript to gene expression
#' @param gtf_file a GTF file
#' @return grlist a GRangesList object
#' @export
read.gtf <- function(gtf_file){
  if (missing(gtf_file)){
    stop('A GTF file is required.')
  }else{
    data=read.delim(gtf_file,header=FALSE)
    colnames(data) <- c("seqname","source","type","start","end","score","strand","frame","attribute")
    data$strand[data$strand=='.'] <- '*'
    data$gene_id = gsub('gene_id (.*); tra.*','\\1',data$attribute)
    data$txid=gsub('.*transcript_id (.*?);.*', '\\1',data$attribute)
    #geneData=unique(data[,c('txid','gene_id','refMatch')]))
    geneData=unique(data[,c('txid','gene_id')])
    grlist <- makeGRangesListFromDataFrame(
      data[data$type=='exon',],split.field='txid', names.field='txid',keep.extra.columns = TRUE)
    geneData=(unique(data[,c('txid','gene_id')]))
    mcols(grlist) <- DataFrame(geneData[(match(names(grlist), geneData$txid)),])
  }
  return (grlist)
}