#' Outputs a GTF dataframe for the nanorna-bam nextflow pipeline 
#' @title transcript to gene expression
#' @param se a summarizedExperiment object from \code{\link{bambu}}
#' @export
SEtoGTF <- function(se){
  if (missing(se)){
    stop("The summarizedExperiment object from bambu is required.")
  }else{
    exp_dat <-cbind(as.data.frame(rowData(se))[,c(1,2,4)],data.frame(assays(se)$counts,assays(se)$CPM))
    df <- merge(as.data.frame(rowRanges(se)),exp_dat,by.x="group_name",by.y="TXNAME",all=TRUE)
    df$GENEID <- paste('gene_id "',df$GENEID,'";',sep= '')
    df$group_name <- paste('transcript_id "',df$group_name,'";',sep= '')
    df$exon_rank <- paste('exon_number "',df$exon_rank,'";',sep= '')
    gtf_exon <- data.frame(seqname=df$seqnames, source= "Bambu",feature= "exon",
                           start=df$start,end=df$end,score=".",strand=df$strand,frame=".",
                           attributes= paste(df$GENEID,df$group_name,df$exon_rank))
    df <- df[!duplicated(df$group_name),]
    df$GIS_MCF7_directcDNA_Replicate1_genome.1 <- paste('CPM "',df$GIS_MCF7_directcDNA_Replicate1_genome.1,'";',sep= '')
    gtf_trns <- data.frame(seqname=df$seqnames, source= "Bambu",feature= "transcript",
                           start=df$start,end=df$end,score=".",strand=df$strand,frame=".",
                           attributes= paste(df$GENEID,df$group_name,df$GIS_MCF7_directcDNA_Replicate1_genome.1))
    gtf <- rbind(gtf_trns,gtf_exon)
  }
  return(gtf)
}

