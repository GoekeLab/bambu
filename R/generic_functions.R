
myPerformance <- function(labels, scores, decreasing = TRUE){
  labels <- labels[order(scores, decreasing = decreasing)]
  results <- list()
  results[['TPR']] <- cumsum(labels)/sum(labels)  # TP/(TP+FP); True Positive Rate;Sensitivity; recall
  results[['FPR']] <- cumsum(!labels)/sum(!labels)  # FP/(FP+TN); False Positive Rate;1-Specificity
  results[['precision']] <- cumsum(labels)/(1:length(labels))  # TP/(TP+FP); positive predictive value;precision
  results[['AUC']] <- sum(results[['TPR']][!duplicated( results[['FPR']],fromLast=TRUE)]/sum(!duplicated( results[['FPR']],fromLast=TRUE)))
  return(results)
}



plotGRangesListByNames<-function(grangesList)
{

  plotTrack = unlist(grangesList)
  plotTrack$grouping <- names(plotTrack)
  plotTrack=AnnotationTrack(plotTrack,group=plotTrack$grouping, id=plotTrack$grouping,fill=2,col.line=1,fontcolor.group=1,fontcolor.item=1)
  plotTracks(list(plotTrack),showId=TRUE)
}
makeTrackFromGrangesList <- function(grangesList)
{

  plotTrack = unlist(grangesList)
  plotTrack$grouping <- names(plotTrack)
  plotTrack=AnnotationTrack(plotTrack,group=plotTrack$grouping, id=plotTrack$grouping,fill=2,col.line=1,fontcolor.group=1,fontcolor.item=1)
  return(plotTrack)
}



grangesListToBed<-function(x) # note: 0 based coordinates
{
  require(GenomicRanges)
  # useful: generate UCSC custom track with the following steps

  # txdbEnsembl91 <- loadDb('/mnt/ont/s3.ontdata.store.genome.sg/annotations/Grch38/ensembl-91/Homo_sapiens.GRCh38.91.annotations-txdb.sqlite')
  # exByTxEns=exonsBy(txdbEnsembl91,by='tx', use.names=TRUE)
  # exByTxSorted <- sort(exByTxEns)
  # genesBed <- grangesListToBed(exByTxSorted)
  # genesBed[,1] <- as.character(genesBed[,1])
  # seqlevelsStyle(genesBed[,1])<-'UCSC'
  # genesBed <- genesBed[order(genesBed$chrom,genesBed$chromStart),]
  # genesBed$score[genesBed$strand=='*'] <- 500
  # genesBed$strand[genesBed$strand=='*'] <- '+'
  #  write.table(genesBed[grepl('chr',genesBed[,1]),],file='/mnt/ont/ensembl-91-ucscChrStyle.bed',col.names=F,row.names=F,sep='\t',quote=F)
  # bedToBigBed  -type=bed12 ensembl-91-ucscChrStyle.bed s3.ontdata.store.genome.sg/annotations/Grch38/ensembl-91/Homo_sapiens.GRCh38.dna_sm.primary_assembly.chrSizes.ucsc.txt ensembl-91-ucscChrStyle.bigbed
  #aws s3 cp ensembl-91-ucscChrStyle.bigbed s3://ucsc-trackdata.store.genome.sg/jonathan-ucsc/Grch38/ --profile ucsc-trackdata.store.genome.sg

  #track type=bigBed name="Ensembl-91" description="Ensembl-91" bigDataUrl=http://s3.ap-southeast-1.amazonaws.com/ucsc-trackdata.store.genome.sg/jonathan-ucsc/Grch38/ensembl-91-ucscChrStyle.bigbed

  #track type=bigBed name="GIS_CRC1685_directRNA_Rep1-Run1-buildTX" description="GIS_CRC1685_directRNA_Rep1-Run1-buildTX" bigDataUrl=http://s3.ap-southeast-1.amazonaws.com/ucsc-trackdata.store.genome.sg/jonathan-ucsc/Grch38/GIS_CRC1685_directRNA_Rep1-Run1_buildTx_ucscChrStyle.bigbed

  #chrom - The name of the chromosome (e.g. chr3, chrY, chr2_random) or scaffold (e.g. scaffold10Gil Sharon1,*, Nikki Jamie Cruz1, Dae-Wook Kang2,3,21, Michael J. Gandal4,5,6,7, Bo Wang1, Young-Mo Kim8, Erika M. Zink8, Cameron P. Casey8, Bryn C. Taylor9, Christianne J. Lane10, Lisa M. Bramer11, Nancy G. Isern8, David W. Hoyt8, Cecilia Noecker12, Michael J. Sweredoski1, Annie Moradian1, Elhanan Borenstein12,13,14,15,16, Janet K. Jansson8, Rob Knight17,18,19, Thomas O. Metz8, Carlos Lois1, Daniel H. Geschwind4,5,6, Rosa Krajmalnik-Brown2,3, and Sarkis K. Mazmanian1,20,*671).
  # chromStart - The starting position of the feature in the chromosome or scaffold. The first base in a chromosome is numbered 0.
  # chromEnd - The ending position of the feature in the chromosome or scaffold. The chromEnd base is not included in the display of the feature. For example, the first 100 bases of a chromosome are defined as chromStart=0, chromEnd=100, and span the bases numbered 0-99.
  # The 9 additional optional BED fields are:
  #
  # name - Defines the name of the BED line. This label is displayed to the left of the BED line in the Genome Browser window when the track is open to full display mode or directly to the left of the item in pack mode.
  # score - A score between 0 and 1000. If the track line useScore attribute is set to 1 for this annotation data set, the score value will determine the level of gray in which this feature is displayed (higher numbers = darker gray). This table shows the Genome Browser's translation of BED score values into shades of gray:
  # strand - Defines the strand. Either "." (=no strand) or "+" or "-".
  # thickStart - The starting position at which the feature is drawn thickly (for example, the start codon in gene displays). When there is no thick part, thickStart and thickEnd are usually set to the chromStart position.
  # thickEnd - The ending position at which the feature is drawn thickly (for example the stop codon in gene displays).
  # itemRgb - An RGB value of the form R,G,B (e.g. 255,0,0). If the track line itemRgb attribute is set to "On", this RBG value will determine the display color of the data contained in this BED line. NOTE: It is recommended that a simple color scheme (eight colors or less) be used with this attribute to avoid overwhelming the color resources of the Genome Browser and your Internet browser.
  # blockCount - The number of blocks (exons) in the BED line.
  # blockSizes - A comma-separated list of the block sizes. The number of items in this list should correspond to blockCount.
  # blockStarts - A comma-separated list of block starts. All of the blockStart positions should be calculated relative to chromStart. The number of items in this list should correspond to blockCount.

  xUnlist <- unlist(range(x))

  bedData = data.frame(as.character(seqnames(xUnlist)),start(xUnlist)-1,end(xUnlist),names(xUnlist),rep(1000,length(xUnlist)),as.character(strand(xUnlist)),start(xUnlist)-1,end(xUnlist),rep(0,length(xUnlist)),elementNROWS(x),unlist(lapply(width(x),paste,collapse=',')),unlist(lapply((start(x)-min(start(x))),paste,collapse=','))) #, stringsAsFactors=FALSE
  colnames(bedData) <-  c('chrom','chromStart','chromEnd','name','score','strand','thickStart','thickEnd','itemRgb','blockCount','blockSizes','blockStarts')
  return(bedData)

  #
  #   xUnlist <- as.data.frame(unlist(range(x)))
  #
  #   bedData = data.frame(as.character(xUnlist[,1]),xUnlist[,2]-1,xUnlist[,3],rownames(xUnlist),rep(1000,nrow(xUnlist)),xUnlist[,5],xUnlist[,2]-1,xUnlist[,3],rep(0,nrow(xUnlist)),elementNROWS(x),unlist(lapply(width(x),paste,collapse=',')),unlist(lapply((start(x)-min(start(x))),paste,collapse=','))) #, stringsAsFactors=FALSE
  #   colnames(bedData) <-  c('chrom','chromStart','chromEnd','name','score','strand','thickStart','thickEnd','itemRgb','blockCount','blockSizes','blockStarts')
  #   return(bedData)
}
