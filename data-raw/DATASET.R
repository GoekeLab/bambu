## code to prepare `sysdata.rda` dataset goes here


data1 <- data.frame(
    tx_id = c( 1,"1Start", 2, "2Start"),
    read_class_id = c(1,"1Start", 2,"2Start"),
    nobs = c(10, 50, 50,10),
    tx_len = c(546,546,2356,2356),
    rc_width = c(300,540,1800,2300),
    minEquiRC = rep(1,4),
    gene_id = 1
)

data2 <- data.frame(
    tx_id = c( 1,"1Start", 2, "2Start", 1, 2),
    read_class_id = c(1,"1Start", 2,"2Start", "1.2","1.2"),
    nobs = c(10, 50, 50,10, 500,500),
    tx_len = c(546,546,2356,2356, 546,2356),
    rc_width = c(300,540,1800,2300, 200, 200),
    minEquiRC = rep(1,6),
    gene_id = 2
)

data3 <- data.frame(
    tx_id = c("1Start",1,2, 2, "2Start", 1, 2),
    read_class_id = c("1Start.1.2","1Start.1.2","1Start.1.2",
        2,"2Start", "1.2","1.2"),
    nobs = c(500, 500, 500,10, 0,50,50),
    tx_len = c(546,546,2356,2356,2356, 546,2356),
    rc_width = c(540,540,540,1800,2300, 200, 200),
    minEquiRC = c(NA,NA,1,1,1,NA,1),
    gene_id = 3
)

data4 <- data.frame(
    tx_id = c("1Start",1,2, 2, "2Start", 1, 2),
    read_class_id = c("1Start.1.2","1Start.1.2","1Start.1.2",
        2,"2Start", "1.2","1.2"),
    nobs = c(50, 50, 50,100, 500,20,20),
    tx_len = c(546,546,2356,2356,2356, 546,2356),
    rc_width = c(540,540,540,1800,2300, 200, 200),
    minEquiRC = c(NA,NA,1,1,1,NA,1),
    gene_id = 4
)

data5 <- data.frame(
    tx_id = c(1,"1Start",2, "2Start", "1Start","2Start"),
    read_class_id = c(1,"1Start",
        2,"2Start", "1Start.2Start","1Start.2Start"),
    nobs = c(5, 50, 10,60, 200,200),
    tx_len = c(546,546,2356,2356, 546,2356),
    rc_width = c(1700,2200,1800,2300, 2000, 2000),
    minEquiRC = rep(1,6),
    gene_id = 5
)



estOutput_woBC <- lapply(seq_len(5), function(s) {
    est <- bambu.quantDT(readClassDt = get(paste0("data", s)), 
        emParameters = list(degradationBias = FALSE, 
            maxiter = 10000, conv = 10^(-2), minvalue = 10^(-8)))
})

estOutput_wBC <- lapply(seq_len(5), function(s) {
    est <- bambu.quantDT(readClassDt = get(paste0("data", s)))
})



## expected output for test isore

seReadClass1 <- readRDS(system.file("extdata", "seReadClassUnstranded_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.rds", package = "bambu"))
extendedAnnotations <- readRDS(system.file("extdata", "extendedAnnotationGranges_txdbGrch38_91_chr9_1_1000000.rds", package = "bambu"))
gr <- readRDS(system.file("extdata", "annotationGranges_txdbGrch38_91_chr9_1_1000000.rds", package = "bambu"))


seWithDistExpected <- isore.estimateDistanceToAnnotations(
    seReadClass = seReadClass1,
    annotationGrangesList = extendedAnnotations,
    min.exonDistance = 35
)


## expected seGeneOutput

se <- readRDS(system.file("extdata", "seOutput_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.rds", package = "bambu"))
seExtended <- readRDS(system.file("extdata", "seOutputExtended_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.rds", package = "bambu"))
seCombined <- readRDS(system.file("extdata", "seOutputCombined_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.rds", package = "bambu"))
seCombinedExtended <- readRDS(system.file("extdata", "seOutputCombinedExtended_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.rds", package = "bambu"))


seGeneExpected <- transcriptToGeneExpression(se)
seExtendedGeneExpected <- transcriptToGeneExpression(seExtended)
seCombinedGeneExpected <- transcriptToGeneExpression(seCombined)
seCombinedExtendedGeneExpected <- transcriptToGeneExpression(seCombinedExtended)


usethis::use_data(data1, data2, data3, data4, data5,
    estOutput_woBC,
    estOutput_wBC,
    standardJunctionModels_temp,
    seWithDistExpected,
    seGeneExpected, seExtendedGeneExpected,
    seCombinedGeneExpected, seCombinedExtendedGeneExpected,
    internal = TRUE, overwrite = TRUE
)


## inst data creation
rm(list = ls())
gc()

require(GenomicAlignments) ##readGAlignments
require(AnnotationDbi)#loadDb
require(data.table)#fast large dataset manipulation
require(readxl)


require(ggplot2)
require(RColorBrewer)
require(gridExtra)


cat('Setting working directory')
wkdir <- ''

## get gene List
se <- readRDS("seOutput2020-04-30_updated_wBC.rds")
tx <- rowRanges(se[[1]])
gene <- rowRanges(se[[2]])

geneTx <- as.data.table(mcols(tx))

genecounts <- assays(se[[2]])$counts
txcounts <- assays(se[[1]])$counts

gr <- GRanges(seqnames = "9",
             #ranges = IRanges(1, 1000000),
             ranges = IRanges(1,200000), # for example run in bambu.R only, for demonstration
             strand = "+")

hit <- findOverlaps(gene, gr, ignore.strand = TRUE)


genevec <- names(gene[queryHits(hit)])
txvec <- geneTx[GENEID %in% genevec]$TXNAME

geneList.file <- paste0(wkdir,"/geneList.txt")
write.table(genevec, file = geneList.file, sep = '\t',
            row.names = FALSE, col.names = FALSE)



## get gtf
gtf.file <- "Homo_sapiens.GRCh38.91.gtf"
new_gtf.file <- paste0(wkdir, "/Homo_sapiens.GRCh38.91_chr",as.character(seqnames(gr)),
                       "_",start(gr),"_",end(gr),".gtf")
system2(paste0("grep -f  ",geneList.file," ",gtf.file," > ",new_gtf.file))

## make txdb for chr9
gtf.file <- system.file("extdata", "Homo_sapiens.GRCh38.91_chr9_1_1000000.gtf", package = "bambu")
txdb <- makeTxDbFromGFF(gtf.file, format = "gtf",
dataSource="Homo_sapiens.GRCh38.91_chr9_1_1000000.gtf",
organism="Homo sapiens",
taxonomyId=NA,
chrominfo=NULL,
miRBaseBuild=NA,
metadata=NULL
)
saveDb(txdb, file="./inst/extdata/Homo_sapiens.GRCh38.91.annotations-txdb_chr9_1_1000000.sqlite")




## generate bam file
runname <- "GIS_A549_directRNA_Rep5-Run1"
bamFile <- dir(paste0("/mnt/s3_ontdata.store.genome.sg/Nanopore/03_Mapping/Grch38/minimap2-2.17-directRNA/",runname),
               pattern = ".bam$", full.names = TRUE)
outBam <- paste0(wkdir, "SGNex_A549_directRNA_replicate5_run1.bam")
system2(paste0('samtools-1.8 view -b ',bamFile,' "',
              paste0(as.character(seqnames(gr))),':',start(gr),'-',end(gr), '" > ', outBam))
system2(paste0("samtools-1.8 index ",outBam))

## generate .fa file
fa.file <- paste0(wkdir,"bambu/inst/extdata/Homo_sapiens.GRCh38.dna_sm.primary_assembly_chr9.fa.gz")
outFa <- paste0(wkdir,"Homo_sapiens.GRCh38.dna_sm.primary_assembly_chr9_1_1000000.fa")
system2(paste0("samtools-1.8 faidx ",fa.file," ",paste0(as.character(seqnames(gr))),":",start(gr),"-",end(gr)," > ", outFa))
system2(paste0("zcat ",fa.file," | head -1"))


## generate annotation rds file
txdb <- loadDb(system.file("extdata", "Homo_sapiens.GRCh38.91.annotations-txdb_chr9_1_1000000.sqlite", package = "bambu"))
gr <- prepareAnnotations(txdb)
saveRDS(gr, file = "./inst/extdata/annotationGranges_txdbGrch38_91_chr9_1_1000000.rds", compress = "xz")


## generate readGrgList file
test.bam <- system.file("extdata", "SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.bam", package = "bambu")
readGrgList <- prepareDataFromBam(Rsamtools::BamFile(test.bam))
saveRDS(readGrgList, file = "./inst/extdata/readGrgList_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.rds", compress = "xz")


## generate read class files

readGrgList <- readRDS(system.file("extdata", "readGrgList_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.rds", package = "bambu"))
gr <- readRDS(system.file("extdata", "annotationGranges_txdbGrch38_91_chr9_1_1000000.rds", package = "bambu"))
genomeSequence <- system.file("extdata", "Homo_sapiens.GRCh38.dna_sm.primary_assembly_chr9_1_1000000.fa", package = "bambu")

genomeSequence <- checkInputSequence(genomeSequence)
refSeqLevels <-  GenomeInfoDb::seqlevels(genomeSequence)
if (!all(GenomeInfoDb::seqlevels(readGrgList) %in% refSeqLevels)) {
    message("not all chromosomes from reads present in reference genome 
        sequence, reads without reference chromosome sequence are dropped")
        readGrgList <- GenomeInfoDb::keepSeqlevels(readGrgList,
            value =  refSeqLevels,
            pruning.mode = "coarse")
    # reassign Ids after seqlevels are dropped
    mcols(readGrgList)$id <- seq_along(readGrgList) 
}
if (!all(GenomeInfoDb::seqlevels(annotations) %in% refSeqLevels)) {
    message("not all chromosomes from annotations present in reference genome 
    sequence, annotations without reference chrosomomse sequence are dropped")
    annotations <- GenomeInfoDb::keepSeqlevels(annotations,
        value = refSeqLevels,pruning.mode = "coarse")
}
unlisted_junctions <- unlistIntrons(readGrgList, use.ids = TRUE)
uniqueJunctions <- isore.constructJunctionTables(unlisted_junctions, 
    annotations = gr,genomeSequence)

seReadClassUnstranded <- isore.constructReadClasses(readGrgList = readGrgList,
        unlisted_junctions, 
        uniqueJunctions,
        runName = "SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000",
        annotations = gr,
        stranded = FALSE,
        verbose = FALSE)
GenomeInfoDb::seqlevels(seReadClassUnstranded) <- refSeqLevels
saveRDS(seReadClassUnstranded, file = "./inst/extdata/seReadClassUnstranded_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.rds", compress = "xz")





seReadClassStranded <- isore.constructReadClasses(readGrgList = readGrgList,
        unlisted_junctions, 
        uniqueJunctions,
        runName = "SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000_Stranded",
        annotations = gr,
        stranded = TRUE,
        verbose = FALSE)
GenomeInfoDb::seqlevels(seReadClassStranded) <- refSeqLevels
saveRDS(seReadClassStranded, file = "./inst/extdata/seReadClassStranded_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.rds", compress = "xz")


## generate seIsoReCombined
seReadClass1 <- readRDS(system.file("extdata", "seReadClassUnstranded_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.rds", package = "bambu"))
seReadClass2 <- readRDS(system.file("extdata", "seReadClassStranded_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.rds", package = "bambu"))

seIsoReRef <- isore.combineTranscriptCandidates(readClassSe=seReadClass1,
                                                readClassSeRef = NULL,
                                                stranded = FALSE,
                                                verbose = FALSE)

seIsoReCombined <- isore.combineTranscriptCandidates(readClassSe=seReadClass2,
                                                     readClassSeRef = seIsoReRef,
                                                     stranded = FALSE,
                                                     verbose = FALSE)

saveRDS(seIsoReRef, file = "./inst/extdata/seIsoReRef_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.rds", compress = "xz")

saveRDS(seIsoReCombined, file = "./inst/extdata/seIsoReCombined_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.rds", compress = "xz")


## extendedAnnotations
seIsoReCombined <- readRDS(system.file("extdata", "seIsoReCombined_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.rds", package = "bambu"))
gr <- readRDS(system.file("extdata", "annotationGranges_txdbGrch38_91_chr9_1_1000000.rds", package = "bambu"))

extendedAnnotations <- isore.extendAnnotations(se=seIsoReCombined,
                                                       annotationGrangesList=gr,
                                                       remove.subsetTx=TRUE,
                                                       min.readCount=2,
                                                       min.readFractionByGene=0.05,
                                                       min.sampleNumber=1,
                                                       min.exonDistance=35,
                                                       min.exonOverlap=10,
                                                       prefix='',
                                                       verbose=FALSE)
saveRDS(extendedAnnotations, file = "./inst/extdata/extendedAnnotationGranges_txdbGrch38_91_chr9_1_1000000.rds", compress = "xz")


## generate se output
test.bam <- system.file("extdata", "SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.bam", package = "bambu")
fa.file <- system.file("extdata", "Homo_sapiens.GRCh38.dna_sm.primary_assembly_chr9_1_1000000.fa", package = "bambu")
gr <- readRDS(system.file("extdata", "annotationGranges_txdbGrch38_91_chr9_1_1000000.rds", package = "bambu"))


set.seed(1234)
seOutput = bambu(reads = test.bam, annotations = gr, genome = fa.file, opt.em = list(degradationBias = FALSE), discovery = FALSE)
saveRDS(seOutput, file = "./inst/extdata/seOutput_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.rds", compress = "xz")

# set.seed(1234)
# seOutputExtended = bambu(reads = test.bam,  annotations =  gr, genome = fa.file, opt.em = list(bias = FALSE), discovery = TRUE)
# saveRDS(seOutputExtended, file = "./inst/extdata/seOutputExtended_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.rds", compress = "xz")


set.seed(1234)
seOutputCombined = bambu(reads =  Rsamtools::BamFileList(c(test.bam, test.bam), yieldSize = 1000),  annotations =  gr, genome = fa.file, discovery = FALSE)
saveRDS(seOutputCombined, file = "./inst/extdata/seOutputCombined_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.rds", compress = "xz")


set.seed(1234)
seOutputCombinedExtended = bambu(reads =  Rsamtools::BamFileList(c(test.bam, test.bam), yieldSize = 1000),  annotations =  gr, genome = fa.file, discovery = TRUE)
saveRDS(seOutputCombinedExtended, file = "./inst/extdata/seOutputCombinedExtended_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.rds", compress = "xz")


seReadClass1 <- system.file("extdata", "seReadClassUnstranded_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.rds", package = "bambu")
gr <- readRDS(system.file("extdata", "annotationGranges_txdbGrch38_91_chr9_1_1000000.rds", package = "bambu"))
set.seed(1234)
seOutputExtended <- bambu(rcFile = seReadClass1, annotations = gr, opt.em = list(degradationBias = FALSE), discovery = TRUE)
saveRDS(seOutputExtended, file = "./inst/extdata/seOutputExtended_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.rds", compress = "xz")


set.seed(1234)
seOutputCombined2 <- bambu(rcFile = c(seReadClass1, seReadClass1), annotations = gr, discovery = FALSE)
saveRDS(seOutputCombined2, file = "./inst/extdata/seOutputCombined2_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.rds", compress = "xz")

