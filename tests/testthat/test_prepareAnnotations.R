context("Prepare annotations")

test_that("prepareAnnotations of txdb object is a GRangesList", {
    txdb <- AnnotationDbi::loadDb(system.file("extdata", "Homo_sapiens.GRCh38.91.annotations-txdb_chr9_1_1000000.sqlite", package = "bambu"))

    expectedGR <- readRDS(system.file("extdata", "annotationGranges_txdb2Grch38_91_chr9_1_1000000.rds", package = "bambu"))

    gr <- prepareAnnotations(x = txdb)

    expect_equal(gr, expectedGR)
    expect_s4_class(gr, class = "CompressedGRangesList")
    expect_named(mcols(gr), c("TXNAME", "GENEID", "txid", "eqClass", "eqClassById"))
})


test_that("prepareAnnotations of genome library is a GRangesList", {
    library(TxDb.Hsapiens.UCSC.hg38.knownGene)
    gr <- prepareAnnotations(TxDb.Hsapiens.UCSC.hg38.knownGene)

    expect_s4_class(gr, class = "CompressedGRangesList")
    expect_named(mcols(gr), c("TXNAME", "GENEID", "txid", "eqClass", "eqClassById"))
})


test_that("prepareAnnotations from gtf file is GRangesList", {
    gtf.file <- system.file("extdata", "Homo_sapiens.GRCh38.91_chr9_1_1000000.gtf", package = "bambu")
    
    expectedGR <- readRDS(system.file("extdata", "annotationGranges_txdbGrch38_91_chr9_1_1000000.rds", package = "bambu"))

    gr <- prepareAnnotations(x = gtf.file)

    expect_equal(gr[order(names(gr))], expectedGR[order(names(expectedGR))])
    expect_s4_class(gr, class = "CompressedGRangesList")
    expect_named(mcols(gr), c("TXNAME", "GENEID", "txid", "eqClass", "eqClassById"))
})


test_that("positive/negative strand gives ascending/descending exon_rank and descending/ascending exon_endRank", {
    gtf.file <- system.file("extdata", "Homo_sapiens.GRCh38.91_chr9_1_1000000.gtf", package = "bambu")
    gr <- prepareAnnotations(x = gtf.file)

    check <- all(sapply(seq_along(gr), function(i){
        if (runValue(strand(gr[[i]])) == "-"){
            # Check whether the ranking follows the direction of strand
            q1 <- all(mcols(gr[[i]])$exon_rank == sort(seq_along(gr[[i]]), decreasing=TRUE))
            q2 <- all(mcols(gr[[i]])$exon_endRank == seq_along(gr[[i]]))
            return(q1 & q2)

        } else {
            q1 <- all(mcols(gr[[i]])$exon_endRank == sort(seq_along(gr[[i]]), decreasing=TRUE))
            q2 <- all(mcols(gr[[i]])$exon_rank == seq_along(gr[[i]]))
            return(q1 & q2)
        }
    }))

    expect_true(check)
})


test_that("txid is a subset of eqClassById", {
    gtf.file <- system.file("extdata", "Homo_sapiens.GRCh38.91_chr9_1_1000000.gtf", package = "bambu")
    gr <- prepareAnnotations(x = gtf.file)
    check <- all(sapply(mcols(gr)$txid, function(x){
        x %in% mcols(gr)$eqClassById[[x]]}))
    expect_true(check)      
})


test_that("eqClassById is as it claims", {
    gtf.file <- system.file("extdata", "Homo_sapiens.GRCh38.91_chr9_1_1000000.gtf", package = "bambu")
    gr <- prepareAnnotations(x=gtf.file)
    
    # Store the junction coordinates of each transcript
    junc_coord <- lapply(seq_along(gr), function(i){
        start_coord <- start(gr[[i]])
        end_coord <- end(gr[[i]])
        tx_coord <- c(rbind(start_coord, end_coord))
        
        if (length(gr[[i]]) == 1){
            return(tx_coord)
        } else {
            return(tx_coord[-c(1,length(tx_coord))])
        }
    })
    
    # Check whether eqClassById matches
    txIntx <- function(txid){
        txlist <- c()
        for (i in seq_along(gr)){
            if (grepl(paste(junc_coord[[txid]], collapse="."), paste(junc_coord[[i]], collapse="."))){
                txlist <- c(txlist, i)
            }
        }
        return(all(txlist == mcols(gr)$eqClassById[[txid]]))
    }
  
    check <- all(sapply(seq_along(gr), txIntx))
    expect_true(check)
})


