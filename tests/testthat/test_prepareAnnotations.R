context("Prepare annotations")

test_that("prepareAnnotations of txdb object is a GRangesList", {
    txdb <- AnnotationDbi::loadDb(system.file("extdata", "Homo_sapiens.GRCh38.91.annotations-txdb_chr9_1_1000000.sqlite", package = "bambu"))

    expectedGR <- readRDS(system.file("extdata", "annotationGranges_txdbGrch38_91_chr9_1_1000000.rds", package = "bambu"))

    gr <- prepareAnnotations(x = txdb)

    expect_equal(gr, expectedGR)
    expect_s4_class(gr, class = "CompressedGRangesList")
    expect_named(mcols(gr), c("TXNAME", "GENEID", "eqClass"))
})


test_that("prepareAnnotations of genome library is a GRangesList", {
    library(TxDb.Hsapiens.UCSC.hg38.knownGene)
    gr <- prepareAnnotations(TxDb.Hsapiens.UCSC.hg38.knownGene)

    expect_s4_class(gr, class = "CompressedGRangesList")
    expect_named(mcols(gr), c("TXNAME", "GENEID", "eqClass"))
})




test_that("prepareAnnotations from gtf file is GRangesList", {
    gtf.file <- system.file("extdata", "Homo_sapiens.GRCh38.91_chr9_1_1000000.gtf", package = "bambu")

    expectedGR <- readRDS(system.file("extdata", "annotationGranges_txdbGrch38_91_chr9_1_1000000.rds", package = "bambu"))

    gr <- prepareAnnotations(x = gtf.file)

    expect_equal(gr[order(names(gr))], expectedGR[order(names(expectedGR))])
    expect_s4_class(gr, class = "CompressedGRangesList")
    expect_named(mcols(gr), c("TXNAME", "GENEID", "eqClass"))
})
