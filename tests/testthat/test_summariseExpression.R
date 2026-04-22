context("Summarise expression")

test_that("summariseByExpression correctly reduces transcript counts to exon and gene counts", {
    # TX1: exons at 1-100, 200-300 -> GENE1, count = 10
    # TX2: exons at 1-100, 400-500 -> GENE1, shares first exon with TX1, count = 5
    # TX3: exon at 600-700         -> GENE2, count = 8
    #
    # Expected exon counts:
    #   chr1:1-100   (TX1 + TX2) -> 15
    #   chr1:200-300 (TX1 only)  -> 10
    #   chr1:400-500 (TX2 only)  ->  5
    #   chr1:600-700 (TX3 only)  ->  8
    #
    # Expected gene counts:
    #   GENE1 (TX1 + TX2) -> 15
    #   GENE2 (TX3 only)  ->  8

    txRanges <- GRangesList(
        TX1 = GRanges("chr1", IRanges(c(1, 200), c(100, 300)), strand = "+"),
        TX2 = GRanges("chr1", IRanges(c(1, 400), c(100, 500)), strand = "+"),
        TX3 = GRanges("chr1", IRanges(600, 700), strand = "+")
    )
    mcols(txRanges)$TXNAME             <- c("TX1", "TX2", "TX3")
    mcols(txRanges)$GENEID             <- c("GENE1", "GENE1", "GENE2")
    mcols(txRanges)$novelTranscript    <- c(FALSE, FALSE, FALSE)
    mcols(txRanges)$txClassDescription <- c("annotation", "annotation", "annotation")

    counts <- matrix(c(10, 5, 8), nrow = 3, ncol = 1,
        dimnames = list(c("TX1", "TX2", "TX3"), "sample1"))

    se <- SummarizedExperiment(
        assays    = list(counts = counts),
        rowRanges = txRanges
    )

    # exon level
    seExon <- summariseByExpression(se, type = "exon")
    exonCounts <- as.matrix(assays(seExon)$counts)
    exonRanges <- rowRanges(seExon)

    shared  <- which(start(exonRanges) == 1   & end(exonRanges) == 100)
    tx1only <- which(start(exonRanges) == 200 & end(exonRanges) == 300)
    tx2only <- which(start(exonRanges) == 400 & end(exonRanges) == 500)
    tx3only <- which(start(exonRanges) == 600 & end(exonRanges) == 700)

    expect_equal(exonCounts[shared,  1], 15)
    expect_equal(exonCounts[tx1only, 1], 10)
    expect_equal(exonCounts[tx2only, 1],  5)
    expect_equal(exonCounts[tx3only, 1],  8)

    # gene level
    seGene <- summariseByExpression(se, type = "gene")
    geneCounts <- as.matrix(assays(seGene)$counts)

    expect_equal(geneCounts["GENE1", 1], 15)
    expect_equal(geneCounts["GENE2", 1],  8)
})
