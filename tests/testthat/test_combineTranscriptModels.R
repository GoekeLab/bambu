context("Combine transcript models across samples")

## Builds a minimal read class SE carrying exactly the columns
## extractFeaturesFromReadClassSE() consumes. Synthetic rather than a stored
## fixture so the sample set can be varied per test, and so the intron chains
## overlap across samples, which is what the combine has to reduce over.
makeReadClassSe <- function(nClass, seed, nShared = 6L) {
    set.seed(seed)
    ## The first nShared chains are the same in every sample; the rest are
    ## sample specific. That gives both keys seen once and keys seen many times.
    chainId <- c(seq_len(nShared), seed * 1000L + seq_len(nClass - nShared))
    starts <- 1000L * chainId
    exonsByRc <- GenomicRanges::GRangesList(lapply(seq_len(nClass), function(i) {
        GenomicRanges::GRanges(seqnames = "chr9",
            ranges = IRanges::IRanges(
                start = c(starts[i], starts[i] + 400L),
                end = c(starts[i] + 200L, starts[i] + 600L)),
            strand = "+")
    }))
    rd <- S4Vectors::DataFrame(
        chr.rc = factor(rep("chr9", nClass), levels = "chr9"),
        strand.rc = factor(rep("+", nClass), levels = c("+", "-", "*")),
        intronStarts = as.character(starts + 201L),
        intronEnds = as.character(starts + 399L),
        confidenceType = rep("highConfidenceJunctionReads", nClass),
        readCount = sample(seq_len(50L), nClass, replace = TRUE),
        geneReadProp = runif(nClass),
        txScore = runif(nClass),
        txScore.noFit = runif(nClass),
        numExons = rep(2L, nClass))
    ## A missing score has to survive the reduction as NA, so seed some.
    rd$txScore[sample(nClass, max(1L, nClass %/% 10L))] <- NA_real_
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = matrix(rd$readCount, ncol = 1)),
        rowRanges = exonsByRc)
    SummarizedExperiment::rowData(se) <- rd
    se
}

test_that("upperMedianByGroup reproduces readCountWeightedMedian", {
    set.seed(1)
    for (trial in seq_len(50)) {
        ngroup <- sample(seq_len(6), 1)
        n <- sample(seq(ngroup, 40), 1)
        groupIndex <- sort(sample(seq_len(ngroup), n, replace = TRUE))
        groupIndex <- match(groupIndex, sort(unique(groupIndex)))
        values <- sample(seq_len(500), n, replace = TRUE)
        weights <- sample(seq_len(20), n, replace = TRUE)
        got <- upperMedianByGroup(groupIndex, values, weights,
            max(groupIndex))
        want <- vapply(seq_len(max(groupIndex)), function(g) {
            keep <- groupIndex == g
            dt <- data.table(v = values[keep], w = weights[keep])
            readCountWeightedMedian(dt, "v", "w")
        }, numeric(1))
        expect_equal(got, want)
    }
})

test_that("the sparse combine is identical to the dense combine", {
    ## Three samples is the smallest set that takes the sparse path.
    for (nSample in c(3L, 5L)) {
        readClassList <- lapply(seq_len(nSample), function(i)
            makeReadClassSe(20L + i, seed = i))
        bpParameters <- BiocParallel::SerialParam()
        seed <- 42L
        ## combineSplicedTranscriptModels draws from the RNG itself, so both
        ## paths have to start from the same seed to see the same grouping.
        set.seed(seed)
        sparse <- combineSplicedTranscriptModels(readClassList, bpParameters,
            min.readCount = 2, min.readFractionByGene = 0.05,
            min.txScore.multiExon = 0, min.txScore.singleExon = 1,
            verbose = FALSE)
        set.seed(seed)
        indexList <- sample(rep(seq_len(max(ceiling(nSample / 10),
            min(BiocParallel::bpworkers(bpParameters),
                round(nSample / 2)))), length.out = nSample))
        indexList <- splitAsList(seq_len(nSample), indexList)
        dense <- combineSplicedTranscriptModelsDense(readClassList, indexList,
            bpParameters, min.readCount = 2, min.readFractionByGene = 0.05,
            min.txScore.multiExon = 0, min.txScore.singleExon = 1)
        expect_identical(setDT(sparse), setDT(dense))
    }
})

test_that("the sparse combine defers when a sample repeats an intron chain", {
    readClassList <- lapply(seq_len(3), function(i)
        makeReadClassSe(20L, seed = i))
    ## Duplicating a key inside one sample would fan out in the dense join, so
    ## the sparse path must decline rather than return a different answer.
    rd <- SummarizedExperiment::rowData(readClassList[[1]])
    rd$intronStarts[2] <- rd$intronStarts[1]
    rd$intronEnds[2] <- rd$intronEnds[1]
    SummarizedExperiment::rowData(readClassList[[1]]) <- rd
    indexList <- splitAsList(seq_len(3), c(1L, 1L, 2L))
    expect_null(combineSplicedTranscriptModelsSparse(readClassList, indexList,
        BiocParallel::SerialParam(), min.readCount = 2,
        min.readFractionByGene = 0.05, min.txScore.multiExon = 0,
        min.txScore.singleExon = 1))
})
