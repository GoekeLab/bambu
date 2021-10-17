context("txRange, full length transcript prediction")

test_that("txRange generates a gene and transcript score",{
    se <- readRDS(system.file("extdata", "test_se.rds", package = "bambu"))
    genomeSequence <- readRDS(system.file("extdata", "test_genomeSequence.rds", 
        package = "bambu"))
    annotations <- readRDS(system.file("extdata", "test_annotations.rds", 
        package = "bambu"))
    se = scoreReadClasses(se, genomeSequence, annotations)
    expect_is(rowData(se)$geneScore, class = 'numeric')
    expect_is(rowData(se)$txScore, class = 'numeric')
})

# test_that("addRowData adds all the correct rowData features",{
#     se <- readRDS(system.file("extdata", "test_se.rds", package = "bambu"))
#     genomeSequence <- readRDS(system.file("extdata", "test_genomeSequence.rds", 
#         package = "bambu"))
#     annotations <- readRDS(system.file("extdata", "test_annotations.rds", 
#         package = "bambu")) 
#     se = addRowData(se, genomeSequence, annotations)
#     expect_is(rowData(se)$numExons, class = 'numeric')
#     expect_is(rowData(se)$equal, class = 'logical')
#     expect_is(rowData(se)$GENEID, class = 'character')
#     expect_is(rowData(se)$novelGene, class = 'logical')
#     expect_is(rowData(se)$totalGeneReadProp, class = 'numeric')
#     expect_is(rowData(se)$numAstart, class = 'integer')
#     expect_is(rowData(se)$numAend, class = 'integer')
#     expect_is(rowData(se)$numTstart, class = 'integer')
#     expect_is(rowData(se)$numTend, class = 'integer')
# })

test_that("calculateGeneProportion",{
    se <- readRDS(system.file("extdata", "test_se.rds", package = "bambu"))
    annotations <- readRDS(system.file("extdata", "test_annotations.rds", 
        package = "bambu"))
    seExpected = readRDS(system.file("extdata", 
        "SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000_readClassSe.rds", 
        package = "bambu"))
    rowData(se)$GENEID = assignGeneIds(rowRanges(se), annotations)
    se = calculateGeneProportion(se)
    expect_equal(rowData(se)$totalGeneReadProp, 
        rowData(seExpected)$totalGeneReadProp)
})

# test_that("isReadClassEqual",{
#     se <- readRDS(system.file("extdata", "test_se.rds", package = "bambu"))
#     annotations <- readRDS(system.file("extdata", "test_annotations.rds", 
#         package = "bambu"))
#     seExpected = readRDS(system.file("extdata", 
#         "SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000_readClassSe.rds", 
#         package = "bambu"))
#     expect_equal(isReadClassEqual(rowRanges(se), annotations),
#         rowData(seExpected)$equal)
# })

test_that("countPolyATerminals",{
    se <- readRDS(system.file("extdata", "test_se.rds", package = "bambu"))
    genomeSequence <- readRDS(system.file("extdata", "test_genomeSequence.rds", 
        package = "bambu"))
    seExpected = readRDS(system.file("extdata", 
        "SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000_readClassSe.rds", 
        package = "bambu"))
    se = countPolyATerminals(se, genomeSequence)
    expect_equal(rowData(se)$numAstart,rowData(seExpected)$numAstart)
    expect_equal(rowData(se)$numAend,rowData(seExpected)$numAend)
    expect_equal(rowData(se)$numTstart,rowData(seExpected)$numTstart)
    expect_equal(rowData(se)$numTend,rowData(seExpected)$numTend)
})


# ## same for this one, is it still existing in the code, can't seem to find it
# test_that("getAgnosticFeatures",{
#     #TODO make a expected features dataset to compare to
#     se <- readRDS(system.file("extdata", "test_se.rds", package = "bambu"))
#     features <- as.data.frame(getAgnosticFeatures(se))
#     expect_is(features$numReads, class = 'numeric')
#     expect_is(features$SD, class = 'numeric')
#     expect_is(features$SDend, class = 'numeric')
#     expect_is(features$geneReadProp, class = 'numeric')
#     expect_is(features$tx_strand_bias, class = 'numeric')
#     expect_is(features$numAstart, class = 'numeric')
#     expect_is(features$numAend, class = 'numeric')
#     expect_is(features$numTstart, class = 'numeric')
#     expect_is(features$numTend, class = 'numeric')
#     expect_is(features$transcriptProp, class = 'numeric')
#     expect_is(features$transcriptPropMin, class = 'numeric')
# })

test_that("assignGeneIds",{
    se <- readRDS(system.file("extdata", "test_se.rds", package = "bambu"))
    annotations <- readRDS(system.file("extdata", "test_annotations.rds", 
        package = "bambu"))
    seExpected = readRDS(system.file("extdata", 
        "SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000_readClassSe.rds", 
        package = "bambu"))
    expect_equal(assignGeneIds(rowRanges(se), annotations),
        rowData(seExpected)$GENEID)
})


# ## is this still use?? should it be prepareTranscriptModelFeatures??
# test_that("prepareGeneModelFeatures",{
#     se <- readRDS(system.file("extdata", "test_se.rds", package = "bambu"))
#     geneFeatures <- as.data.frame(prepareGeneModelFeatures(se))
#     expect_is(features$numReadsLog, class = 'numeric')
#     expect_is(features$strand_bias, class = 'numeric')
#     expect_is(features$numRCs, class = 'numeric')
#     expect_is(features$numNonSubsetRCs, class = 'numeric')
#     expect_is(features$numExons, class = 'numeric')
#     expect_is(features$isSpliced, class = 'numeric')
#      expect_is(features$highConfidence, class = 'numeric')
# })
