context("txRange, full length transcript prediction")

test_that("txRange generates a gene and transcript score",{
    se <- readRDS(system.file("extdata", "test_se.rds", 
        package = "bambu"))
    genomeSequence <- readRDS(system.file("extdata", 
        "test_genomeSequence.rds", package = "bambu"))
    annotations <- readRDS(system.file("extdata", 
        "test_annotations.rds", package = "bambu"))
    se = scoreReadClasses(se, genomeSequence, annotations, defaultModels)

    seExpected = readRDS(system.file("extdata", "test_se_scored.rds", 
        package = "bambu"))

    expect_is(rowData(se)$txScore, class = 'numeric')
    expect_equal(rowData(se)$txScore, rowData(seExpected)$txScore)
})


test_that("calculateGeneProportion()",{
    se <- readRDS(system.file("extdata", "test_se.rds", package = "bambu"))
    annotations <- readRDS(system.file("extdata", "test_annotations.rds", 
        package = "bambu"))
    seExpected = readRDS(system.file("extdata", "test_se_scored.rds", 
        package = "bambu"))
    rowData(se)$GENEID = assignGeneIds(rowRanges(se), annotations)
    countsTBL = calculateGeneProportion(counts=mcols(se)$readCount,
                                        geneIds=mcols(se)$GENEID)
    expect_equal(countsTBL$geneReadProp, 
        rowData(seExpected)$geneReadProp)
})

test_that("isReadClassCompatible() classifies RCs correctly",{
    se <- readRDS(system.file("extdata", "test_se.rds", package = "bambu"))
    annotations <- readRDS(system.file("extdata", "test_annotations.rds", 
        package = "bambu"))
    seExpected = readRDS(system.file("extdata", "test_se_scored.rds", 
        package = "bambu"))
    thresholdIndex = which(rowData(se)$readCount>=2)
    compTable = isReadClassCompatible(rowRanges(se[thresholdIndex,]), annotations)
    newRowData = data.frame(equal = compTable$equal,
                            compatible = compTable$compatible)
    rowData(se)[names(newRowData)] = NA
    rowData(se)[thresholdIndex,names(newRowData)] = newRowData
    expect_equal(rowData(se)$equal,rowData(seExpected)$equal)
    expect_equal(rowData(se)$compatible, rowData(seExpected)$compatible)
})

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


test_that("checkFeatures() detects insufficient samples",{
    se <- readRDS(system.file("extdata", "test_se_scored.rds", package = "bambu"))
    thresholdIndex = which(rowData(se)$readCount>=2)
    features <- prepareTranscriptModelFeatures(rowData(se)[thresholdIndex,])
    trainable = checkFeatures(features)
    #normal se
    expect_equal(trainable, FALSE)
    #TODO se has all TRUE labels
    expect_equal(trainable, FALSE)
    #TODOse with no TRUE labels
    expect_equal(trainable, FALSE)
    #TODO se with not enough points
    expect_equal(trainable, FALSE)
    #TODO se with not enough true/false labels
    expect_equal(trainable, FALSE)
})

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
