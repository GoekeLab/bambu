context("txRange, full length transcript prediction")

test_that("txRange generates a gene and transcript score",{
    readClasses <- readRDS(system.file("extdata","readClassesUnstranded_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.rds",package = "bambu"))
    annotations <- readRDS(system.file("extdata","annotationGranges_txdbGrch38_91_chr9_1_1000000.rds",package = "bambu"))
    genomeSequence <- system.file("extdata","Homo_sapiens.GRCh38.dna_sm.primary_assembly_chr9_1_1000000.fa",package = "bambu")
    genomeSequence <- checkInputSequence(genomeSequence)
    
    se = scoreReadClasses(readClasses, genomeSequence, annotations, defaultModels)
  
    seExpected = readRDS(system.file("extdata", "seReadClassUnstranded_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.rds", package = "bambu"))
  
    expect_is(rowData(se)$txScore, class = 'numeric')
    expect_equal(rowData(se)$txScore, rowData(seExpected)$txScore)
})


test_that("calculateGeneProportion()",{
    readClasses <- readRDS(system.file("extdata","readClassesUnstranded_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.rds",package = "bambu"))
    annotations <- readRDS(system.file("extdata","annotationGranges_txdbGrch38_91_chr9_1_1000000.rds",package = "bambu"))
    seExpected = readRDS(system.file("extdata", "seReadClassUnstranded_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.rds", package = "bambu"))
  
    geneIds = assignGeneIds(rowRanges(readClasses), annotations)
    rowData(readClasses)[,names(geneIds)] = geneIds
    countsTBL = calculateGeneProportion(counts=mcols(readClasses)$readCount,
                                      geneIds=mcols(readClasses)$GENEID)
    expect_equal(countsTBL$geneReadProp, 
               rowData(seExpected)$geneReadProp)
})

test_that("isReadClassCompatible() classifies RCs correctly",{
    readClasses <- readRDS(system.file("extdata","readClassesUnstranded_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.rds",package = "bambu"))
    annotations <- readRDS(system.file("extdata","annotationGranges_txdbGrch38_91_chr9_1_1000000.rds",package = "bambu"))
    seExpected = readRDS(system.file("extdata", "seReadClassUnstranded_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.rds", package = "bambu"))
    thresholdIndex = which(rowData(readClasses)$readCount>=2)
    compTable = isReadClassCompatible(rowRanges(readClasses[thresholdIndex,]), annotations)
    newRowData = data.frame(equal = compTable$equal,
                          compatible = compTable$compatible)
    rowData(readClasses)[names(newRowData)] = NA
    rowData(readClasses)[thresholdIndex,names(newRowData)] = newRowData
  expect_equal(rowData(readClasses)$equal,rowData(seExpected)$equal)
  expect_equal(rowData(readClasses)$compatible, rowData(seExpected)$compatible)
})

test_that("countPolyATerminals",{
    readClasses <- readRDS(system.file("extdata","readClassesUnstranded_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.rds",package = "bambu"))
    genomeSequence <- system.file("extdata","Homo_sapiens.GRCh38.dna_sm.primary_assembly_chr9_1_1000000.fa",package = "bambu")
    genomeSequence <- checkInputSequence(genomeSequence)
    seExpected = readRDS(system.file("extdata", "seReadClassUnstranded_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.rds", package = "bambu"))
    thresholdIndex = which(rowData(readClasses)$readCount>=2)
    polyATerminals = countPolyATerminals(rowRanges(readClasses[thresholdIndex,]), genomeSequence)
    newRowData = data.frame(numAstart = polyATerminals$numAstart,
                        numAend = polyATerminals$numAend,
                        numTstart = polyATerminals$numTstart,
                        numTend = polyATerminals$numTend)
    rowData(readClasses)[names(newRowData)] = NA
    rowData(readClasses)[thresholdIndex,names(newRowData)] = newRowData
  
    expect_equal(rowData(readClasses)$numAstart,rowData(seExpected)$numAstart)
    expect_equal(rowData(readClasses)$numAend,rowData(seExpected)$numAend)
    expect_equal(rowData(readClasses)$numTstart,rowData(seExpected)$numTstart)
    expect_equal(rowData(readClasses)$numTend,rowData(seExpected)$numTend)
})

test_that("getTranscriptScore() calculates the correct score", {
    seExpected = readRDS(system.file("extdata", "seReadClassUnstranded_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.rds", package = "bambu"))
    se=seExpected
    rowData(se)$txScore = NULL
    thresholdIndex = which(rowData(se)$readCount >= 2)
    txScore = getTranscriptScore(rowData(se)[thresholdIndex,], model = NULL, defaultModels)
    rowData(se)$txScore = rep(NA,nrow(se))
    rowData(se)$txScore[thresholdIndex] = txScore
    expect_equal(rowData(se)$txScore, rowData(seExpected)$txScore)
})

test_that("prepareTranscriptModelFeatures() produces the correct data types",{
  #TODO make a expected features dataset to compare to
    se = readRDS(system.file("extdata", "seReadClassUnstranded_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.rds", package = "bambu"))
    features <- prepareTranscriptModelFeatures(rowData(se))
    expect_is(features$numReads, class = 'numeric')
    expect_is(features$geneReadProp, class = 'numeric')
    expect_is(features$startSD, class = 'numeric')
    expect_is(features$endSD, class = 'numeric')
    expect_is(features$numAstart, class = 'integer')
    expect_is(features$numAend, class = 'integer')
    expect_is(features$numTstart, class = 'integer')
    expect_is(features$numTend, class = 'integer')
    expect_is(features$tx_strand_bias, class = 'numeric')
})

test_that("checkFeatures() detects insufficient samples",{
    se = readRDS(system.file("extdata", "seReadClassUnstranded_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.rds", package = "bambu"))
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
