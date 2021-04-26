context("txRange, full length transcript prediction")

test_that("scoreReadClasses adds the correct rowData",{
  se<- readRDS(system.file("extdata", "seReadClassUnstranded_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.rds", package = "bambu"))
  gr <- readRDS(system.file("extdata", "annotationGranges_txdbGrch38_91_chr9_1_1000000.rds", package = "bambu"))
  genomeSequence <- system.file("extdata","Homo_sapiens.GRCh38.dna_sm.primary_assembly_chr9_1_1000000.fa",package = "bambu")
  genomeSequence <- checkInputSequence(genomeSequence)
  se = scoreReadClasses(se, genomeSequence, gr)
  expected = read.table(system.file("extdata", "rowData_seOutput_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.csv", package = "bambu"), header = TRUE)
  expect_equal(rowData(se)$numExons, expected$numExons)
  expect_equal(rowData(se)$novel, expected$novel)
})

test_that("calculateGeneProportion",{
  se<- readRDS(system.file("extdata", "seReadClassUnstranded_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.rds", package = "bambu"))
  gr <- readRDS(system.file("extdata", "annotationGranges_txdbGrch38_91_chr9_1_1000000.rds", package = "bambu"))
  rowData(se)$GENEID = assignGeneIds(rowRanges(se), gr)
  countsTBL = calculateGeneProportion(counts=mcols(se)$readCount,
                               geneIds=mcols(se)$GENEID)
  expected = read.table(system.file("extdata", "rowData_seOutput_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.csv", package = "bambu"), header = TRUE)
  expect_equal(countsTBL$geneReadProp, expected$geneReadProp)
  expect_equal(countsTBL$geneReadCount, expected$geneReadCount)
})

test_that("isReadClassCompatible classifies read classes correctly",{
  se<- readRDS(system.file("extdata", "seReadClassUnstranded_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.rds", package = "bambu"))
  gr <- readRDS(system.file("extdata", "annotationGranges_txdbGrch38_91_chr9_1_1000000.rds", package = "bambu"))
  expected = read.table(system.file("extdata", "rowData_seOutput_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.csv", package = "bambu"), header = TRUE)
  thresholdIndex = which(rowData(se)$readCount >= 2)
  compTable <- isReadClassCompatible(rowRanges(se[thresholdIndex,]), 
                                     gr)
  equal = compatible = rep(NA, nrow(se))
  equal[thresholdIndex] = compTable$equal
  compatible[thresholdIndex] = compTable$compatible
  
  expect_equal(equal, expected$equal)
  expect_equal(compatible, expected$compatible)
})

test_that("countPolyATerminals returns the correct number of bases",{
  se<- readRDS(system.file("extdata", "seReadClassUnstranded_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.rds", package = "bambu"))
  genomeSequence <- system.file("extdata","Homo_sapiens.GRCh38.dna_sm.primary_assembly_chr9_1_1000000.fa",package = "bambu")
  genomeSequence <- checkInputSequence(genomeSequence)
  expected = read.table(system.file("extdata", "rowData_seOutput_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.csv", package = "bambu"), header = TRUE)
  thresholdIndex = which(rowData(se)$readCount >= 2)
  polyATerminals = countPolyATerminals(rowRanges(se[thresholdIndex,]), 
                                       genomeSequence)
  numAstart = numAend = numTstart = numTend = rep(NA, nrow(se))
  numAstart[thresholdIndex] = polyATerminals$numAstart 
  numAend[thresholdIndex] = polyATerminals$numAend
  numTstart[thresholdIndex] = polyATerminals$numTstart
  numTend[thresholdIndex] = polyATerminals$numTend
  
  expect_equal(numAstart,expected$numAstart)
  expect_equal(numAend,expected$numAend)
  expect_equal(numTstart,expected$numTstart)
  expect_equal(numTend,expected$numTend)
})

test_that("assignGeneId returns correct gene ids",{
  se <- readRDS(system.file("extdata", "seReadClassUnstranded_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.rds", package = "bambu"))
  gr <- readRDS(system.file("extdata", "annotationGranges_txdbGrch38_91_chr9_1_1000000.rds", package = "bambu"))
  expected = read.table(system.file("extdata", "rowData_seOutput_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.csv", package = "bambu"), header = TRUE)
  expect_equal(assignGeneIds(rowRanges(se), gr), expected$GENEID)
})

test_that('Calculated FDR scores are as expected', {
  scores = as.numeric(readLines(system.file('extdata','calculateFDR_arg_scores.txt', package='bambu')))
  labels = as.numeric(readLines(system.file('extdata','calculateFDR_arg_labels.txt', package='bambu')))
  FDR = calculateFDR(scores, labels)
  
  ref_FDR = unname(unlist(read.csv(system.file('extdata','FDR.csv', package='bambu'))))
  
  expect_equal(FDR, ref_FDR)
})