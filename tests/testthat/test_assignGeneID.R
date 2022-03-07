context("bambu-processReads_utilityConstructReadClasses")

test_that("assignGeneIds",{
  se <- readRDS(system.file("extdata", "test_se.rds", package = "bambu"))
  annotations <- readRDS(system.file("extdata", "test_annotations.rds", 
                                     package = "bambu"))
  seExpected = readRDS(system.file("extdata", "test_se_scored.rds", 
                                   package = "bambu"))
  expect_equal(assignGeneIds(rowRanges(se), annotations)$geneIds,
               rowData(seExpected)$GENEID)
})
