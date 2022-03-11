context("bambu-processReads_utilityConstructReadClasses")

test_that("assignGeneIds",{
  se <- readRDS(system.file("extdata", "test_se.rds", package = "bambu"))
  annotations <- readRDS(system.file("extdata", "test_annotations.rds", 
                                     package = "bambu"))
  seExpected = readRDS(system.file("extdata", "test_se_scored.rds", 
                                   package = "bambu"))
    geneIds = assignGeneIds(rowRanges(se), annotations)
  expect_equal(geneIds$GENEID,
               rowData(seExpected)$GENEID)
  expect_equal(geneIds$novelGene,
               rowData(seExpected)$novelGene)
})
