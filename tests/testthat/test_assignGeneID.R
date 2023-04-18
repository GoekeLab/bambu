context("bambu-processReads_utilityConstructReadClasses")

test_that("assignGeneIds",{
    se <- readRDS(system.file("extdata",
        "seReadClassUnstranded_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.rds",
        package = "bambu"))
    annotations <- readRDS(system.file("extdata","annotationGranges_txdbGrch38_91_chr9_1_1000000.rds",package = "bambu"))

    geneIds = assignGeneIds(rowRanges(se), annotations)
    expect_equal(geneIds$GENEID,
               rowData(se)$GENEID)
    expect_equal(geneIds$novelGene,
               rowData(se)$novelGene)
})
