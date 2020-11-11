context("Visualization of estimates")


test_that("visualization for transcript expression is successful", {
    seCombined <- readRDS(system.file("extdata", "seOutputCombined_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.rds", package = "bambu"))

    colnames(seCombined) <- colData(seCombined)$name <- c("sample1", "sample2")
    set.seed(1)
    assays(seCombined)$CPM[, 2] <- pmax(0, rnorm(length(assays(seCombined)$CPM[, 2]), assays(seCombined)$CPM[, 2], 10))

    # case 1: no group.variable provided
    expect_is(plotBambu(seCombined, type = "heatmap"), "Heatmap")
    expect_true(ggplot2::is.ggplot(plotBambu(seCombined, type = "pca")))
    expect_is(plotBambu(seCombined, type = "annotation", gene_id = unique(rowData(seCombined)$GENEID)[c(4, 6)]), "list")
    expect_is(plotBambu(seCombined, type = "annotation", transcript_id = rownames(seCombined)[c(4, 6)]), "grob")


    # case 2: group.variable is provided
    colData(seCombined)$groupVar <- c("group1", "group2")
    expect_is(plotBambu(seCombined, group.variable = "groupVar", type = "heatmap"), "Heatmap")
    expect_true(ggplot2::is.ggplot(plotBambu(seCombined, group.variable = "groupVar", type = "pca")))
})


test_that("visualization for gene expression  is successful", {
    colnames(seCombinedGeneExpected) <- colData(seCombinedGeneExpected)$name <- c("sample1", "sample2")

    set.seed(1)
    assays(seCombinedGeneExpected)$CPM[, 2] <- pmax(0, rnorm(length(assays(seCombinedGeneExpected)$CPM[, 2]), assays(seCombinedGeneExpected)$CPM[, 2], 10))

    # case 1: no group.variable provided
    expect_is(plotBambu(seCombinedGeneExpected, type = "heatmap"), "Heatmap")
    expect_true(ggplot2::is.ggplot(plotBambu(seCombinedGeneExpected, type = "pca")))
    expect_is(plotBambu(seCombinedGeneExpected, type = "annotation", gene_id = rownames(seCombinedGeneExpected)[c(4, 6)]), "grob")

    # case 2: group.variable is provided
    colData(seCombinedGeneExpected)$groupVar <- c("group1", "group2")
    expect_is(plotBambu(seCombinedGeneExpected, group.variable = "groupVar", type = "heatmap"), "Heatmap")
    expect_true(ggplot2::is.ggplot(plotBambu(seCombinedGeneExpected, group.variable = "groupVar", type = "pca")))
})
