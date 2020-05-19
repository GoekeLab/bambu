context("Visualization of estimates")

test_that("visualization of heatmap is successful",{
  seCombined <- readRDS(system.file("extdata", "seOutputCombined_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.rds", package = "bambu"))

  colnames(seCombined) <- colData(seCombined)$name <- c("sample1","sample2")
  set.seed(1)
  assays(seCombined)$CPM[,2]  <- pmax(0, rnorm(length(assays(seCombined)$CPM[,2]),assays(seCombined)$CPM[,2],10))

  # case 1: no group.variable provided
  expect_is(plot.bambu(seCombined, type = "heatmap"),"Heatmap")
  expect_true(is.ggplot(plot.bambu(seCombined, type = "pca")))

  # case 2: group.variable is provided
  colData(seCombined)$groupVar <- c("group1","group2")
  expect_is(plot.bambu(seCombined, group.variable = "groupVar", type = "heatmap"),"Heatmap")
  expect_true(is.ggplot(plot.bambu(seCombined, group.variable = "groupVar", type = "pca")))

})



