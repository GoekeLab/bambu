context("Isoform quantification")
library(bambu)
library(testthat)


# bambu.quantDT
test_that("generic function of isoform quantification of data.table is list of 2",{
  # test case
  #   1: Overlapping scenario with no read support for each transcript
  #   2: Same scenario as in 1 but no empty class
  #   3: Transcripts not overlapping
  #   4: Second overlapping scenario with read support for each transcript
  #   5: A real example of observed counts with complex overlapping scenario

  ## without bias correction
  lapply(1:5, function(s){
    est <- bambu.quantDT(dt = get(paste0("data",s)),  algo.control=list(ncore = 1))
    expect_type(est, "list")
    expect_equal(est, estOutput_woBC[[s]])
  })

  ## with bias correction
  lapply(1:5, function(s){
    est <- bambu.quantDT(dt = get(paste0("data",s)),  algo.control=list(ncore = 1))
    expect_type(est, "list")
    expect_equal(est, estOutput_wBC[[s]])
  })

})



# bambu.quantISORE
test_that("bambu.quantISORE (isoform quantification of bam file) produces expected output",{

  test.bam <- system.file("extdata", "SGNex_HepG2_cDNAStranded_replicate5_run4_chr9_108865774_109014097.bam", package = "bambu")
  fa.file <- system.file("extdata", "Homo_sapiens.GRCh38.dna_sm.primary_assembly_chr9.fa.gz", package = "bambu")


  txdb <- loadDb(system.file("extdata", "Homo_sapiens.GRCh38.91.annotations-txdb_chr9_108865774_109014097.sqlite", package = "bambu"))

  gr <- readRDS(system.file("extdata", "annotationGranges_txdbGrch38_91_chr9_108865774_109014097.rds", package = "bambu"))

  # test case 1: bambu with single bam file, only using annotations (default option)
  set.seed(1234)
  se = bambu(reads = test.bam,  annotations = txdb, genomeSequence = fa.file, algo.control=list(bias_correction = FALSE, ncore = 1), extendAnnotations = FALSE)
  expect_type(se, "list")
  expect_equal(assays(se[["Transcript"]][,1]),assays(expectedSE[["Transcript"]]))
  expect_equal(assays(se[["Gene"]][,1]),assays(expectedSE[["Gene"]]))

  set.seed(1234)
  se = bambu(reads = test.bam,  annotations =  gr, genomeSequence = fa.file, algo.control=list(bias_correction = FALSE, ncore = 1), extendAnnotations = FALSE)
  expect_type(se, "list")
  expect_equal(assays(se[["Transcript"]][,1]),assays(expectedSE[["Transcript"]]))
  expect_equal(assays(se[["Gene"]][,1]),assays(expectedSE[["Gene"]]))

  set.seed(1234)
  se = bambu(reads = test.bam,  annotations =  gr, genomeSequence = fa.file, algo.control=list(bias_correction = FALSE, ncore = 1), extendAnnotations = TRUE)
  expect_type(se, "list")
  expect_equal(assays(se[["Transcript"]][,1]),assays(expectedSE_extended[["Transcript"]]))
  expect_equal(assays(se[["Gene"]][,1]),assays(expectedSE_extended[["Gene"]]))



  # test case 2: bambu with multiple bam file, only using annotations (default option), yieldSize lower than read count
  set.seed(1234)
  se = bambu(reads =  Rsamtools::BamFileList(c(test.bam, test.bam), yieldSize = 1000),  annotations =  gr, genomeSequence = fa.file, algo.control = list(ncore = 1), extendAnnotations = FALSE)
  expect_type(se, "list")
  expect_equal(assays(se[["Transcript"]][,1]),assays(expectedSE[["Transcript"]]))
  expect_equal(assays(se[["Gene"]][,1]),assays(expectedSE[["Gene"]]))

  # test case 3: bambu with multiple bam file, extending annotations, yieldSize lower than read count
  set.seed(1234)
  se = bambu(reads =  Rsamtools::BamFileList(c(test.bam, test.bam), yieldSize = 1000), annotations =  gr, genomeSequence = fa.file, algo.control = list(ncore = 1), extendAnnotations = TRUE)
  expect_type(se, "list")
  expect_equal(assays(se[["Transcript"]][,1]),assays(expectedSE_extended[["Transcript"]]))
  expect_equal(assays(se[["Gene"]][,1]),assays(expectedSE_extended[["Gene"]]))


})



test_that("bambu.preprocess (isoform quantification of bam file and save readClassFiles) produces expected output",{
  ## ToDo: update data sets for comparison

  test.bam <- system.file("extdata", "SGNex_HepG2_cDNAStranded_replicate5_run4_chr9_108865774_109014097.bam", package = "bambu")
  fa.file <- system.file("extdata", "Homo_sapiens.GRCh38.dna_sm.primary_assembly_chr9.fa.gz", package = "bambu")

  txdb <- loadDb(system.file("extdata", "Homo_sapiens.GRCh38.91.annotations-txdb_chr9_108865774_109014097.sqlite", package = "bambu"))
  gr <- readRDS(system.file("extdata", "annotationGranges_txdbGrch38_91_chr9_108865774_109014097.rds", package = "bambu"))
  outputReadClassDir = tempdir()

  # test case 1: bambu with single bam file, only using annotations (default option)
  set.seed(1234)
  se = bambu(reads = test.bam,  annotations = txdb, genomeSequence = fa.file, algo.control = list(bias_correction = FALSE, ncore = 1), extendAnnotations = FALSE, outputReadClassDir = outputReadClassDir)
  expect_type(se, "list")
  expect_equal(assays(se[["Transcript"]][,1]),assays(expectedSE[["Transcript"]]))
  expect_equal(assays(se[["Gene"]][,1]),assays(expectedSE[["Gene"]]))


  set.seed(1234)
  se = bambu(reads = test.bam,  annotations = gr, genomeSequence = fa.file, algo.control = list(bias_correction = FALSE, ncore = 1), extendAnnotations = FALSE, outputReadClassDir = outputReadClassDir)
  expect_type(se, "list")
  expect_equal(assays(se[["Transcript"]][,1]),assays(expectedSE[["Transcript"]]))
  expect_equal(assays(se[["Gene"]][,1]),assays(expectedSE[["Gene"]]))


  set.seed(1234)
  se = bambu(reads = test.bam,  annotations = gr, genomeSequence = fa.file, algo.control = list(bias_correction = FALSE, ncore = 1), extendAnnotations = TRUE, outputReadClassDir = outputReadClassDir)
  expect_type(se, "list")
  expect_equal(assays(se[["Transcript"]][,1]),assays(expectedSE_extended[["Transcript"]]))
  expect_equal(assays(se[["Gene"]][,1]),assays(expectedSE_extended[["Gene"]]))



  # test case 2: bambu with multiple bam file, only using annotations (default option), yieldSize lower than read count
  set.seed(1234)
  se = bambu(reads = Rsamtools::BamFileList(c(test.bam, test.bam), yieldSize = 1000),  annotations =  gr, genomeSequence = fa.file, algo.control = list(bias_correction = FALSE, ncore = 1), extendAnnotations = FALSE,outputReadClassDir = outputReadClassDir)
  expect_type(se, "list")
  expect_equal(assays(se[["Transcript"]][,1]),assays(expectedSE[["Transcript"]]))
  expect_equal(assays(se[["Gene"]][,1]),assays(expectedSE[["Gene"]]))



  # test case 3: bambu with multiple bam file, extending annotations, yieldSize lower than read count
  set.seed(1234)
  se = bambu(reads = Rsamtools::BamFileList(c(test.bam, test.bam), yieldSize = 1000),  annotations =  gr, genomeSequence = fa.file, algo.control = list(bias_correction = FALSE, ncore = 1), extendAnnotations = TRUE, outputReadClassDir = outputReadClassDir)
  expect_type(se, "list")
  expect_equal(assays(se[["Transcript"]][,1]),assays(expectedSE_extended[["Transcript"]]))
  expect_equal(assays(se[["Gene"]][,1]),assays(expectedSE_extended[["Gene"]]))


})




test_that("bambu.combineQuantify (isoform quantification of saved readClassFiles) produces expected output",{
  ## ToDo: update data sets for comparison
  readclass.file <- system.file("extdata", "SGNex_HepG2_cDNAStranded_replicate5_run4_chr9_108865774_109014097_readClassSe.rds", package = "bambu")
  gr <- readRDS(system.file("extdata", "annotationGranges_txdbGrch38_91_chr9_108865774_109014097.rds", package = "bambu"))


  # test case 1: bambu with single bam file, only using annotations (default option)
  set.seed(1234)
  se = bambu(readclass.file = readclass.file,  annotations = gr, algo.control = list(ncore = 1), extendAnnotations = FALSE)
  expect_type(se, "list")
  expect_equal(assays(se[["Transcript"]][,1]),assays(expectedSE[["Transcript"]]))
  expect_equal(assays(se[["Gene"]][,1]),assays(expectedSE[["Gene"]]))

  set.seed(1234)
  se = bambu(readclass.file = readclass.file,  annotations = gr, algo.control = list(ncore = 1), extendAnnotations = TRUE)
  expect_type(se, "list")
  expect_equal(assays(se[["Transcript"]][,1]),assays(expectedSE_extended[["Transcript"]]))
  expect_equal(assays(se[["Gene"]][,1]),assays(expectedSE_extended[["Gene"]]))


  # test case 2: bambu with multiple bam file, only using annotations (default option), yieldSize lower than read count
  set.seed(1234)
  se = bambu(readclass.file = c(readclass.file, readclass.file),  annotations =  gr, algo.control = list(ncore = 1), extendAnnotations = FALSE)
  expect_type(se, "list")
  expect_equal(assays(se[["Transcript"]][,1]),assays(expectedSE[["Transcript"]]))
  expect_equal(assays(se[["Gene"]][,1]),assays(expectedSE[["Gene"]]))



  # test case 3: bambu with multiple bam file, extending annotations, yieldSize lower than read count
  set.seed(1234)
  se = bambu(readclass.file = c(readclass.file, readclass.file), annotations =  gr, algo.control = list(ncore = 1), extendAnnotations = TRUE)
  expect_type(se, "list")
  expect_equal(assays(se[["Transcript"]][,1]),assays(expectedSE_extended[["Transcript"]]))
  expect_equal(assays(se[["Gene"]][,1]),assays(expectedSE_extended[["Gene"]]))


})


