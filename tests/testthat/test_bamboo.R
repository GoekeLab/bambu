context("Isoform quantification")
library(bamboo)
library(testthat)
# bamboo.quantDT
test_that("generic function of isoform quantification of data.table is list of 2",{
  # test case
  #   1: Overlapping scenario with no read support for each transcript
  #   2: Same scenario as in 1 but no empty class
  #   3: Transcripts not overlapping
  #   4: Second overlapping scenario with read support for each transcript
  #   5: A real example of observed counts with complex overlapping scenario

  ## without bias correction
  lapply(1:5, function(s){
    est <- bamboo.quantDT(dt = get(paste0("data",s)),  algo.control=list(ncore = 1))
    expect_type(est, "list")
    expect_equal(est, estOutput_woBC[[s]])
  })

  ## with bias correction
  lapply(1:5, function(s){
    est <- bamboo.quantDT(dt = get(paste0("data",s)),  algo.control=list(ncore = 1))
    expect_type(est, "list")
    expect_equal(est, estOutput_wBC[[s]])
  })

})



# bamboo.quantISORE
test_that("bamboo.quantISORE (isoform quantification of bam file) produces expected output",{

  test.bam <- system.file("extdata", "SGNex_HepG2_cDNAStranded_replicate5_run4_chr9_108865774_109014097.bam", package = "bamboo")
  fa.file <- system.file("extdata", "Homo_sapiens.GRCh38.dna_sm.primary_assembly_chr9.fa.gz", package = "bamboo")


  samplefile <- "http://s3.ap-southeast-1.amazonaws.com/ucsc-trackdata.store.genome.sg/chenying/bamboo_exampleDataset/Homo_sapiens.GRCh38.91.annotations-txdb.sqlite"
  # download remote sqlite file as a temporary object
  dbfile <- tempfile(fileext=".sqlite")
  download.file(samplefile, dbfile)
  txdb <- AnnotationDbi::loadDb(dbfile)
  gr <- readRDS(system.file("extdata", "annotationGranges_txdbGrch38_91_chr9.rds", package = "bamboo"))

  # test case 1: bamboo with single bam file, only using annotations (default option)
  se = bamboo(reads = test.bam,  annotations = txdb, genomeSequence = fa.file, algo.control=list(bias_correction = FALSE, ncore = 1), extendAnnotations = FALSE)
  expect_type(se, "list")

  se = bamboo(reads = test.bam,  annotations =  gr, genomeSequence = fa.file, algo.control=list(bias_correction = FALSE, ncore = 1), extendAnnotations = FALSE)
  expect_type(se, "list")
  # expect_equal(se, expectedSE)

  se = bamboo(reads = test.bam,  annotations =  gr, genomeSequence = fa.file, algo.control=list(bias_correction = FALSE, ncore = 1), extendAnnotations = TRUE)
  expect_type(se, "list")
  # expect_equal(se, expectedSE_extended)


  # test case 2: bamboo with multiple bam file, only using annotations (default option), yieldSize lower than read count
  se = bamboo(reads = BamFileList(c(test.bam, test.bam), yieldSize = 1000),  annotations =  gr, genomeSequence = fa.file, algo.control = list(ncore = 1), extendAnnotations = FALSE)
  expect_type(se, "list")
  # expect_equal(assays(se[["Transcript"]]),assays(expectedSE[["Transcript"]]))
  # expect_equal(assays(se[["Gene"]])[,1],assays(expectedSE[["Gene"]]))

  # test case 3: bamboo with multiple bam file, extending annotations, yieldSize lower than read count
  se = bamboo(reads = BamFileList(c(test.bam, test.bam), yieldSize = 1000),  annotations =  gr, genomeSequence = fa.file, algo.control = list(ncore = 1), extendAnnotations = TRUE)
  expect_type(se, "list")
  # expect_equal(assays(se[["Transcript"]])[,1],assays(expectedSE_extended[["Transcript"]]))
  # expect_equal(assays(se[["Gene"]])[,1],assays(expectedSE_extended[["Gene"]]))


})



test_that("bamboo.preprocess (isoform quantification of bam file and save readClassFiles) produces expected output",{
  ## ToDo: update data sets for comparison

  test.bam <- system.file("extdata", "SGNex_HepG2_cDNAStranded_replicate5_run4_chr9_108865774_109014097.bam", package = "bamboo")
  fa.file <- system.file("extdata", "Homo_sapiens.GRCh38.dna_sm.primary_assembly_chr9.fa.gz", package = "bamboo")

  samplefile <- "http://s3.ap-southeast-1.amazonaws.com/ucsc-trackdata.store.genome.sg/chenying/bamboo_exampleDataset/Homo_sapiens.GRCh38.91.annotations-txdb.sqlite"
  # download remote sqlite file as a temporary object
  dbfile <- tempfile(fileext=".sqlite")
  download.file(samplefile, dbfile)
  txdb <- AnnotationDbi::loadDb(dbfile)
  gr <- readRDS(system.file("extdata", "annotationGranges_txdbGrch38_91_chr9.rds", package = "bamboo"))
  stranded=FALSE
  outputReadClassDir = tempdir()

  # test case 1: bamboo with single bam file, only using annotations (default option)
  se = bamboo(reads = test.bam,  annotations = txdb, genomeSequence = fa.file, algo.control = list(ncore = 1), extendAnnotations = FALSE, outputReadClassDir = outputReadClassDir)
  expect_type(se, "list")


  se = bamboo(reads = test.bam,  annotations = gr, genomeSequence = fa.file, algo.control = list(ncore = 1), extendAnnotations = FALSE, outputReadClassDir = outputReadClassDir)
  expect_type(se, "list")
  # expect_equal(se, expectedSE)

  se = bamboo(reads = test.bam,  annotations = gr, genomeSequence = fa.file, algo.control = list(ncore = 1), extendAnnotations = TRUE, outputReadClassDir = outputReadClassDir)
  expect_type(se, "list")
  # expect_equal(se, expectedSE_extended)



  # test case 2: bamboo with multiple bam file, only using annotations (default option), yieldSize lower than read count
  se = bamboo(reads = BamFileList(c(test.bam, test.bam), yieldSize = 1000),  annotations =  gr, genomeSequence = fa.file, algo.control = list(ncore = 1), extendAnnotations = FALSE,outputReadClassDir = outputReadClassDir)
  expect_type(se, "list")
  # expect_equal(assays(se[["Transcript"]])[,1],assays(expectedSE[["Transcript"]]))
  # expect_equal(assays(se[["Gene"]])[,1],assays(expectedSE[["Gene"]]))



  # test case 3: bamboo with multiple bam file, extending annotations, yieldSize lower than read count
  se = bamboo(reads = BamFileList(c(test.bam, test.bam), yieldSize = 1000),  annotations =  gr, genomeSequence = fa.file, algo.control = list(ncore = 1), extendAnnotations = TRUE, outputReadClassDir = outputReadClassDir)
  expect_type(se, "list")
  # expect_equal(assays(se[["Transcript"]])[,1],assays(expectedSE_extended[["Transcript"]]))
  # expect_equal(assays(se[["Gene"]])[,1],assays(expectedSE_extended[["Gene"]]))


})




test_that("bamboo.combineQuantify (isoform quantification of saved readClassFiles) produces expected output",{
  ## ToDo: update data sets for comparison
  readclass.file <- system.file("extdata", "SGNex_HepG2_cDNAStranded_replicate5_run4_chr9_108865774_109014097_readClassSe.rds", package = "bamboo")
  gr <- readRDS(system.file("extdata", "annotationGranges_txdbGrch38_91_chr9.rds", package = "bamboo"))


  # test case 1: bamboo with single bam file, only using annotations (default option)
  se = bamboo(readclass.file = readclass.file,  annotations = gr, algo.control = list(ncore = 1), extendAnnotations = FALSE)
  expect_type(se, "list")
  #expect_equal(se, expectedSE)


  se = bamboo(readclass.file = readclass.file,  annotations = gr, algo.control = list(ncore = 1), extendAnnotations = TRUE)
  expect_type(se, "list")
  #expect_equal(se, expectedSE_extended)


  # test case 2: bamboo with multiple bam file, only using annotations (default option), yieldSize lower than read count
  se = bamboo(readclass.file = c(readclass.file, readclass.file),  annotations =  gr, algo.control = list(ncore = 1), extendAnnotations = FALSE)
  expect_type(se, "list")
  #expect_equal(assays(se[["Transcript"]])[,1],assays(expectedSE[["Transcript"]]))
  #expect_equal(assays(se[["Gene"]])[,1],assays(expectedSE[["Gene"]]))



  # test case 3: bamboo with multiple bam file, extending annotations, yieldSize lower than read count
  se = bamboo(readclass.file = c(readclass.file, readclass.file), annotations =  gr, algo.control = list(ncore = 1), extendAnnotations = TRUE)
  expect_type(se, "list")
  #expect_equal(assays(se[["Transcript"]])[,1],assays(expectedSE_extended[["Transcript"]]))
  #expect_equal(assays(se[["Gene"]])[,1],assays(expectedSE_extended[["Gene"]]))


})


