context("Isoform quantification")

# bambu.quantDT
test_that("generic function of isoform quantification of data.table is list of 2",{
  # test case
  #   1: Overlapping scenario with no read support for each transcript
  #   2: Same scenario as in 1 but no empty class
  #   3: Transcripts not overlapping
  #   4: Second overlapping scenario with read support for each transcript
  #   5: A real example of observed counts with complex overlapping scenario

  lapply(1:5, function(s){
    est <- bambu.quantDT(readClassDt = get(paste0("data",s)))
    expect_type(est, "list")
    expect_equal(est, estOutput_woBC[[s]])
  })

  ## with bias correction
  lapply(1:5, function(s){
    est <- bambu.quantDT(readClassDt = get(paste0("data",s)))
    expect_type(est, "list")
    expect_equal(est, estOutput_wBC[[s]])
  })

})




test_that("bambu (isoform quantification of bam file) produces expected output",{

  test.bam <- system.file("extdata", "SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.bam", package = "bambu")
  fa.file <- system.file("extdata", "Homo_sapiens.GRCh38.dna_sm.primary_assembly_chr9_1_1000000.fa", package = "bambu")


  txdb <- loadDb(system.file("extdata", "Homo_sapiens.GRCh38.91.annotations-txdb_chr9_1_1000000.sqlite", package = "bambu"))
  gr <- readRDS(system.file("extdata", "annotationGranges_txdbGrch38_91_chr9_1_1000000.rds", package = "bambu"))


  seExpected <- readRDS(system.file("extdata", "seOutput_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.rds", package = "bambu"))
  seExtendedExpected <- readRDS(system.file("extdata", "seOutputExtended_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.rds", package = "bambu"))
  seCombinedExpected <- readRDS(system.file("extdata", "seOutputCombined_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.rds", package = "bambu"))
  seCombinedExtendedExpected <- readRDS(system.file("extdata", "seOutputCombinedExtended_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.rds", package = "bambu"))


  
  # test case 1: bambu with single bam file, only using annotations (default option)
  set.seed(1234)
  se = bambu(reads = test.bam,  annotations = txdb, genomeSequence = fa.file, emParameters=list(bias = FALSE), extendAnnotations = FALSE)
  expect_s4_class(se, "SummarizedExperiment")
  expect_equal(assays(se),assays(seExpected))

  set.seed(1234)
  se = bambu(reads = test.bam,  annotations =  gr, genomeSequence = fa.file, emParameters=list(bias = FALSE), extendAnnotations = FALSE)
  expect_s4_class(se, "SummarizedExperiment")
  expect_equal(assays(se),assays(seExpected))

  set.seed(1234)
  seExtended = bambu(reads = test.bam,  annotations =  gr, genomeSequence = fa.file, emParameters=list(bias = FALSE), extendAnnotations = TRUE)
  expect_s4_class(seExtended, "SummarizedExperiment")
  expect_equal(assays(seExtended),assays(seExtendedExpected))



  # test case 2: bambu with multiple bam file, only using annotations (default option), yieldSize lower than read count
  set.seed(1234)
  seCombined = bambu(reads =  Rsamtools::BamFileList(c(test.bam, test.bam), yieldSize = 1000),  annotations =  gr, genomeSequence = fa.file, extendAnnotations = FALSE)
  expect_s4_class(seCombined, "SummarizedExperiment")
  expect_equal(seCombined,seCombinedExpected)

  # test case 3: bambu with multiple bam file, extending annotations, yieldSize lower than read count
  set.seed(1234)
  seCombinedExtended = bambu(reads =  Rsamtools::BamFileList(c(test.bam, test.bam), yieldSize = 1000), annotations =  gr, genomeSequence = fa.file, extendAnnotations = TRUE)
  expect_s4_class(seCombinedExtended, "SummarizedExperiment")
  expect_equal(seCombinedExtended,seCombinedExtendedExpected)
  
})



test_that("bambu (isoform quantification of bam file and save readClassFiles) produces expected output",{
  ## ToDo: update data sets for comparison

  test.bam <- system.file("extdata", "SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.bam", package = "bambu")
  fa.file <- system.file("extdata", "Homo_sapiens.GRCh38.dna_sm.primary_assembly_chr9_1_1000000.fa", package = "bambu")

  txdb <- loadDb(system.file("extdata", "Homo_sapiens.GRCh38.91.annotations-txdb_chr9_1_1000000.sqlite", package = "bambu"))
  gr <- readRDS(system.file("extdata", "annotationGranges_txdbGrch38_91_chr9_1_1000000.rds", package = "bambu"))
  readClass.outputDir = tempdir()


  seExpected <- readRDS(system.file("extdata", "seOutput_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.rds", package = "bambu"))
  seExtendedExpected <- readRDS(system.file("extdata", "seOutputExtended_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.rds", package = "bambu"))
  seCombinedExpected <- readRDS(system.file("extdata", "seOutputCombined_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.rds", package = "bambu"))
  seCombinedExtendedExpected <- readRDS(system.file("extdata", "seOutputCombinedExtended_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.rds", package = "bambu"))

  # test case 1: bambu with single bam file, only using annotations (default option)
  set.seed(1234)
  se = bambu(reads = test.bam,  annotations = txdb, genomeSequence = fa.file, emParameters = list(bias = FALSE), extendAnnotations = FALSE, readClass.outputDir = readClass.outputDir)
  expect_s4_class(se, "SummarizedExperiment")
  expect_equal(se,seExpected)
  
  
  set.seed(1234)
  se = bambu(reads = test.bam,  annotations = gr, genomeSequence = fa.file, emParameters = list(bias = FALSE), extendAnnotations = FALSE, readClass.outputDir = readClass.outputDir)
  expect_s4_class(se, "SummarizedExperiment")
  expect_equal(se, seExpected)
  
  
  set.seed(1234)
  seExtended = bambu(reads = test.bam,  annotations = gr, genomeSequence = fa.file, emParameters = list(bias = FALSE), extendAnnotations = TRUE, readClass.outputDir = readClass.outputDir)
  expect_s4_class(seExtended, "SummarizedExperiment")
  expect_equal(seExtended,seExtendedExpected)
  
  
  
  # test case 2: bambu with multiple bam file, only using annotations (default option), yieldSize lower than read count
  set.seed(1234)
  seCombined = bambu(reads = Rsamtools::BamFileList(c(test.bam, test.bam), yieldSize = 1000),  annotations =  gr, genomeSequence = fa.file,  extendAnnotations = FALSE,readClass.outputDir = readClass.outputDir)
  expect_s4_class(seCombined, "SummarizedExperiment")
  expect_equal(seCombined,seCombinedExpected)
  
  
  
  # test case 3: bambu with multiple bam file, extending annotations, yieldSize lower than read count
  set.seed(1234)
  seCombinedExtended = bambu(reads = Rsamtools::BamFileList(c(test.bam, test.bam), yieldSize = 1000),  annotations =  gr, genomeSequence = fa.file, extendAnnotations = TRUE, readClass.outputDir = readClass.outputDir)
  expect_s4_class(seCombinedExtended, "SummarizedExperiment")
  expect_equal(seCombinedExtended, seCombinedExtendedExpected)
  
})




test_that("bambu (isoform quantification of saved readClassFiles) produces expected output",{
  ## ToDo: update data sets for comparison
  seReadClass1 <- system.file("extdata", "seReadClass_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.rds", package = "bambu")
  gr <- readRDS(system.file("extdata", "annotationGranges_txdbGrch38_91_chr9_1_1000000.rds", package = "bambu"))


  seExpected <- readRDS(system.file("extdata", "seOutput_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.rds", package = "bambu"))
  seExtendedExpected <- readRDS(system.file("extdata", "seOutputExtended_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.rds", package = "bambu"))
  seCombinedExpected <- readRDS(system.file("extdata", "seOutputCombined_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.rds", package = "bambu"))
  seCombinedExtendedExpected <- readRDS(system.file("extdata", "seOutputCombinedExtended_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.rds", package = "bambu"))


  # test case 1: bambu with single bam file, only using annotations (default option)
  set.seed(1234)
  se = bambu(readClass.file = seReadClass1,  annotations = gr, extendAnnotations = FALSE)
  expect_s4_class(se, "SummarizedExperiment")
  expect_equal(se,seExpected)

  set.seed(1234)
  seExtended = bambu(readClass.file = seReadClass1,  annotations = gr, extendAnnotations = TRUE)
  expect_s4_class(seExtended, "SummarizedExperiment")
  expect_equal(seExtended,seExtendedExpected)


  # test case 2: bambu with multiple bam file, only using annotations (default option), yieldSize lower than read count
  set.seed(1234)
  seCombined = bambu(readClass.file = c(seReadClass1, seReadClass1),  annotations =  gr, extendAnnotations = FALSE)
  expect_s4_class(seCombined, "SummarizedExperiment")
  expect_equal(seCombined,seCombinedExpected)


  # test case 3: bambu with multiple bam file, extending annotations, yieldSize lower than read count
  set.seed(1234)
  seCombinedExtended = bambu(readClass.file = c(seReadClass1, seReadClass1), annotations =  gr, extendAnnotations = TRUE)
  expect_s4_class(seCombinedExtended, "SummarizedExperiment")
  expect_equal(seCombinedExtended,seCombinedExtendedExpected)


})


