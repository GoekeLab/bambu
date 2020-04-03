context("Isoform quantification")
library(bamboo)


test_that("generic function of isoform quantification of data.table is list of 2",{
  expect_type(bamboo.quantDT(dt = example_data[[1]],  algo.control=list(ncore = 1)), "list")
  expect_output(bamboo.quantDT(dt = example_data[[1]], algo.control=list(ncore = 1)),"Finished EM estimation in")
  expect_type(bamboo.quantDT(dt = example_data[[2]], algo.control=list(ncore = 1)), "list")
  expect_output(bamboo.quantDT(dt = example_data[[2]],  algo.control=list(ncore = 1)),"Finished EM estimation in")
  expect_type(bamboo.quantDT(dt = example_data[[3]], algo.control=list(ncore = 1)), "list")
  expect_output(bamboo.quantDT(dt = example_data[[3]], algo.control=list(ncore = 1)),"Finished EM estimation in")
  expect_type(bamboo.quantDT(dt = example_data[[4]], algo.control=list(ncore = 1)), "list")
  expect_output(bamboo.quantDT(dt = example_data[[4]], algo.control=list(ncore = 1)),"Finished EM estimation in")
  expect_type(bamboo.quantDT(dt = example_data[[5]], algo.control=list(ncore = 1)), "list")
  expect_output(bamboo.quantDT(dt = example_data[[5]], algo.control=list(ncore = 1)),"Finished EM estimation in")
})


test_that("isoform quantification of summarizedExperiment is summrizedExperiment",{
  example_txdbTablesList <- readRDS(url("http://s3.ap-southeast-1.amazonaws.com/ucsc-trackdata.store.genome.sg/chenying/bamboo_exampleDataset/HomoSapiens_ensembl91_Grch38_txdbTablesList.rds"))
  example_quantOutput <- readRDS(url("http://s3.ap-southeast-1.amazonaws.com/ucsc-trackdata.store.genome.sg/chenying/bamboo_exampleDataset/example_quantOutput.rds"))

  expect_s4_class(bamboo(se = example_se, txdbTablesList = example_txdbTablesList, algo.control=list(ncore = 1)),"SummarizedExperiment")
  expect_output(bamboo(se = example_se,txdbTablesList = example_txdbTablesList,algo.control=list(ncore = 1)),"Finished EM estimation in")
  expect_equal(bamboo(se = example_se,txdbTablesList = example_txdbTablesList,algo.control=list(ncore = 1)),example_quantOutput)
})

test_that("bamboo.quantISORE (isoform quantification of bam file) produces expected output",{
   ## ToDo: update data sets for comparison
  test.bam <- system.file("extdata", "GIS_HepG2_cDNAStranded_Rep5-Run4_chr9_108865774_109014097.bam", package = "bamboo")
  fa.file <- system.file("extdata", "Homo_sapiens.GRCh38.dna_sm.primary_assembly_chr9.fa.gz", package = "bamboo")
  gr <- prepareAnnotations(txdb) ## saveRDS(prepareAnnotations(txdb), file='annotationGranges.rds', compress = 'xz')
  stranded=FALSE

  # needs to be updated
  #example_quantOutput <- readRDS(url("http://s3.ap-southeast-1.amazonaws.com/ucsc-trackdata.store.genome.sg/chenying/bamboo_exampleDataset/example_quantOutput.rds")) ## ToDo: should be updated

  # test case 1: bamboo with single bam file, only using annotations (default option)
  se = bamboo(bam.file = test.bam,  txdb = txdb, fa.file = fa.file, algo.control=list(ncore = 1), extendAnnotations=F)
  se = bamboo(bam.file = test.bam,  annotationGrangesList =  gr, fa.file = fa.file, algo.control=list(ncore = 1), extendAnnotations=F)
  se = bamboo(bam.file = test.bam,  annotationGrangesList =  gr, fa.file = fa.file, algo.control=list(ncore = 1), extendAnnotations=T)


  expect_s4_class(se, "SummarizedExperiment")
 # expect_equal(assays(se),assays(example_quantOutput))
 # expect_equal(se,se.isore)

  # test case 2: bamboo with multiple bam file, only using annotations (default option), yieldSize lower than read count
  se = bamboo(bam.file = BamFileList(c(test.bam, test.bam), yieldSize = 1000),  annotationGrangesList =  gr, fa.file = fa.file, algo.control=list(ncore = 1), extendAnnotations=FALSE)
  expect_s4_class(se, "SummarizedExperiment")
#  expect_equal(assays(se)[,1],assays(example_quantOutput))

  # test case 3: bamboo with multiple bam file, extending annotations, yieldSize lower than read count
  se = bamboo(bam.file = BamFileList(c(test.bam, test.bam), yieldSize = 1000),  annotationGrangesList =  gr, fa.file = fa.file, algo.control=list(ncore = 1), extendAnnotations=TRUE)
  expect_s4_class(se, "SummarizedExperiment")
 # expect_equal(assays(se)[,1],assays(example_quantOutput))


})





test_that("bamboo.preprocess (isoform quantification of bam file and save readClassFiles) produces expected output",{
  ## ToDo: update data sets for comparison
  samplefile <- "http://s3.ap-southeast-1.amazonaws.com/ucsc-trackdata.store.genome.sg/chenying/bamboo_exampleDataset/Homo_sapiens.GRCh38.91.annotations-txdb.sqlite"
  # download remote sqlite file as a temporary object
  dbfile <- tempfile(fileext=".sqlite")
  download.file(samplefile, dbfile)
  txdb <- AnnotationDbi::loadDb(dbfile)
  test.bam <- system.file("extdata", "GIS_HepG2_cDNAStranded_Rep5-Run4_chr9_108865774_109014097.bam", package = "bamboo")
  fa.file <- system.file("extdata", "Homo_sapiens.GRCh38.dna_sm.primary_assembly_chr9.fa.gz", package = "bamboo")
  gr <- prepareAnnotations(txdb) ## saveRDS(prepareAnnotations(txdb), file='annotationGranges.rds', compress = 'xz')
  stranded=FALSE
  outputReadClassDir = tempdir()
  # needs to be updated
  #example_quantOutput <- readRDS(url("http://s3.ap-southeast-1.amazonaws.com/ucsc-trackdata.store.genome.sg/chenying/bamboo_exampleDataset/example_quantOutput.rds")) ## ToDo: should be updated

  # test case 1: bamboo with single bam file, only using annotations (default option)
  se = bamboo(bam.file = test.bam,  txdb = txdb, fa.file = fa.file, algo.control=list(ncore = 1), extendAnnotations=F, outputReadClassDir = outputReadClassDir )
  se = bamboo(bam.file = test.bam,  annotationGrangesList =  gr, fa.file = fa.file, algo.control=list(ncore = 1), extendAnnotations=F, outputReadClassDir = outputReadClassDir)
  se = bamboo(bam.file = test.bam,  annotationGrangesList =  gr, fa.file = fa.file, algo.control=list(ncore = 1), extendAnnotations=T, outputReadClassDir = outputReadClassDir)


  expect_s4_class(se, "SummarizedExperiment")
  # expect_equal(assays(se),assays(example_quantOutput))
  # expect_equal(se,se.isore)

  # test case 2: bamboo with multiple bam file, only using annotations (default option), yieldSize lower than read count
  se = bamboo(bam.file = BamFileList(c(test.bam, test.bam), yieldSize = 1000),  annotationGrangesList =  gr, fa.file = fa.file, algo.control=list(ncore = 1), extendAnnotations=FALSE,outputReadClassDir = outputReadClassDir)
  expect_s4_class(se, "SummarizedExperiment")
  #  expect_equal(assays(se)[,1],assays(example_quantOutput))

  # test case 3: bamboo with multiple bam file, extending annotations, yieldSize lower than read count
  se = bamboo(bam.file = BamFileList(c(test.bam, test.bam), yieldSize = 1000),  annotationGrangesList =  gr, fa.file = fa.file, algo.control=list(ncore = 1), extendAnnotations=TRUE, outputReadClassDir = outputReadClassDir)
  expect_s4_class(se, "SummarizedExperiment")
  # expect_equal(assays(se)[,1],assays(example_quantOutput))


})


