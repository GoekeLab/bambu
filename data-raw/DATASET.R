## code to prepare `example_data` dataset goes here
library(data.table)

example_data <- list()

a <- 500

example_data[[1]] <-  data.table(tx_id = c(1,2,3,1,2,1,2,3),
                                 read_class_id = c(1,2,3,rep(4,2),rep(5,3)),
                                 nobs = c(0,0,0,rep(10,2),rep(a,3)),
                                 gene_id = 1)

example_data[[2]] <- data.table(tx_id = c(1,2,1,2,3),
                                read_class_id = c(rep(4,2),rep(5,3)),
                                nobs = c(rep(10,2),rep(a,3)),
                                gene_id = 2)

example_data[[3]] <- data.table(tx_id = c(1,2),
                                read_class_id = c(1,2),
                                nobs = c(5,10),
                                gene_id = 3)

example_data[[4]] <- data.table(tx_id = c(1,2,3,1,2,1,2,3),
                                read_class_id = c(1,2,3,rep(4,2),rep(5,3)),
                                nobs = c(100,5,500,rep(20,2),rep(5,3)),
                                gene_id = 4)

example_data[[5]] <- data.table(tx_id = c(1:6,2,3,1,2,5,1,2,3,5,1:5),
                                read_class_id = c(1:6,rep(7,2),rep(8,3),rep(9,4),rep(10,5)),
                                nobs = c(2,0,1,0,0,1,rep(2,2),rep(282,3),rep(64,4),rep(5,5)),
                                gene_id = 5)

samplefile <- system.file("extdata","Homo_sapiens.GRCh38.91.annotations-txdb.sqlite", package = "bamboo")
txdb <- AnnotationDbi::loadDb(samplefile)
example_txdbTablesList <- prepareAnnotations(txdb)

example_se <- isore(bamFile = test.bam,  txdbTablesList = example_txdbTablesList,genomeFA = fa.file)

test.bam <- system.file("extdata", "GIS_HepG2_cDNAStranded_Rep5-Run4_chr9_108865774_109014097.bam", package = "bamboo")
fa.file <- system.file("extdata", "Homo_sapiens.GRCh38.dna_sm.primary_assembly_chr9.fa.gz", package = "bamboo")
example_quantOutput <- bamboo(test.bam, fa.file = fa.file, txdbTablesList = example_txdbTablesList,algo.control=list(ncore = 1))

usethis::use_data("example_data","example_txdbTablesList","example_se","example_quantOutput")
