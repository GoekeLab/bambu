Sys.setenv("R_TESTS" = "")

library(testthat)
library(bambu)
library(BiocManager)
install("BSgenome.Hsapiens.NCBI.GRCh38")
install("TxDb.Hsapiens.UCSC.hg38.knownGene")

library(BSgenome.Hsapiens.NCBI.GRCh38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)



test_check("bambu")
