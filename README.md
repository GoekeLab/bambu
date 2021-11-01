<img src="figures/transparent-bambu.png" title="Bambu" alt="Bambu">

# bambu: reference-guided transcript discovery and quantification for long read RNA-Seq data

[![Lifecycle: maturing](https://img.shields.io/badge/lifecycle-stable-green.svg)](https://www.tidyverse.org/lifecycle/#stable)
[![Maintained?](https://img.shields.io/badge/Maintained%3F-Yes-brightgreen)](https://github.com/GoekeLab/bambu/graphs/contributors)
[![Install](https://img.shields.io/badge/Install-Github-brightgreen)](#installation)
[![R build status](https://github.com/goekelab/bambu/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/goekelab/bambu/actions)  
[![GitHub issues](https://img.shields.io/github/issues/goekelab/bambu)](https://github.com/goekelab/bambu/issues) 
[![GitHub pulls](https://img.shields.io/github/issues-pr/goekelab/bambu)](https://github.com/goekelab/bambu/pulls) 
[![BioC status](http://bioconductor.org/shields/build/release/bioc/bambu.svg)](http://bioconductor.org/checkResults/release/bioc-LATEST/bambu/)
[![BioC dev status](http://www.bioconductor.org/shields/build/devel/bioc/bambu.svg)](https://bioconductor.org/checkResults/devel/bioc-LATEST/bambu)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![CodeFactor](https://www.codefactor.io/repository/github/goekelab/bambu/badge)](https://www.codefactor.io/repository/github/goekelab/bambu)
[![codecov](https://codecov.io/gh/GoekeLab/bambu/branch/master/graph/badge.svg?token=PMeRi0r1tj)](https://codecov.io/gh/GoekeLab/bambu)


***bambu*** is a R package for multi-sample transcript discovery and quantification using long read RNA-Seq data. You can use ***bambu*** after read alignment to obtain expression estimates for known and novel transcripts and genes. The output from ***bambu*** can directly be used for visualisation and downstream analysis such as differential gene expression or transcript usage.


---

### Content

  - [Installation](#installation)
  - [General usage](#general-usage)
  - [Use precalculated annotation objects](#use-precalculated-annotation-objects)
  - [Advanced options](#advanced-options)
  - [Details on the output](#details-on-the-output)
  - [Complementary functions](#complementary-functions)
  - [Release History](#release-history)
  - [Citation](#citation)
  - [Contributors](#contributors)


### Installation

You can install ***bambu*** from bioconductor:
```rscript
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("bambu")
```

---

### General Usage 

The default mode to run ***bambu** is using a set of aligned reads (bam files), reference genome annotations (gtf file, TxDb object, or bambuAnnotation object), and reference genome sequence (fasta file or BSgenome). ***bambu*** will return a summarizedExperiment object with the genomic coordinates for annotated and new transcripts and transcript expression estimates.

We highly recommend to use the same annotations that were used for genome alignment. If you have a gtf file and fasta file you can run ***bambu*** with the following options:

```rscript
test.bam <- system.file("extdata", "SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.bam", package = "bambu")
  
fa.file <- system.file("extdata", "Homo_sapiens.GRCh38.dna_sm.primary_assembly_chr9_1_1000000.fa", package = "bambu")

gtf.file <- system.file("extdata", "Homo_sapiens.GRCh38.91_chr9_1_1000000.gtf", package = "bambu")

bambuAnnotations <- prepareAnnotations(gtf.file)

se <- bambu(reads = test.bam, annotations = bambuAnnotations, genome = fa.file)

```
**Transcript discovery only (no quantification)**

```rscript
bambu(reads = test.bam, annotations = txdb, genome = fa.file, quant = FALSE)
```

**Quantification of annotated transcripts and genes only (no transcript/gene discovery)**

```rscript
bambu(reads = test.bam, annotations = txdb, genome = fa.file, discovery = FALSE)
```

**Large sample number/ limited memory**     
For larger sample numbers we recommend to write the processed data to a file:

```rscript
bambu(reads = test.bam, rcOutDir = "./bambu/", annotations = bambuAnnotations, genome = fa.file)
```

For very large samples (>100 million reads) where memory is limiting we recommend running Bambu in lowMemory mode:

```rscript
bambu(reads = test.bam, annotations = bambuAnnotations, genome = fa.file, lowMemory = TRUE)
```

---


### Use precalculated annotation objects

You can also use precalculated annotations.

If you plan to run ***bambu*** more frequently, we recommend to save the bambuAnnotations object.

The bambuAnnotation object can be calculated from a *.gtf* file:

```rscript
annotations <- prepareAnnotation(gtf.file)
```

From *TxDb* object

```rscript
annotations <- prepareAnnotations(txdb)
```

---

### Advanced Options

**More stringent filtering thresholds imposed on potential novel transcripts**    
 
- Keep novel transcripts with min 5 read count in at least 1 sample: 

```rscript
bambu(reads, annotations, genome, opt.discovery = list(min.readCount = 5))
```

- Keep novel transcripts with min 5 samples having at least 2 counts:

```rscript
bambu(reads, annotations, genome, opt.discovery = list(min.sampleNumber = 5))
```

- Filter out transcripts with relative abundance within gene lower than 10%: 

```rscript
bambu(reads, annotations, genome, opt.discovery = list(min.readFractionByGene = 0.1))
```

- Set novel transcript discovery rate to 50% of the detected transcripts (lower is more): 

```rscript
bambu(reads, annotations, genome, NDR = 0.5)
```

**Quantification without bias correction**     

 The default estimation automatically does bias correction for expression estimates. However, you can choose to perform the quantification without bias correction.

```rscript
bambu(reads, annotations, genome, opt.em = list(bias = FALSE))
```

**Parallel computation**      
 ***bambu***  allows parallel computation.  

```rscript
bambu(reads, annotations, genome, ncore = 8)
```

See [our page](https://goekelab.github.io/bambu/) for a complete step-by-step workflow and manual on how to customize other condictions.

---

### Details on the output 

***bambu*** will output different results depending on whether *quant* mode is on. 

By default, *quant* is set to TRUE, 
so ***bambu*** will generate a *SummarizedExperiment* object that contains the transcript expression estimates.  

* access transcript expression estimates by ***counts()***, including a list of variables: counts, CPM, fullLengthCount, partialLengthCounts, and uniqueCounts, and theta
    + counts: expression estimates
    + CPM: sequencing depth normalized estimates
    + fullLengthCounts: estimates of read counts mapped as full length reads for each transcript
    + partialLengthCounts: estimates of read counts mapped as partial length reads for each transcript
    + uniqueCounts: counts of reads that are uniquely mapped to each transcript
    + theta: raw estimates
* access annotations that are matched to the transcript expression estimates by ***rowRanges()***
* access transcript to gene id map by ***rowData()***, *eqClass* that defines the equivalent class transcripts is also reported

In the case when *quant* is set to FALSE, i.e., only transcript discovery is performed, 
***bambu*** will report the *grangeslist* of the extended annotations

### Complementary functions

**Transcript expression to gene expression**

```rscript
transcriptToGeneExpression(se)
```

**Visualization**

 You can visualize the novel genes/transcripts using ***plotBambu*** function 

```rscript
plotBambu(se, type = "annotation", gene_id)

plotBambu(se, type = "annotation", transcript_id)
```

- ***plotBambu*** can also be used to visualize the clustering of input samples on gene/transcript expressions

```rscript
plotBambu(se, type = "heatmap") # heatmap 

plotBambu(se, type = "pca") # PCA visualization
```

- ***plotBambu*** can also be used to visualize the clustering of input samples on gene/transcript expressions with grouping variable

```rscript
plotBambu(se, type = "heatmap", group.var) # heatmap 

plotBambu(se, type = "pca", group.var) # PCA visualization
```

**Write bambu outputs to files**

- ***writeBambuOutput*** will generate three files, including a *.gtf* file for the extended annotations, and two *.txt* files for the expression counts at transcript and gene levels.

```rscript
writeBambuOutput(se, path = "./bambu/")
```
---

### Release History

**bambu version 1.99.0**

Release date: 2021-10-18

Major Changes:

- Implemented a machine learning model to estimate transcript-level novel discovery rate
- Implemented full length estimates, partial length estimates and unique read counts in final output
- Improved the performance when extending annotations with simplified code
- Improved the performance when large amounts of annotations are missing.
- Implemented a lowMemory option to reduce the memory requirements for very large samples (>100 million reads)


Minor fixes:

- remove the use of get() which looks into environment variables (prone to crashes if a variable of the same name exists) 
and directly references the functions that should be used instead. 
- bug fix when a fa file ois provdied as a string variable to non-windows system
- bug fix when no single exon read class in provided samples
- nug fix when no splice overlaps found between read class and annotations

**bambu version 1.0.2**

Release date: 2020-11-10

- bug fix for author name display
- bug fix for calling fasta file and bam file from ExperimentHub
- update NEWS file 

**bambu version 1.0.0**

Release date: 2020-11-06

- bug fix for parallel computation to avoid bplapply

**bambu version 0.99.4**

Release date: 2020-08-18 

- remove codes using seqlevelStyle to allow customized annotation
- update the requirement of R version and ExperimentHub version

**bambu version 0.3.0**     

Release date: 2020-07-27 

- bambu now runs on windows with a fasta file
- update to the documentation (vignette)
- prepareAnnotations now works with TxDb or gtf file
- minor bug fixes

**bambu version 0.2.0**

Release date: 2020-06-18

**bambu version 0.1.0**

Release date: 2020-05-29 

### Citation
A manuscript describing bambu is currently in preparation. If you use bambu for your research, please cite using the following doi: 10.18129/B9.bioc.bambu. Please specificy that you are using a pre-publication release.

### Contributors

This package is developed and maintained by [Ying Chen](https://github.com/cying111), [Andre Sim](https://github.com/andredsim), [Yuk Kei Wan](https://github.com/yuukiiwa), and  [Jonathan Goeke](https://github.com/jonathangoeke) at the Genome Institute of Singapore. If you want to contribute, please leave an issue. Thank you.

<img src="figures/bambu_design_highres.png" title="Bambu" alt="Bambu">
