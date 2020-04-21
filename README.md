<a href="https://raw.githubusercontent.com/GoekeLab/bamboo/master/figures/transparent-bamboo.png?token=AGA7DTCYEJWLSUX5XFFRHCS6UZYAE"><img src="https://raw.githubusercontent.com/GoekeLab/bamboo/master/figures/transparent-bamboo.png?token=AGA7DTCYEJWLSUX5XFFRHCS6UZYAE" title="Bamboo" alt="Bamboo"></a>

# bamboo: reference-guided transcript discovery and quantification for long read RNA-Seq data


bamboo is a R package for multi-sample transcript discovery and quantification using long read RNA-Seq data. You can use bamboo after read alignment to obtain expression estimates for known and novel transcripts and genes. The output from bamboo can directly be used for visualisation and downstream analysis such as differential gene expression or transcript usage.



[![bamboo](https://img.shields.io/badge/bamboo-v0.9.0-brightgreen)](https://github.com/GoekeLab/bamboo) [![Maintained?](https://img.shields.io/badge/Maintained%3F-Yes-blue)](https://gemnasium.com/badges/badgerbadgerbadger)  [![Install](https://img.shields.io/badge/Install-Github-brightgreen)](https://github.com/badges/badgerbadgerbadger/issues) 
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

---


## Installation

You can install bamboo from github:

```rscript
install.packages("devtools")
devtools::install_github("GoekeLab/bamboo")
```
---

## General Usage 

The default mode to run bamboo is using a set of aligned reads (bam files), reference genome annotations (gtf file, TxDb object, or bambooAnnotation object), and reference genome sequence (fasta file or BSgenome). Bamboo will return a summarizedExperiment object with the genomic coordinates for annotated and new transcripts and genes and with the transcript expression estimates: 
 ```rscript
 
bamboo(reads, annotations, genomeSequence,...)

library(bamboo)
test.bam <- system.file("extdata", "GIS_HepG2_cDNAStranded_Rep5-Run4_chr9_108865774_109014097.bam", package = "bamboo")

se <- bamboo(reads = test.bam, annotations = "TxDb.Hsapiens.UCSC.hg38.knownGene", genomeSequence = "BSgenome.Hsapiens.NCBI.GRCh38")
       
```


We highly recommend to use the same annotations that were used for genome alignment. If you have a gtf file and fasta file you can run bamboo with the following options:

```rscript
test.bam <- system.file("extdata", "GIS_HepG2_cDNAStranded_Rep5-Run4_chr9_108865774_109014097.bam", package = "bamboo")

 fa.file <- system.file("extdata", "Homo_sapiens.GRCh38.dna_sm.primary_assembly_chr9.fa.gz", package = "bamboo")
gtf.file <- [...]
bambooAnnotations <- prepareAnnotationsFromGtf(gtf file)
se <- bamboo(reads = test.bam, annotations=bambooAnnotations, genomeSequence ="BSgenome.Hsapiens.NCBI.GRCh38")

```
---


## Other Use Cases
**Use precalculated annotation objects**

If you plan to run bamboo more frequently, you can store a bambooAnnotations object.

```rscript
library(bamboo)
test.bam <- system.file("extdata", "GIS_HepG2_cDNAStranded_Rep5-Run4_chr9_108865774_109014097.bam", package = "bamboo")

hg38 <- readRDS(url(‘...’))

se <- bamboo(reads = test.bam, annotations=hg38, genomeSequence = fa.file)

```

The bambooAnnotation object can be calculated from a gtf file:
```rscript
annotations <- prepareAnnotationFromGTF(gtf.file)
```

From TxDb object
```rscript
annotations <- prepareAnnotations(txdb)
```

You can also download annotations for recent genome releases here:
-  hg38:  ftp://ftp.ensembl.org/pub/release-99/gtf/homo_sapiens          
-  mm10:  ftp://ftp.ensembl.org/pub/release-99/gtf/mus_musculus     

**Quantification of annotated transcripts and genes only (no transcript/gene discovery)**

```rscript
bamboo(reads = test.bam, annotations = txdb, genomeSequence = fa.file, extendAnnotations = FALSE)
```

**Large sample number/ limited memory**     
For larger sample numbers we recommend to write the processed data to a file:
```rscript
bamboo(reads, outputReadClassFolder, genomeSequence, annotations)
```

---

## Advanced Options

**More strigent filtering thresholds imposed on potential novel transcripts**    
- For example   
> Keep novel transcripts with min 5 read count in at least 1 sample:  
```rscript
bamboo(reads, genomeSequence, annotations, ir.control = list(min.readCount = 5))
```

> Keep novel transcripts with min 5 samples having at least 2 counts:

```rscript
bamboo(reads, genomeSequence, annotations, ir.control = list(min.sampleNumber = 5))
```

> Filter out transcripts with relative abundance within gene lower than 10%: 
```rscript
bamboo(reads, genomeSequence, annotations, ir.control = list(min.sampleNumber = 5))
```

See more details in the manual.

**Quantification without bias correction**     
> The default estimation automatically does bias correction for expression estimates. However, you can choose to perform the quantification without bias correction.    
```rscript
bamboo(reads, genomeSequence, annotations, algo.control(bias_correction = FALSE))
```

**Parallel computation**      
> bamboo also allows parallel computation for EM.    
```rscript
bamboo(reads, genomeSequence, annotations, algo.control(ncore = 8))
```
---

## Contributors

This package is developed and maintained by [Jonathan Goeke](https://github.com/jonathangoeke) and [Ying Chen](https://github.com/cying111) at Genome Institute of Singapore. If you want to contribute, please leave an issue. Thank you.

## License

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

- **[GPL v3](https://www.gnu.org/licenses/gpl-3.0)**
