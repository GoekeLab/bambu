<a href="https://raw.githubusercontent.com/GoekeLab/bamboo/master/figures/transparent-bamboo.png?token=AGA7DTCYEJWLSUX5XFFRHCS6UZYAE"><img src="https://raw.githubusercontent.com/GoekeLab/bamboo/master/figures/transparent-bamboo.png?token=AGA7DTCYEJWLSUX5XFFRHCS6UZYAE" title="Bamboo" alt="Bamboo"></a>

# Bamboo: genomic alignments based long read isoform reconstruction and quantification

This R package helps build transcript models from read-to-genome alignments while discovering novel transcripts and perform transcript quantification for long read RNA-Seq data. 



[![bamboo](https://img.shields.io/badge/bamboo-version0.1.0-brightgreen)](https://github.com/GoekeLab/bamboo) [![Maintained?](https://img.shields.io/badge/Maintained%3F-Yes-blue)](https://gemnasium.com/badges/badgerbadgerbadger)  [![Install](https://img.shields.io/badge/Install-Github%20R%20package-brightgreen)](https://github.com/badges/badgerbadgerbadger/issues) [![License](https://img.shields.io/badge/license-GPL--3.0-brightgreen)](http://badges.mit-license.org)

---


## Installation

- You may install it from github, or clone this repo to your local machine using `https://github.com/GoekeLab/bamboo.git`

>  install the dev package first if not installed 

```rscript
install.packages("devtools")
```

> now install bamboo packages

```rscript
library(devtools)
install_github("GoekeLab/bamboo")
```
---

## Usage 

We provide two modes of transcript quantification. 
 - The first mode is to quantify transcript using the provided annotation information. In this mode, all reads will be matched to provided annotation and contribute to the quantification.
 - The second mode will extend annotation with new annotations combined across provided samples. Quantification will then be performed based on the extended annotations. 
 
 Below we will present an example of the basic usage of bamboo function:
 ```rscript
 # Load bam file
 test.bam <- system.file("extdata", "GIS_HepG2_cDNAStranded_Rep5-Run4_chr9_108865774_109014097.bam", package = "bamboo")
 
 # Load fa file
  fa.file <- system.file("extdata", "Homo_sapiens.GRCh38.dna_sm.primary_assembly_chr9.fa.gz", package = "bamboo")
 
 Load txdb object
samplefile <- "http://s3.ap-southeast-1.amazonaws.com/ucsc-trackdata.store.genome.sg/chenying/bamboo_exampleDataset/Homo_sapiens.GRCh38.91.annotations-txdb.sqlite"

# download remote sqlite file as a temporary object
dbfile <- tempfile(fileext=".sqlite")
download.file(samplefile, dbfile)
txdb <- AnnotationDbi::loadDb(dbfile)

# Perform bamboo 
bamboo(bam.file = test.bam, txdb = txdb, genomeSequence = fa.file)

# Extend annotations and perform bamboo on extended annotations
bamboo(bam.file = test.bam, txdb = txdb, genomeSequence = fa.file, extendAnnotations = TRUE)
        
```


---

## Contributors

This package is developed and maintained by [Jonathan Goeke](https://github.com/jonathangoeke) and [Ying Chen](https://github.com/cying111). If you want to contribute, please leave an issue. Thank you.

## License

[![License](http://img.shields.io/:license-mit-blue.svg?style=flat-square)](http://badges.mit-license.org)

- **[GPL-3.0](https://opensource.org/licenses/GPL-3.0)**
