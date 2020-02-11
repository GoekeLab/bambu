## code to prepare `example_data` dataset goes here
library(data.table)

example_data <- list()

a <- 500

example_data[[1]] <-  data.table(tx_sid = c(1,2,3,1,2,1,2,3),
                                 read_class_id = c(1,2,3,rep(4,2),rep(5,3)),
                                 nobs = c(0,0,0,rep(10,2),rep(a,3)),
                                 gene_sid = 1)

example_data[[2]] <- data.table(tx_sid = c(1,2,1,2,3),
                                read_class_id = c(rep(4,2),rep(5,3)),
                                nobs = c(rep(10,2),rep(a,3)),
                                gene_sid = 2)

example_data[[3]] <- data.table(tx_sid = c(1,2),
                                read_class_id = c(1,2),
                                nobs = c(5,10),
                                gene_sid = 3)

example_data[[4]] <- data.table(tx_sid = c(1,2,3,1,2,1,2,3),
                                read_class_id = c(1,2,3,rep(4,2),rep(5,3)),
                                nobs = c(100,5,500,rep(20,2),rep(5,3)),
                                gene_sid = 4)

usethis::use_data("example_data")
