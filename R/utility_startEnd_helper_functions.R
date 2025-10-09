#helper function to identify and quantify tss and tes after Tx discovery
getTss <- function(grList, width = 20){
  starts <- as.integer(endoapply(start(grList), function(x) x[1]))
  ends <- as.integer(endoapply(end(grList), function(x) x[length(x)]))
  strands <- as.character(getStrandFromGrList(grList))
  seqnames_list <- as.character(getChrFromGrList(grList))
  tss <- ifelse(strands == "+", starts - ceiling(width/2), ends - ceiling(width/2))
  tssGranges <- GRanges(
    seqnames = seqnames_list,
    ranges = IRanges(start = tss, width = width),
    strand = strands
  )
  return(tssGranges)
}

getTes <- function(grList, width = 20){
  starts <- as.integer(endoapply(start(grList), function(x) x[1]))
  ends <- as.integer(endoapply(end(grList), function(x) x[length(x)]))
  strands <- as.character(getStrandFromGrList(grList))
  seqnames_list <- as.character(getChrFromGrList(grList))
  tes <- ifelse(strands == "+", ends - ceiling(width/2), starts - ceiling(width/2))
  tesGranges <- GRanges(
    seqnames = seqnames_list,
    ranges = IRanges(start = tes, width = width),
    strand = strands
  )
  return(tesGranges)
}

cluster_and_assign <- function(gr, prefix) {
  clusters <- reduce(gr, min.gapwidth = 50, ignore.strand = FALSE, with.revmap = TRUE)
  names(clusters) <- paste0(prefix, seq_along(clusters))
  revmap_vec <- unlist(clusters$revmap)
  cluster_vec <- rep(names(clusters), lengths(clusters$revmap))
  cluster_ids <- rep(NA_character_, length(gr))
  cluster_ids[revmap_vec] <- cluster_vec
  list(ids = cluster_ids, clusters = clusters)
}

# Main function
assignGlobalTssTesId <- function(annotations) {
  tss_res <- cluster_and_assign(getTss(annotations), "BambuTss")
  mcols(annotations)$globleTssId <- tss_res$ids
  metadata(annotations)$tss_clusters <- tss_res$clusters
  tes_res <- cluster_and_assign(getTes(annotations), "BambuTes")
  mcols(annotations)$globleTesId <- tes_res$ids
  metadata(annotations)$tes_clusters <- tes_res$clusters
  return(annotations)
}