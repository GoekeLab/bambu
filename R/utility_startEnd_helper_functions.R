#helper function to identify and quantify tss and tes after Tx discovery
getTss <- function(grList, width = 20){
  starts <- as.integer(endoapply(start(grList), function(x) x[1]))
  ends <- as.integer(endoapply(end(grList), function(x) x[length(x)]))
  strands <- as.character(getStrandFromGrList(grList))
  seqnames_list <- as.character(getChrFromGrList(grList))
  tss <- ifelse(strands == "+", starts - floor(width/2), ends - floor(width/2))
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
  tes <- ifelse(strands == "+", ends - floor(width/2), starts - floor(width/2))
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

# Helper: cluster TES anchors (sorted by realEnd) within a single chr/strand group.
# A new cluster is forced when:
#   (a) gap to previous anchor exceeds threshold, OR
#   (b) the same sample appears with a *different* tes_key
#       (different sampleTesId within one sample must not share mulSamTesId)
cluster_tes_anchors <- function(df, threshold) {
  n <- nrow(df)
  if (n == 0L) return(integer(0))
  group_num <- integer(n)
  cur_grp <- 0L
  cur_tes_keys <- df$tes_key[1]
  cur_samples  <- df$sample_idx[1]
  for (k in seq_len(n)[-1]) {
    dist_gap <- df$realEnd[k] - df$realEnd[k - 1]
    same_sample_conflict <- dist_gap > 0 &&
      (df$sample_idx[k] %in% cur_samples) &&
      !(df$tes_key[k] %in% cur_tes_keys)
    if (is.na(dist_gap) || dist_gap > threshold || same_sample_conflict) {
      cur_grp      <- cur_grp + 1L
      cur_tes_keys <- df$tes_key[k]
      cur_samples  <- df$sample_idx[k]
    } else {
      cur_tes_keys <- c(cur_tes_keys, df$tes_key[k])
      cur_samples  <- c(cur_samples, df$sample_idx[k])
    }
    group_num[k] <- cur_grp
  }
  return(group_num)
}

updateTesIdAcrossSamples <- function(readClassList, threshold = 10) {
  # Collect all samples into one dataframe
  master_df <- lapply(seq_along(readClassList), function(i) {
    rd <- as.data.frame(rowData(readClassList[[i]]))
    rd$sample_idx <- i
    rd$original_row <- seq_len(nrow(rd))
    return(rd)
  }) %>% bind_rows()

  master_df <- master_df %>%
    mutate(
      realEnd = ifelse(strand.rc == "+", end.rc, start.rc),
      tes_key = ifelse(!is.na(sampleTesId), as.character(sampleTesId),
                       paste0("end_", realEnd))
    )

  anchors <- master_df %>%
    group_by(chr.rc, strand.rc, sample_idx, tes_key) %>%
    summarise(realEnd = median(realEnd), .groups = "drop") %>%
    arrange(chr.rc, strand.rc, realEnd) %>%
    group_by(chr.rc, strand.rc) %>%
    group_modify(~ { .x$group_num <- cluster_tes_anchors(.x, threshold); .x }) %>%
    ungroup() %>%
    group_by(chr.rc, strand.rc, group_num) %>%
    mutate(tesId_count = n_distinct(sample_idx)) %>%
    ungroup() %>%
    mutate(
      mulSamTesId = paste0(chr.rc, "_", strand.rc, "_tes_", group_num),
      mulSamTesId = ifelse(tesId_count == 1,
                           paste0(mulSamTesId, ".s", sample_idx), mulSamTesId)
    ) %>%
    select(chr.rc, strand.rc, sample_idx, tes_key, mulSamTesId)

  master_df <- master_df %>%
    left_join(anchors, by = c("chr.rc", "strand.rc", "sample_idx", "tes_key")) %>%
    arrange(sample_idx, original_row)

  updated_list <- lapply(seq_along(readClassList), function(i) {
    sample_data <- master_df %>% filter(sample_idx == i) %>% arrange(original_row)
    rowData(readClassList[[i]])$mulSamTesId <- sample_data$mulSamTesId
    return(readClassList[[i]])
  })
  return(updated_list)
}
