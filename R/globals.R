## quiets concerns of R CMD check re:no binding global varaible or function
if (getRversion() >= "2.15.1") {
    utils::globalVariables(c(
        "GENEID",
        "PC1", "PC2", "TXNAME",
        "annotationTxId", "chr", "chr.rc",
        "compatible", "confidenceType", "endMatch",
        "eqClass", "eqClassReadCount", "equal",
        "estimates", "exon_rank", "feature",
        "geneCount", "geneId", "gene_id", "gene_sid",
        "groupVar", "group_name", "intronEnds",
        "intronStarts", "newGeneClass", "newTxClass",
        "nobs_stored", "queryElementsOutsideMaxDist",
        "queryHits.x", "queryHits.y", "readClassId",
        "readCount", "read_class_id", "read_class_sid",
        "read_class_sid_stored", "rowMaxs", "rowMins",
        "runname", "seqlengths", "startMatch", "strand.rc",
        "subjectCount", "subjectElementsOutsideMaxDist",
        "subjectHits.y", "sum_nobs", "txId", "txNumberFiltered",
        "tx_id", "tx_name", "tx_sid", "uniqueEndLengthQuery",
        "uniqueLengthQuery", "uniqueLengthSubject",
        "uniqueStartLengthQuery", "value", "valueGene", "valueGeneCPM",
        "variable", "junctionMatchList", "readClass.file",
        "start.ptm", "verbose", ".insertGaps", ".","subjectList"
    ))
}
