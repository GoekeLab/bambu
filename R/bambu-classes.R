#' @import methods
#' @importFrom Matrix Matrix
#' @importFrom data.table data.table
#' @importClassesFrom Matrix Matrix
#' @importClassesFrom data.table data.table

# Valid values for metadata(se)$seType, which identifies the SE variant:
#   EMCounts     — transcript-level SE with EM-estimated counts (main bambu output)
#   geneCounts   — gene-level SE derived by collapsing transcript counts
#   uniqueCounts — transcript-level SE with uniquely-assigned counts only, no EM
SE_TYPES <- c(
    EMCounts     = "EMCounts",
    geneCounts   = "geneCounts",
    uniqueCounts = "uniqueCounts"
)

# quantData holds per-sample intermediate results from the assignDist step.
# distTable and readToTranscriptMap are optional: populated only when
# returnDistTable=TRUE or trackReads=TRUE respectively, and lifted into
# metadata(countsSe) at the end of bambu().
setClass("quantData",
    slots = c(
        sampleData          = "data.frame",   # per-sample metadata (id, sampleName)
        readClassDt         = "data.table",   # read-class-level count and assignment data
        incompatibleCounts  = "sparseMatrix", # counts of reads incompatible with any annotation
        distTable           = "ANY",          # read-class-to-transcript compatibility table (DataFrame or NULL)
        readToTranscriptMap = "ANY"           # per-read assignment to transcripts (tibble or NULL)
    ),
    prototype = list(
        sampleData  = data.frame(id = integer(), sampleName = character()),
        readClassDt = data.table::data.table()
    ),
    validity = function(object) {
        errs <- character()
        if (!all(c("id", "sampleName") %in% names(object@sampleData)))
            errs <- c(errs, "sampleData must have columns 'id' and 'sampleName'")
        if (!is.null(object@distTable) && !is(object@distTable, "DataFrame"))
            errs <- c(errs, "distTable must be a DataFrame or NULL")
        if (!is.null(object@readToTranscriptMap) && !is(object@readToTranscriptMap, "tbl_df"))
            errs <- c(errs, "readToTranscriptMap must be a tibble or NULL")
        if (length(errs)) errs else TRUE
    })

#' Construct a quantData object
#' @noRd
constructQuantData <- function(sampleData, readClassDt,
                               incompatibleCounts,
                               distTable           = NULL,
                               readToTranscriptMap = NULL) {
    new("quantData",
        sampleData          = sampleData,
        readClassDt         = readClassDt,
        incompatibleCounts  = incompatibleCounts,
        distTable           = distTable,
        readToTranscriptMap = readToTranscriptMap)
}

#' @export
setGeneric("getSampleData",          function(x) standardGeneric("getSampleData"))
#' @noRd
setGeneric("getReadClassDt",         function(x) standardGeneric("getReadClassDt"))
#' @noRd
setGeneric("getIncompatibleCounts",  function(x) standardGeneric("getIncompatibleCounts"))
#' @noRd
setGeneric("getDistTable",           function(x) standardGeneric("getDistTable"))
#' @noRd
setGeneric("getReadToTranscriptMap", function(x) standardGeneric("getReadToTranscriptMap"))

#' @exportMethod getSampleData
setMethod("getSampleData",          "quantData", function(x) x@sampleData)
#' @noRd
setMethod("getReadClassDt",         "quantData", function(x) x@readClassDt)
#' @noRd
setMethod("getIncompatibleCounts",  "quantData", function(x) x@incompatibleCounts)
#' @noRd
setMethod("getDistTable",           "quantData", function(x) x@distTable)
#' @noRd
setMethod("getReadToTranscriptMap", "quantData", function(x) x@readToTranscriptMap)
