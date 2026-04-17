#' @import methods
#' @importFrom Matrix Matrix
#' @importFrom data.table data.table
#' @importClassesFrom Matrix Matrix
#' @importClassesFrom data.table data.table
setClass("quantData",
    slots = c(
        sampleData          = "data.frame",
        readClassDt         = "data.table",
        incompatibleCounts  = "sparseMatrix",
        distTable           = "ANY",
        readToTranscriptMap = "ANY"
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

#' @noRd
setGeneric("getSampleData",          function(x) standardGeneric("getSampleData"))
#' @noRd
setGeneric("getReadClassDt",         function(x) standardGeneric("getReadClassDt"))
#' @noRd
setGeneric("getIncompatibleCounts",  function(x) standardGeneric("getIncompatibleCounts"))
#' @noRd
setGeneric("getDistTable",           function(x) standardGeneric("getDistTable"))
#' @noRd
setGeneric("getReadToTranscriptMap", function(x) standardGeneric("getReadToTranscriptMap"))

#' @noRd
setMethod("getSampleData",          "quantData", function(x) x@sampleData)
#' @noRd
setMethod("getReadClassDt",         "quantData", function(x) x@readClassDt)
#' @noRd
setMethod("getIncompatibleCounts",  "quantData", function(x) x@incompatibleCounts)
#' @noRd
setMethod("getDistTable",           "quantData", function(x) x@distTable)
#' @noRd
setMethod("getReadToTranscriptMap", "quantData", function(x) x@readToTranscriptMap)
