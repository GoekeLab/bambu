#' @import methods
#' @importFrom Matrix Matrix
#' @importFrom data.table data.table
#' @importClassesFrom Matrix Matrix
#' @importClassesFrom data.table data.table
#' @export
setClass("quantData",
    slots = list(
        sampleData = "data.frame",
        uniqueCounts = "Matrix",
        readClassDt = "data.table",
        incompatibleCountMatrix = "Matrix",
        sampleNames = "character",
        incompatibleCounts = "ANY",
        nonuniqueCounts = "Matrix",
        distTable = "ANY",
        readToTranscriptMap = "ANY"
    )
)

#' @export
setMethod("$", signature(x = "quantData"), function(x, name) {
    slot(x, name)
})
