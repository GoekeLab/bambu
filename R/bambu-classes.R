#' @import methods
#' @importFrom Matrix Matrix
#' @importFrom data.table data.table
#' @importClassesFrom Matrix Matrix
#' @importClassesFrom data.table data.table
#' @export
setClass("quantData",
    slots = list(
        sampleData = "data.frame",
        readClassDt = "data.table",
        incompatibleCounts = "ANY",
        nonuniqueCounts = "ANY",
        distTable = "ANY",
        readToTranscriptMap = "ANY"
    )
)

#' @export
setMethod("$", signature(x = "quantData"), function(x, name) {
    slot(x, name)
})
