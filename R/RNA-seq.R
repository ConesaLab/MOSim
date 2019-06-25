#' @include Simulator.R
NULL

# setMethod("initialize", signature="SimRNAseq", function(.Object, ...) {
#     .Object <- callNextMethod(.Object, ...)
#
#     return(.Object)
# })
#
# setMethod("initializeData", signature="SimRNAseq", function(object, simulation) {
#
#     message("Initializing RNA-seq data.")
#
#     # Call parent method AFTER correcting zero counts
#     object <- callNextMethod()
#
#     return(object)
# })
#
#
# setMethod("postSimulation", signature="SimRNAseq", function(object, simulation) {
#     object <- callNextMethod()
#
#     return(object)
# })

setMethod("IDfromGenes", signature="SimRNAseq", function(object, geneNames) {
    return(geneNames)
})

setMethod("IDtoGenes", signature="SimRNAseq", function(object, idNames) {
    return(idNames)
})
