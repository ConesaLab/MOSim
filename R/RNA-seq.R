#' @include Simulator.R
NULL

setMethod("initialize", signature="SimRNAseq", function(.Object, ...) {
    .Object <- callNextMethod(.Object, ...)

    return(.Object)
})

setMethod("initializeData", signature="SimRNAseq", function(object, simulation) {

    message("Initializing RNA-seq data.")

    # Replace expressed genes with a value of zero with random values from
    # zero to max value.
    # TODO: disabled this because randomCounts will be generated later
    # message("Setting random counts for expressed genes with 0 base counts.")
    #
    # zeroIndexes <- (! object@data) &
    #     rownames(object@data) %in% simulation@simSettings$geneSamples$expr
    #
    # # TODO: make sure to assign a value greater than 0 or trust this?
    # object@data[zeroIndexes] <- sample.int(max(object@data), size = sum(zeroIndexes))

    # Call parent method AFTER correcting zero counts
    object <- callNextMethod()

    # # Ensure change between groups for DE genes with flat profiles
    # if (simulation@numberGroups > 1) {
    #     # onlyFlat <- (length(simulation@times) > 1)
    #
    #     # Count range to exclude comparing with the reference (group 1)
    #     excludeMargin <- 10
    #
    #     # Select only DE genes
    #     exprGenesDE <- simulation@simSettings$geneSamples$DE
    #
    #     profileSubset <- dplyr::filter(simulation@simSettings$geneProfiles[[class(object)]], ID %in% exprGenesDE)
    #     dataSubset <- object@data[as.character(profileSubset$ID), ]
    #
    #     # Determine the rows with the same condition on all columns
    #     # Change only base count values when the profiles are flat on all conditions,
    #     # otherwise they are already DE genes.
    #     sameCond <- apply(apply(dplyr::select(profileSubset, starts_with("Group")), 2, '==', 'flat'), 1, all)
    #
    #     # sameCond <-
    #     #     apply(
    #     #         apply(
    #     #             dplyr::select(profileSubset, starts_with("Group"), - Group1),
    #     #             2,
    #     #             '==',
    #     #             dplyr::select(profileSubset, Group1)
    #     #         ),
    #     #         1,
    #     #         all
    #     #     )
    #
    #     if ((nSameCount <- sum(sameCond))) {
    #         message(sprintf("Changing count values on %d DE genes with the same profile on all groups.", nSameCount))
    #
    #         # Replace counts on a random set of groups, assigning a random count excluding a range around group 1 counts
    #         countRange <- 0:object@max
    #
    #         dataSubset[sameCond, ] <- t(apply(dataSubset[sameCond, ], 1, function(row) {
    #             colChange <- sample(2:simulation@numberGroups, size = sample.int(simulation@numberGroups - 1, 1))
    #
    #             # Remove similar counts to first group, excluding negative numbers.
    #             excludeCounts <- max(row[1] - excludeMargin, 0):(row[1] + excludeMargin)
    #
    #             # Position 1: count 0
    #             row[colChange]  <- sample(countRange[-(excludeCounts + 1)], length(colChange), replace = TRUE)
    #
    #             return(row)
    #         }))
    #
    #         object@data[rownames(dataSubset), ] <- dataSubset
    #     }
    # }

    return(object)
})


setMethod("postSimulation", signature="SimRNAseq", function(object, simulation) {
    object <- callNextMethod()

    # Apply only to RNA-seq, not child classes
    # TODO: disabled this
    # if (! object@regulator) {
    #     # Set all non-expressed genes counts to zero
    #     nonExprGenes <- simulation@simSettings$geneSamples$nonExpr
    #
    #     message(sprintf("Adding %d non-expressed genes to simulated data with zero counts.", length(nonExprGenes)))
    #
    #     object@simData[nonExprGenes, ] <- 0
    # }

    return(object)
})

setMethod("IDfromGenes", signature="SimRNAseq", function(object, geneNames) {
    return(geneNames)
})

setMethod("IDtoGenes", signature="SimRNAseq", function(object, idNames) {
    return(idNames)
})
