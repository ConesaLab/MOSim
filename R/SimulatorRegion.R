#' @include Simulator.R
NULL

setMethod("initialize", signature="SimulatorRegion", function(.Object, ...) {
    .Object <- callNextMethod(.Object, ...)

    message("Configuring simulator ", .Object@name)

    # Split rownames using an specified character, with the convention:
    # chr_start_end
    splitRegions <- function(rLocs) {
        rLocs <- unique(rLocs)

        rData <- strsplit(rLocs, .Object@splitChar)

        rData <- matrix(unlist(rData), ncol = length(rData[[1]]), byrow = TRUE)

        return(data.frame(
            chr = rData[, 1],
            start = as.numeric(rData[, 2]),
            end = as.numeric(rData[, ncol(rData)]),
            ID = rLocs # Keep the rownames as a column for easier joins
        , row.names = rLocs, stringsAsFactors = FALSE))
    }

    # Check rownames format
    if (is.declared(.Object@data)) {

        rLocs <- rownames(.Object@data)

        # If not null, the rownames must follow the scheme:
        # chr_start_<end>
        #
        # Being end optional when there's no associated range, only the
        # location of 1 pb.
        if (is.null(rLocs))
            stop("Rownames on data are required on simulators with region names.")

        .Object@locs <- splitRegions(rLocs)
    } else if (is.declared(.Object@locs)) {

        if (! is.data.frame(.Object@locs)) {
            .Object@locs <- splitRegions(.Object@locs)
        }

        # if (class(.Object@locs) == 'character') {
        #     # Filename
        #     if (length(.Object@locs) == 1) {
        #         .Object@locs <- import.bed(.Object@locs)
        #     } else {
        #         # Make a list from the vector
        #         .Object@locs <- .splitRegions(.Object@locs)
        #     }
        # }

        # if (class(.Object@locs) == 'list') {
        #     .Object@locs <-  GRanges(seqnames=Rle(.Object@locs$chr),
        #                             ranges=IRanges(
        #                                 start=.Object@locs$start,
        #                                 end=.Object@locs$end),
        #                             strand=Rle(strand('*')))
        # }

        # .Object@locs <- sort(.Object@locs)
        # .Object@chGRanges <- split(.Object@locs, seqnames(.Object@locs))
        # .Object@nLocs <- length(.Object@locs)

        # smessage(.Object@nLocs, " locations from ", length(.Object@chGRanges), " seqnames.")
    } else if (is.declared(.Object@idToGene)) {
        rLocs <- .Object@idToGene$ID

        # If not null, the rownames must follow the scheme:
        # chr_start_<end>
        #
        # Being end optional when there's no associated range, only the
        # location of 1 pb.
        if (is.null(rLocs))
            stop("Rownames on data are required on simulators with region names.")

        .Object@locs <- splitRegions(rLocs)
    } else {
        stop("No regions provided.")

        # .Object@chGRanges <- list(NULL)
        # .Object@locs <- NULL
    }

    return(.Object)
})

# setValidity("SimulatorRegion", function(object) {
#     return(getValidity(getClassDef("Simulator"))(object))
# })


# .SimulatorRegion.regionNames <- function(object, chrNumber, start, end=NULL) {
#     end <- if(is.null(end)) start else end
#
#     paste(chrNumber, start, end, sep="_")
# }

# .SimulatorRegion.getLengths <- function(object) {
#     width(object@locs)
# }


# setMethod("regionNames", signature="SimulatorRegion", .SimulatorRegion.regionNames)
# setMethod("locGRanges", signature="SimulatorRegion", .SimulatorRegion.locGRanges)
# setMethod("nearestGenes", signature="SimulatorRegion", .SimulatorRegion.nearestGenes)
# setMethod("getLengths", signature="SimulatorRegion", .SimulatorRegion.getLengths)
