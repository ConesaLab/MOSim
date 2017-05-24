#' @include AllGeneric.R
NULL

#' @rdname initialize-methods
#' @aliases initialize,Simulator
setMethod("initialize", signature="Simulator", function(.Object, ...) {
    .Object <- callNextMethod(.Object, ...)

    # Casting to dataframe
    .Object@idToGene <- as.data.frame(.Object@idToGene, stringsAsFactors = FALSE)
    .Object@data <- as.data.frame(.Object@data, stringsAsFactors = FALSE)

    # Convert data frames to characters or numbers
    .Object@idToGene[] <- lapply(.Object@idToGene, as.character)
    .Object@data[] <- lapply(.Object@data, as.numeric)

    # # If required, adjust the total number of features to simulate on every omic
    # # sim@totalFeatures <- min(nrow(sim@data), sim@totalFeatures)
    # # TODO: if the object is already initialized this will not be changed
    # # according to the new passed options
    # if (nrow(.Object@data) > .Object@totalFeatures) {
    #
    #     if (.Object@regulator) {
    #         # TODO: limit the sample to those regulators affecting the possibly
    #         # limited gene list?
    #     }
    #
    #     # TODO: keep copy of the original data?
    #     .Object@data <- dplyr::sample_n(.Object@data, .Object@totalFeatures)
    # }

    return(.Object)
})

#' initializeData
#'
#' Generates the initial sample of data. It copies the initial data for every
#' group of the experimental design, and ensures to change the count values
#' for those genes (and associated regulators with an active effect on it)
#' that are marked as DE with only flat profiles in all groups.
#'
#' @param object Instance of \linkS4class{Simulator} class.
#' @param simulation Initialized instance of \linkS4class{Simulation} class
#'
#' @return An object of class \linkS4class{Simulator} with the @data slot filled
#' with the correct initial data structure and values.
#' @export
#'
setMethod("initializeData", signature="Simulator", function(object, simulation) {

    # Adjust sequencing depth
    object@depth <- object@depth*10^6

    # Set min/max values based on provided data
    # TODO: calculate this AFTER adjusting the total number of features to simulate?
    if (is.declared(object@data)) {
        object@min <- min(object@data)
        object@max <- max(object@data)

        dataPerc <- quantile(object@data$Counts, object@minQuantile)

        # Minimum value to change when simulating data
        # TODO: this does not use indexes anymore but numeric positions instead
        object@increment <- (dataPerc[2] - dataPerc[1]) * 0.1
    }


    # Clone the base counts for every group
    if (ncol(object@data) < simulation@numberGroups) {
        object@data <- object@data[, rep.int(1, simulation@numberGroups), drop = FALSE]
    }

    # Important: column names pattern "Counts.Group" is required .
    colnames(object@data) <- paste('Counts.Group', seq(simulation@numberGroups))

    if (! is.null(flatProfiles <- simulation@simSettings$geneProfiles$FlatGroups)) {

        # Structure: Gene_ID | Group_1 | ... | Group_N
        # onlyFlat <- (length(simulation@times) > 1)

        # Count range to exclude comparing with the reference (group 1)
        excludeMargin <- 10

        # Different treatments between RNA-seq and the rest of simulators
        if ( ! object@regulator) {
            # Select genes
            # dataSubset <- object@data[flatProfiles$ID, ]

            # Iterate over the rows on a unique matrix with the structure:
            # Gene_ID | Profile_Group_1 | ... | Counts_Group_1 | ...
            object@data[flatProfiles$ID, ] <- t(apply(cbind(flatProfiles, object@data[flatProfiles$ID, ]), 1, function(idRow) {
                # Types of profiles (exclude first column)
                idProfiles <- idRow[seq(2, length.out = simulation@numberGroups)]
                # Count values (columns after profiles)
                idCounts <- as.numeric(idRow[seq(simulation@numberGroups + 2, length.out = simulation@numberGroups)])

                # Profiles
                return(ifelse(idProfiles == 'flat',
                       # If the profile is flat, keep the current value
                       idCounts,
                       ifelse(idProfiles == 'enhancer',
                              # If the profile is enhancer, generate a new count value in the range
                              # (actual_value + excludeMargin) to maxValue
                              rep.int(runif(1, min = min(idCounts[1] + excludeMargin, object@max), max = object@max), simulation@numberGroups),
                              # If the profile is repressor, generate a new count value in the range
                              # minValue to (actual_value - excludeMargin)
                              rep.int(runif(1, min = object@min, max = max(idCounts[1] - excludeMargin, object@min)), simulation@numberGroups)
                              )
                ))
            }))
        } else {
            # In case of regulators, the initial value will depend on the effect
            # on the associated gene.
            profileSubset <- dplyr::rename(flatProfiles, Gene = ID) %>%
                dplyr::inner_join(simulation@simSettings$geneProfiles[[class(object)]][, c('ID', 'Effect', 'Gene')], by = c("Gene" = "Gene")) %>%
                dplyr::select(ID, Effect, dplyr::starts_with("Group")) %>%
                dplyr::filter(! is.na(Effect))

            if (nrow(profileSubset)) {
                object@data[profileSubset$ID, ] <- t(apply(cbind(profileSubset, object@data[profileSubset$ID, ]), 1, function(idRow) {
                    # Regulator effect
                    regEffect <- idRow[2]
                    # Types of profiles (skip first 2 columns: ID and Effect)
                    idProfiles <- idRow[seq(3, length.out = simulation@numberGroups)]
                    # Count values (following columns)
                    idCounts <- as.numeric(idRow[seq(simulation@numberGroups + 3, length.out = simulation@numberGroups)])

                    # Profiles
                    return(ifelse(idProfiles == 'flat',
                                  # If the profile is flat, keep the current value
                                  idCounts,
                                  # If the regulator effect is enhancer, keep the same behaviour,
                                  # if not, assign the opposite.
                                  ifelse(idProfiles == regEffect,
                                         # If the profile is enhancer, generate a new count value in the range
                                         # (actual_value + excludeMargin) to maxValue
                                         rep.int(runif(1, min = min(idCounts[1] + excludeMargin, object@max), max = object@max), simulation@numberGroups),
                                         # If the profile is repressor, generate a new count value in the range
                                         # minValue to (actual_value - excludeMargin)
                                         rep.int(runif(1, min = object@min, max = max(idCounts[1] - excludeMargin, object@min)), simulation@numberGroups)
                                  )
                    ))
                }))
            }
        }
    }

    return(object)
})

#' @rdname simulate-methods
#' @aliases simulate,Simulation
setMethod("simulate", signature="Simulator", function(object, simulation) {

    message(sprintf("Starting simulation of %s.", object@name))

    object <- initializeData(object, simulation)

    makeReplicates <- function(counts, groupInfo, timeInfo) {

        message(sprintf("\t- Making replicates for group %d on time %s.", groupInfo, timeInfo))

        # Allowing for some noise in the NB mean
        nbNoise <- if (exists('NB', object@noiseParams)) object@noiseParams$NB else object@noiseParams[[1]]

        # if (nbNoise > 0) {
        #     message("\t(NB noise enabled ", nbNoise, ")")
        #     counts.noise <- t(sapply(counts, function(x) { x + c(-1,1)*nbNoise*x }))
        #     counts.noise[counts.noise < 0] <- 0
        #
        #     mu.noise <- apply(counts.noise, 1, function(x) { runif(1, x[1], x[2]) })
        # } else {
        #     message("\t(NB noise disabled)")
        #     mu.noise <- counts.noise <- counts
        # }

        # message("\t(NB noise disabled)")
        # logging::loginfo("\t(NB noise disabled)")

        mu.noise <- counts.noise <- counts

        # Transform to CPM (rgamma tables based on CPM bins)
        mu.noise.cpm <- 10^6 * mu.noise / sum(mu.noise)

        # TODO: disable cpm
        # mu.noise.cpm <- mu.noise

        # Replace NA for 0
        mu.noise.cpm[is.na(mu.noise.cpm)] <- 0.1
        mu.noise.cpm[mu.noise.cpm < 0.1] <- 0.1

        # Counts (not CPM)
        mu.noise[is.na(mu.noise)] <- 0.1
        mu.noise[mu.noise < 0.1] <- 0.1

        # Calculate stdev using gamma tables
        stdev.cpm <- sapply(mu.noise.cpm, function(x) {
            # Select row of table using mean column as intervals
            # TODO: revert the +1 and max
            # rgammaIndex <- min(which.max(x <= object@gammaTable$mean) + 1, length(object@gammaTable$mean))
            rgammaIndex <- which.max(x <= object@gammaTable$mean)

            rgamma(1, shape=object@gammaTable$shape[rgammaIndex], scale=object@gammaTable$scale[rgammaIndex])
        })

        # TODO: remove
        # stdev.which <- sapply(mu.noise.cpm, function(x) which.max(x <= object@gammaTable$mean))
        stdev <- sum(mu.noise) * stdev.cpm/1E6

        # Apply to every row:
        #   Param x: [mu.noise.cpm, stdev]
        # temp <- t(apply(cbind(mu.noise.cpm, stdev.cpm), 1, function (x) {
        temp <- t(matrix(apply(cbind(mu.noise, stdev), 1, function (x) {
            if (x[1] == x[2]^2) {
                replis <- rpois(n = simulation@numberReps, lambda = x[1])
            } else {
                replis <- rnbinom(n = simulation@numberReps, size = x[1]^2/abs(x[2]^2 - x[1]),  mu = x[1])
            }

            return(replis)
        }), nrow = simulation@numberReps))

        # TODO: reverse CPM here?
        # Adjust to depth
        # Note: temp is a matrix, take into account if this is modified with
        # dynamic values.

        # TODO: change this
        temp.depth <- apply(temp, 2, function(x) x * object@depth / sum(x))
        # temp.depth <- temp*object@depth/1E6
        # temp.depth <- temp

        object@debugInfo <<- simDebug(
            object,
            simulation,
            method = sprintf("makeReplicates[group=%d,time=%d]", groupInfo, timeInfo),
            counts.noise,
            mu.noise,
            mu.noise.cpm,
            stdev,
            stdev.cpm,
            temp,
            temp.depth
            # my.size.nb
        )

        return(temp.depth)
    }

    # Profile table of simulator
    simProfiles <- simulation@simSettings$geneProfiles[[class(object)]]

    # Select only the active rows of the regulator profile table, discarding
    # the duplicated rows (all active rows for a given regulator will have the
    # same profile among groups by now).
    if (object@regulator) {
        # If the object is pregenerated (e.g. methylation) then keep the profiles
        # as they are.
        if (object@pregenerated) {
            simProfiles <- dplyr::select(simProfiles, ID, dplyr::starts_with("Group")) %>% dplyr::distinct()
        } else {
            simProfiles <- dplyr::filter(simProfiles, ! is.na(Effect)) %>%
                dplyr::select(ID, dplyr::starts_with("Group")) %>%
                dplyr::distinct_()

            # Add remaining IDs present on data with a flat profile
            simProfiles <- rbind(simProfiles, do.call(cbind, setNames(append(
                # ID
                list(rownames(object@data)[! rownames(object@data) %in% as.character(simProfiles$ID)]),
                # Groups
                rep('flat', simulation@numberGroups)
            ), colnames(simProfiles))))
        }
    }


    # If data is already generated (i.e. methylation simulator) skip some
    # steps like generating random counts or replicates, adjusting also the
    # parameters for mapply.

    # Pass only "group" columns
    columnsParam <- dplyr::select(simProfiles, dplyr::starts_with("Group"))
    # Keep track of iteration
    iterationParam <- seq(simulation@numberGroups)
    # Pass IDs
    idsParam <- simProfiles[, rep('ID', simulation@numberGroups), drop = FALSE]

    if (object@pregenerated) {
        # Repeat for every replicate
        # columnsParam <- do.call(cbind, replicate(simulation@numberReps, columnsParam, simplify = FALSE))
        # columnsParam <- do.call(cbind, lapply(columnsParam, function(x) replicate(simulation@numberReps, x, simplify=TRUE)))
        columnsParam <- columnsParam[, rep(colnames(columnsParam), each = simulation@numberReps), drop = FALSE]
        # iterationParam <- rep(iterationParam, times = simulation@numberReps)
        iterationParam <- rep(iterationParam, each = simulation@numberReps)
        idsParam <- simProfiles[, rep('ID', simulation@numberGroups * simulation@numberReps), drop = FALSE]
    }

    # Generate time series (if any) with replicates
    object@simData <- data.frame(
        mapply(
            function(counts, profiles, group, ids) {
                if (! object@pregenerated) {
                    message(sprintf("- Simulating count values for group %d.", group))
                } else {
                    message(sprintf("[Pregenerated replicate] Simulating count values for group %d.", group))
                }

                # Create a matrix of <times> columns and <ID> rows, with every
                # row being a vector of the associated profile.
                timeVectors <-
                    matrix(
                        unlist(simulation@simSettings$profiles[as.character(profiles)], use.names = FALSE),
                        ncol = length(simulation@times),
                        byrow = TRUE
                    )

                # Generate base counts using the formulas:
                # X = initial_counts + noise [for flat]
                # X = m + lambda(M - m) + noise [for lambda in [0,1]]
                # X = M + lambda(M - m) + noise [for lambda in [-1,0]]

                # Delegate creation of specific simulation parameters
                simulateParams <- simulateParams(object, simulation, counts, profiles, group, ids)

                randomCounts <- simulateParams$randomCounts
                noiseValues <- simulateParams$noiseValues
                M <- simulateParams$M
                m <- simulateParams$m

                # Rows with flat profile
                indexFlat <- (profiles == "flat")

                # Simulate values
                # Make sure that the result M - m is at least
                # (P90 - P10) * 0.1 of the initial data
                simData <-
                    ifelse(grepl('repression', profiles, fixed = TRUE), M, m) +
                    (ifelse((M-m) < object@increment, object@increment, M - m) * timeVectors) +
                    noiseValues

                # Overwrite flat rows with the correct value
                # TODO: use randomCounts or counts?
                # simData[indexFlat, ] <- counts[indexFlat] + noiseValues[indexFlat, ]
                simData[indexFlat, ] <- randomCounts[indexFlat] + noiseValues[indexFlat, ]

                object@debugInfo <<- simDebug(
                    object,
                    simulation,
                    method = sprintf("simData[group=%d]", group),
                    simData
                )

                # If the object is flagged as pregenerated it should already have the
                # replicates (i.e. methyl seek), so skip this step.
                if (! object@pregenerated) {
                    simData <-
                        mapply(
                            makeReplicates,
                            as.data.frame(simData, stringsAsFactors = FALSE),
                            group,
                            simulation@times,
                            SIMPLIFY = FALSE,
                            USE.NAMES = FALSE
                        )
                }

                return(simData)
            },
            # Count column for each group in the correct order
            data.frame(object@data[as.character(simProfiles$ID), , drop = FALSE], stringsAsFactors = FALSE),
            # Pass only "group" columns
            columnsParam, #dplyr::select(simProfiles, starts_with("Group")),
            # Keep track of iteration
            iterationParam, #seq(simulation@numberGroups),
            # Pass the row IDs
            idsParam,
            SIMPLIFY = FALSE
        ), row.names = simProfiles$ID, stringsAsFactors = FALSE)

    # Ensure that the columns have the proper pattern "Counts.GroupX"
    colsPerGroup <- length(simulation@times) * simulation@numberReps

    colnames(object@simData) <- paste0("Counts.Group.",
                                       rep(seq(simulation@numberGroups), each = colsPerGroup),
                                       ".", seq(colsPerGroup))

    # Final modifications (if required)
    object <- postSimulation(object, simulation)

    return(object)
})

# .Simulator.adjustDepth <- function(object, data) {
#     if (object@depthAdjust) {
#         message("Adjusting to sequencing depth ", object@depth, " on ", object@name)
#
#         data <- apply(data, 2, function(col) round(object@depth * col/sum(col), object@depthRound))
#     }
#
#     data
# }

#' postSimulation
#'
#' Method to make final modifications required by the different omics. Like rounding
#' the final count values or renaming columns for making them more human-readable.
#'
#' For internal use only.
#'
#' @param object Instance of class \linkS4class{Simulator}.
#' @param simulation Instance of class \linkS4class{Simulation}.
#'
#' @return Object of class \linkS4class{Simulation} with the simulated data (@simData)
#' correctly formatted.
#'
setMethod("postSimulation", signature="Simulator", function(object, simulation) {

    if (ncol(object@simData) < simulation@numberGroups * simulation@numberReps * length(simulation@times))
        stop("Invalid number of columns after simulation. Please, contact package maintainer!")

    # Change number of columns with pattern "Counts.Group"
    # GroupX.TimeY.RepN
    colnames(object@simData)[grep("Counts.Group", colnames(object@simData))] <-
        paste0(
            "Group",
            rep(
                1:simulation@numberGroups,
                each = length(simulation@times) * simulation@numberReps
            ),
            ".Time",
            rep(simulation@times,
                each = simulation@numberReps),
            ".Rep",
            1:simulation@numberReps
        )

    # simData[simData < object@minValue] <- object@minValue
    #
    # if (length(object@maxValue)) {
    #     simData[simData > object@maxValue] <- object@maxValue
    # }

    # object@simData <- .Simulator.adjustDepth(object, object@simData)

    message("Rounding ", object@name, " count values.")

    # Round only counts columns
    roundCols <- grep(".Time", colnames(object@simData))

    object@simData[, roundCols] <- round(object@simData[, roundCols], digits = object@roundDigits)

    return(object)
})


#' IDfromGenes
#'
#' Returns the regulator IDs associated to a particular set of genes identifiers.
#'
#' @param object Instance of class \linkS4class{Simulator}.
#' @param geneNames Character vector of the genes to look for in the association table.
#' @param simplify Return only the genes IDs or a table containing both types of identifiers.
#'
#' @return Depending on \emph{simplify} parameter it would be a character vector or a
#' data frame containing columns \emph{ID} and \emph{Gene}.
#' @export
#'
setMethod("IDfromGenes", signature="Simulator", function(object, geneNames, simplify = TRUE) {
    regTable <- object@idToGene

    selCols <- if (simplify) c('ID') else c('ID', 'Gene')
    selRows <- regTable[, 'Gene'] %in% geneNames

    return(regTable[selRows, selCols])
})

#' IDtoGenes
#'
#' Returns the gene identifiers associated to a particular set of regulator IDs.
#'
#' @param object Instance of class \linkS4class{Simulator}.
#' @param idNames Character vector of the regulator IDs to look for in the association table.
#' @param simplify Return only the regulator IDs or a table containing both types of identifiers.
#'
#' @return Depending on \emph{simplify} parameter it would be a character vector or a
#' data frame containing columns \emph{ID} and \emph{Gene}.
#' @export
#'
setMethod("IDtoGenes", signature="Simulator", function(object, idNames, simplify = TRUE) {
    regTable <- object@idToGene

    selCols <- if (simplify) c('ID') else c('ID', 'Gene')
    selRows <- regTable[, 'ID'] %in% idNames

    return(regTable[selRows, selCols])
})

#' simulateParams
#'
#' For internal use. Creates the values to be replaced in the formulas used
#' to simulate the profile values, including noise.
#'
#' @param object Instance of class \linkS4class{Simulator}.
#' @param simulation Instance of class \linkS4class{Simulation}.
#' @param counts Counts taken from the initial data.
#' @param profiles Types of profiles to simulate.
#' @param ids IDs associated. Needed in some simulators.
#'
#' @return A list containing the following items:
#' \describe{
#'      \item{randomCounts}{numeric vector containing random count values.}
#'      \item{noiseValues}{numeric vector containing noise values generated with the noise function and parameters specified.}
#'      \item{m}{numeric vector of lower values comparing original and random counts.}
#'      \item{M}{numeric vector of maximum values comparing original and random counts.}
#' }
#' @export
#'
setMethod("simulateParams", signature="Simulator", function(object, simulation, counts, profiles, ids) {
    # Generate random counts
    randomCounts <-
        runif(length(counts), min = object@min, max = object@max)

    # Order [m, M]
    m <- pmin(counts, randomCounts)
    M <- pmax(counts, randomCounts)

    # Generate noise values. Different values for each time.
    noiseValues <-
        replicate(length(simulation@times),
                  do.call(
                      object@noiseFunction,
                      append(list("n" = length(counts)),
                             object@noiseParams[names(object@noiseParams) != 'NB'])
                  ))

    return(list(
        'randomCounts' = randomCounts,
        'noiseValues' = noiseValues,
        'm' = m,
        'M' = M
    ))
})

#' adjustProfiles
#'
#' For internal use. Allows every omic class to adjust the profiles matrix
#' (simulation settings) to use, according to its needs.
#'
#' @param object Instance of class \linkS4class{Simulator}.
#' @param simulation Instance of class \linkS4class{Simulation}.
#' @param profiles Data frame containing the profile type associated to each
#' ID.
#' @param step It accepts two options \emph{Effect} and \emph{Groups} depending
#' on the part of the process where this method is called. Check 'Simulation.R'.
#'
#' @return By default it returs the profiles data frame without modifications.
#'
setMethod("adjustProfiles", signature="Simulator", function(object, simulation, profiles, step) {
    # By default return the same
    return(profiles)
})

setMethod("show", signature="Simulator", function(object) {
    print(object@simData)
})

setMethod("simSettings", signature="Simulator", function(object) {
    cat(sprintf("Simulation settings of class %s:\n", class(object)))
    cat(sprintf("- Depth: %d\n", object@depth))
})

#' simDebug
#'
#' Helper function to trace how the different values are modified during
#' the execution of the algorithm. It is only intended for testing purposes
#' and as a final user it should not be used.
#'
#' @param object Instance of class \linkS4class{Simulator}.
#' @param simulation Instance of class \linkS4class{Simulation}.
#' @param method Name of the method calling the function.
#' @param ... Values to keep track of.
setMethod("simDebug", signature="Simulator", function(object, simulation, method, ...) {

    if (simulation@debug) {
        values <- list(...)

        valuesNames <-
            stri_trim(unlist(strsplit(
                stri_match_first(deparse(substitute(list(...)), width.cutoff = 500),
                                 regex = "list\\((.*)\\)")[, 2], ","))
                )

        object@debugInfo[[method]] <-
            append(object@debugInfo[[method]],
                   setNames(
                       values,
                       if (is.null(names(values)))
                           valuesNames
                       else
                           ifelse(names(values) == "", valuesNames, names(values))
                   ))
    }

    return(object@debugInfo)
})


setValidity("Simulator", function(object) {
    errors <- c()

    if (! is.declared(object@data) && ! object@pregenerated)
        errors <- c(errors, "Every omic needs to have the initial data set.")

    if (object@regulator && ! is.declared(object@idToGene))
        errors <- c(errors, "Regulators must provide the association gene list.")

    return(if(length(errors)) errors else TRUE)
})
