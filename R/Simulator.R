#' @include AllGeneric.R
#' @importFrom rlang .data
NULL

#' @rdname Generics
setMethod("initialize", signature="Simulator", function(.Object, ...) {
    .Object <- callNextMethod(.Object, ...)

    # Casting to dataframe
    .Object@idToGene <- as.data.frame(.Object@idToGene, stringsAsFactors = FALSE)
    .Object@data <- as.data.frame(.Object@data, stringsAsFactors = FALSE)

    # Convert data frames to characters or numbers
    .Object@idToGene[] <- lapply(.Object@idToGene, as.character)
    .Object@data[] <- lapply(.Object@data, as.numeric)

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
#' @keywords internal
#' @rdname Generics
#'
setMethod("initializeData", signature="Simulator", function(object, simulation) {

    # Adjust sequencing depth
    object@depth <- object@depth*10^6

    if (! is.declared(object@data)) {
        stop(sprintf("Data not declared on simulator %s.", object@name))
    }

    # Adjust counts to depth
    object@data$Counts <- object@depth * object@data$Counts / sum(object@data$Counts)

    # Set min/max values based on provided data
    object@min <- min(object@data)
    object@max <- max(object@data)

    nonzero.org.counts <- object@data[object@data > 0]
    no.outliers.counts <- object@data[object@data < stats::quantile(object@data$Counts, c(0.99))]
    min.value.nonzero <- stats::quantile(nonzero.org.counts, simulation@minMaxQuantile[1])

    dataPerc <- stats::quantile(nonzero.org.counts, object@minMaxQuantile)

    range.no.outliers <- range(no.outliers.counts)
    minDE.value <- stats::quantile(nonzero.org.counts, c(simulation@minMaxQuantile[1]))
    maxDE.value <- stats::quantile(nonzero.org.counts, c(simulation@minMaxQuantile[2]))
    maxDE.value <- max(no.outliers.counts)


    # TODO: this might fail when data does have more than one column (pregenerated?)
    object@minMaxDist <- list(
        "min" = dataPerc[1],
        "max" = dataPerc[2],
        "minDE" = minDE.value,
        "maxDE" = maxDE.value,
        "minValue" = min.value.nonzero
    )

    # Generate a random counts series
    omicSettings <- simulation@simSettings$geneProfiles[[class(object)]]

    featureSamples <- simulation@simSettings$featureSamples[[class(object)]]

    featuresDE <- featureSamples$DE

    features.Flat <- simulation@simSettings$geneProfiles$FlatGroups[[class(object)]]$ID

    featuresDE.notFlat <- setdiff(featuresDE, features.Flat)

    featuresNONDE <- featureSamples$nonDE

    featuresNonDE.noised <- featureSamples$noiseNonDE

    featuresNONDE.diff <- setdiff(featuresNONDE, featuresNonDE.noised)

    object@randData <- sample(object@data$Counts)
    names(object@randData) <- rownames(object@data)

    # Adjust randomCounts to depth
    # object@randData <- object@depth * object@randData / sum(object@randData)

    swapCountValues <- function(countVector, idsToChange, idsNonDE, meanValue, operation = NULL) {
        # Select non-DE IDs to swap. Select the maximum amount available
        swapID <- sample(as.character(idsNonDE),
                         size = min(length(idsToChange), length(idsNonDE)),
                         replace = FALSE)

        # Take a sample to replace
        swapID.bad <- sample(as.character(idsToChange),
                             size = length(swapID),
                             replace = FALSE)

        # Good values to use
        swapValues <- countVector[swapID]

        # Set bad values into the Non-DE features
        countVector[swapID] <- countVector[swapID.bad]

        # Assign good values to DE features
        countVector[swapID.bad] <- swapValues

        # In case any value has not been changed due to
        # not having enough replacement values, use random distribution.
        remainingIds <- setdiff(idsToChange, swapID.bad)

        if (length(remainingIds)) {
            countVector[remainingIds] <- stats::rnorm(n = length(remainingIds),
                                                             mean = meanValue,
                                                             sd = 1)
        }

        return(countVector)
    }

    swapCountValuesDiff <- function(countVector, initialVector, idsToChange, nonDEIds, limitMax, limitMin, limit, absMaxLimit, absMinLimit) {

        if (length(idsToChange)) {
            minorCounts <- setNames(pmin(initialVector[idsToChange], countVector[idsToChange]),
                                    idsToChange)

            majorCounts <- setNames(pmax(initialVector[idsToChange], countVector[idsToChange]),
                                    idsToChange)

            diffMask <- (countVector[idsToChange] > initialVector[idsToChange])

            diffCounts <- setNames(majorCounts - minorCounts, idsToChange)

            limitCriteria.min <- setNames(minorCounts * limitMin, idsToChange)
            limitCriteria.max <- setNames(minorCounts * limitMax, idsToChange)

            vectorToChange <- ifelse(diffMask, 'countVector', 'initialVector')

            limitSelValue <-  if (limit == 'max') limitMax else limitMin

            diffToLimit <- setNames(minorCounts * limitSelValue - majorCounts, idsToChange)

            # Order the ids for the more restrictive (lower range) to the most.
            idsToChange <- idsToChange[order(limitCriteria.max - limitCriteria.min)]

            # Repeat the process for each feature, as the pool will change with every situation.
            for (featureID in idsToChange) {
                minValue <- limitCriteria.min[featureID]
                maxValue <- limitCriteria.max[featureID]

                vectorVariable <- get(vectorToChange[featureID])

                matchIds <- names(vectorVariable)[vectorVariable >= minValue & vectorVariable <= maxValue]

                availableIds <- intersect(matchIds, nonDEIds)

                if (length(availableIds)) {

                    swapID <- sample(availableIds, size = 1)

                    swap.Value <- vectorVariable[swapID]
                    feature.Value <- vectorVariable[featureID]

                    vectorVariable[featureID] <- swap.Value
                    vectorVariable[swapID] <- feature.Value
                } else {
                    # If no match has been found, generate a new random value.
                    newValue <- min(max(stats::runif(1,
                                                     min = minValue,
                                                     max = maxValue), absMinLimit), absMaxLimit)

                    vectorVariable[featureID] <- newValue
                }

                # As no pointers exist in R, replace original vector value with the new modified one.
                assign(vectorToChange[featureID], vectorVariable)
            }
        }

        outputList <- list(
            'data' = initialVector,
            'randData' = countVector
        )

        return(outputList)
    }

    # Determine which DE non-flat genes have been assigned a non-expressed status AND have an original
    # value lower than that, then swap their values with non-DE genes
    minDEassignments.rand <- names(object@randData)[object@randData < minDE.value] # Keep IDS from randomCounts separate
    minDEassignments.org <- rownames(object@data)[object@data$Counts < minDE.value]


    badAssignmentDE.min.org <- intersect(featuresDE, minDEassignments.org)


    availableIdsNONDE.rand <- setdiff(featuresNONDE, minDEassignments.rand)
    availableIdsNONDE.min.org <- setdiff(featuresNONDE.diff, minDEassignments.org)


    # First check: make sure they have a minimum value
    object@data[] <- swapCountValues(setNames(object@data$Counts, rownames(object@data)),
                                     badAssignmentDE.min.org,
                                     availableIdsNONDE.min.org,
                                     minDE.value)

    maxDEassignments.org <- rownames(object@data)[object@data$Counts > maxDE.value]
    badAssignmentDE.max.org <- intersect(featuresDE, maxDEassignments.org)

    availableIds <- names(object@data)[object@data < maxDE.value & object@data > minDE.value]
    availableIdsNONDE.max.org <- intersect(featuresNONDE.diff, availableIds)

    object@data[] <- swapCountValues(setNames(object@data$Counts, rownames(object@data)),
                                     badAssignmentDE.max.org,
                                     availableIdsNONDE.max.org,
                                     maxDE.value)


    # TODO: remove this temp mod
    index.DE <- names(object@randData) %in% featuresDE.notFlat
    index.noised <- names(object@randData) %in% featuresNonDE.noised

    minNoiseValue <- 100

    noisedMin.org <- rownames(object@data)[object@data$Counts < minNoiseValue]
    noisedMin.rand <- names(object@randData)[object@randData < minNoiseValue]

    badAssignmentNoiseNonDE.min.org <- intersect(featuresNonDE.noised, noisedMin.org)
    badAssignmentNoiseNonDE.min.rand <- intersect(featuresNonDE.noised, noisedMin.rand)

    availableIdsNONDE.min.org <- setdiff(featuresNONDE.diff, noisedMin.org)
    availableIdsNONDE.min.rand <- setdiff(featuresNONDE.diff, noisedMin.rand)

    object@data[] <- swapCountValues(setNames(object@data$Counts, rownames(object@data)),
                                     badAssignmentNoiseNonDE.min.org,
                                     availableIdsNONDE.min.org,
                                     minNoiseValue)

    object@randData[] <- swapCountValues(object@randData,
                                         badAssignmentNoiseNonDE.min.rand,
                                         availableIdsNONDE.min.rand,
                                         minNoiseValue)

    minDEassignments.rand <- names(object@randData)[object@randData < minDE.value]
    badAssignmentDE.min.rand <- intersect(featuresDE, minDEassignments.rand)
    badAssignmentNoiseNonDE.min.rand <- intersect(featuresNonDE.noised, minDEassignments.rand)
    availableIdsNONDE.min.rand <- setdiff(featuresNONDE.diff, minDEassignments.rand)

    object@randData[] <- swapCountValues(object@randData,
                                         badAssignmentDE.min.rand,
                                         availableIdsNONDE.min.rand,
                                         minDE.value)

    maxDEassignments.rand <- names(object@randData)[object@randData > maxDE.value]
    badAssignmentDE.max.rand <- intersect(featuresDE, maxDEassignments.rand)

    availableIds <- names(object@randData)[object@randData < maxDE.value & object@randData > minDE.value]
    availableIdsNONDE.max.rand <- intersect(featuresNONDE.diff, availableIds)

    object@randData[] <- swapCountValues(object@randData,
                                         badAssignmentDE.max.rand,
                                         availableIdsNONDE.max.rand,
                                         maxDE.value)

    retrieveDiffIndex <- function(instance, multiplier) {
        random.M <- pmax(instance@randData, instance@data$Counts)
        random.m <- pmin(instance@randData, instance@data$Counts)

        diff.out <- (random.M > random.m * multiplier)

        return(diff.out)
    }

    # Maximum difference between M and m
    min.diff <- object@minMaxFC[1]
    max.diff <- object@minMaxFC[2]

    index.match.maxdiff.available <- retrieveDiffIndex(object, max.diff)
    index.match.maxdiff <- (index.match.maxdiff.available & index.DE)

    # Swap values with others to meet the constraints at the time we maintain a similar distribution
    # of the initial sample.

    # Ids to change
    badAssignmentDE.rand.maxdiff <- names(object@randData)[index.match.maxdiff]

    # Pool of available ids to meet the condition.
    availableIdsNONDE.maxdiff.rand <- setdiff(featuresNONDE.diff, names(index.match.maxdiff.available))

    correctedCounts <- swapCountValuesDiff(object@randData,
                                           setNames(object@data$Counts, rownames(object@data)),
                                           badAssignmentDE.rand.maxdiff,
                                           featuresNONDE.diff,
                                           max.diff,
                                           min.diff,
                                           'max',
                                           maxDE.value,
                                           minDE.value
                                           )

    object@data[] <- correctedCounts$data
    object@randData <- correctedCounts$randData

    # Repeat the process with minimum diffs
    index.match.mindiff.available <- ! (retrieveDiffIndex(object, min.diff))
    index.min.match <- (index.match.mindiff.available & index.DE)

    badAssignmentDE.rand.mindiff <- names(object@randData)[index.min.match]

    availableIdsNONDE.mindiff.rand <- setdiff(featuresNONDE.diff, names(index.match.mindiff.available))

    correctedCounts <- swapCountValuesDiff(object@randData,
                                           setNames(object@data$Counts, rownames(object@data)),
                                           badAssignmentDE.rand.mindiff,
                                           featuresNONDE.diff,
                                           max.diff,
                                           min.diff,
                                           'min',
                                           maxDE.value,
                                           minDE.value)

    object@data[] <- correctedCounts$data
    object@randData <- correctedCounts$randData

    # MIN/MAX FOR NOISED NON-DEG
    max.diff.noised <- 1.1
    min.diff.noised <- 1.1

    index.match.maxdiff.available.noised <- retrieveDiffIndex(object, max.diff.noised)
    index.max.match.noised <- (index.match.maxdiff.available.noised & index.noised)

    badAssignmentDE.rand.maxdiff.noised <- names(object@randData)[index.max.match.noised]

    availableIdsNONDE.maxdiff.rand.noised <- setdiff(featuresNONDE.diff, names(badAssignmentDE.rand.maxdiff.noised))

    correctedCounts <- swapCountValuesDiff(object@randData,
                                           setNames(object@data$Counts, rownames(object@data)),
                                           badAssignmentDE.rand.maxdiff.noised,
                                           setdiff(featuresNONDE.diff, featuresNonDE.noised),
                                           max.diff.noised,
                                           min.diff.noised,
                                           'max',
                                           maxDE.value,
                                           minNoiseValue)

    object@data[] <- correctedCounts$data
    object@randData <- correctedCounts$randData

    # Clone the base counts for every group
    if (ncol(object@data) < simulation@numberGroups) {
        object@data <- object@data[, rep.int(1, simulation@numberGroups), drop = FALSE]
    }

    # Important: column names pattern "Counts.Group" is required .
    colnames(object@data) <- paste('Counts.Group', seq(simulation@numberGroups))

    if (! is.null(flatProfiles <- simulation@simSettings$geneProfiles$FlatGroups$SimRNAseq)) {

        # Structure: Gene_ID | Group_1 | ... | Group_N
        # onlyFlat <- (length(simulation@times) > 1)

        # Count range to exclude comparing with the reference (group 1)
        excludeMargin <- 10

        # Match positions
        flatIndexes <- match(flatProfiles$ID, rownames(object@data))

        # Different treatments between RNA-seq and the rest of simulators
        if ( ! object@regulator) {

            # Profiles + Counts
            profileCounts <- cbind(flatProfiles, object@data[flatIndexes, ])

            profileColumns <- dplyr::starts_with("Group", vars = colnames(profileCounts))
            countColumns <- dplyr::starts_with("Counts.Group", vars = colnames(profileCounts))

            # Iterate over the rows on a unique matrix with the structure:
            # Gene_ID | Profile_Group_1 | ... | Counts_Group_1 | ...
            object@data[flatIndexes, ] <- t(apply(profileCounts, 1, function(idRow) {
                # Types of profiles (exclude first column)
                idProfiles <- as.character(idRow[profileColumns])

                # Count values (columns after profiles)
                idCounts <- as.numeric(idRow[countColumns])

                min.flat.value.enhancer <- idCounts[1] * object@minMaxFC[1]
                max.flat.value.enhancer <- idCounts[1] * object@minMaxFC[2]

                min.flat.value.repressor <- idCounts[1] / object@minMaxFC[2]
                max.flat.value.repressor <- idCounts[1] / object@minMaxFC[1]

                # Profiles
                tmpCounts <- ifelse(idProfiles == 'flat',
                                    # If the profile is flat, keep the current value
                                    idCounts,
                                    ifelse(idProfiles == 'enhancer',
                                           # If the profile is enhancer, generate a new count value in the range
                                           # (actual_value + excludeMargin) to maxValue
                                           rep.int(stats::runif(1,
                                                                min = min(min.flat.value.enhancer, object@max),
                                                                max = min(max.flat.value.enhancer, object@max)),
                                                   simulation@numberGroups),
                                           # If the profile is repressor, generate a new count value in the range
                                           # minValue to (actual_value - excludeMargin)
                                           rep.int(stats::runif(1,
                                                                min = max(min.flat.value.repressor, object@min),
                                                                max = max(max.flat.value.repressor, object@min)),
                                                   simulation@numberGroups)
                                    )
                )

                # Check the differences
                # referenceGroup <- head(tmpCounts[idProfiles == 'flat'], 1)
                # diffValues <- abs(tmpCounts - referenceGroup)
                #
                # # For enhancer there shouldn't be any problem
                # invalidDiffs <- (diffValues < min.diff.flat) & (idProfiles == 'repressor')
                #
                # if (any(invalidDiffs)) {
                #     # Calculate the amount to increase the reference group
                #     diffRange <- min.diff.flat - min(diffValues)
                #
                #     tmpCounts[idProfiles == 'flat'] <- tmpCounts[idProfiles == 'flat'] + diffRange
                # }

                return(tmpCounts)
            }))
        } else {
            # In case of regulators, the initial value will depend on the effect
            # on the associated gene.
            profileSubset <- dplyr::rename(flatProfiles, Gene = .data$ID) %>%
                dplyr::inner_join(simulation@simSettings$geneProfiles[[class(object)]][, c('ID', 'Effect', 'Gene')], by = c("Gene" = "Gene")) %>%
                dplyr::select(.data$ID, .data$Effect, dplyr::starts_with("Group")) %>%
                dplyr::filter(! is.na(.data$Effect))

            # TODO: change object@data with object@randData?
            if (nrow(profileSubset)) {
                # CAREFUL: the counts value NEED to have the same amount as columns as number of groups, otherwise a NA
                # will be returned hence assigning the maximum value.
                # Profiles + Counts
                profileCounts <- cbind(profileSubset, object@data[profileSubset$ID, ])

                effectColumn <- dplyr::matches("Effect", vars = colnames(profileCounts))
                profileColumns <- dplyr::starts_with("Group", vars = colnames(profileCounts))
                countColumns <- dplyr::starts_with("Counts.Group", vars = colnames(profileCounts))

                object@data[profileSubset$ID, ] <- t(apply(profileCounts, 1, function(idRow) {
                    # Regulator effect
                    regEffect <- idRow[effectColumn]
                    # Types of profiles (skip first 2 columns: ID and Effect)
                    idProfiles <- as.character(idRow[profileColumns])
                    # Count values (following columns)
                    idCounts <- as.numeric(idRow[countColumns])

                    min.flat.value.enhancer <- idCounts[1] * 3
                    max.flat.value.enhancer <- idCounts[1] * 10

                    min.flat.value.repressor <- idCounts[1] / 10
                    max.flat.value.repressor <- idCounts[1] / 3

                    # Profiles
                    # TODO: change to use randData
                    return(ifelse(idProfiles == 'flat',
                                  # If the profile is flat, keep the current value
                                  idCounts,
                                  # If the regulator effect is enhancer, keep the same behaviour,
                                  # if not, assign the opposite.
                                  ifelse(idProfiles == regEffect,
                                         # If the profile is enhancer, generate a new count value in the range
                                         # (actual_value + excludeMargin) to maxValue
                                         rep.int(stats::runif(1,
                                                              min = min(min.flat.value.enhancer, object@max),
                                                              max = min(max.flat.value.enhancer, object@max)),
                                                 simulation@numberGroups),
                                         # If the profile is repressor, generate a new count value in the range
                                         # minValue to (actual_value - excludeMargin)
                                         rep.int(stats::runif(1,
                                                              min = max(min.flat.value.repressor, object@min),
                                                              max = max(max.flat.value.repressor, object@min)),
                                                 simulation@numberGroups)
                                  )
                    ))
                }))
            }
        }

        # Adjust again to depth for each column
        object@data <- apply(object@data, 2, function(x) x * object@depth / sum(x))
        object@randData <-  object@randData * object@depth / sum(object@randData)

        # Replace with a mean value between the groups
        object@data[featuresNONDE, ] <- apply(object@data[featuresNONDE, , drop = FALSE], 1, mean)
    }

    return(object)
})

#' @rdname Generics
#' @keywords internal
setMethod("simulate", signature="Simulator", function(object, simulation) {
    message(sprintf("Starting simulation of %s.", object@name))

    object <- initializeData(object, simulation)

    makeReplicates <- function(counts, groupInfo, timeInfo, profileInfo, verbose = TRUE) {
        if (verbose) {
            message(sprintf("\t- Making replicates for group %d on time %s.", groupInfo, timeInfo))
        }

        counts[is.na(counts) | counts < 0.1] <- 0.1

        varest <- 10^object@replicateParams$a * (counts + 1)^object@replicateParams$b - 1

        # Force a minimum variance
        varest[varest < 0.03] <- 0.03

        # Apply to every row:
        #   Param x: [mu.noise.cpm, stdev]
        temp <- t(matrix(apply(cbind(counts, varest), 1, function (x) {
            if (x[1] >= x[2]) {
                replis <- stats::rpois(n = simulation@numberReps, lambda = x[1])
            } else {
                replis <- stats::rnbinom(n = simulation@numberReps, size = x[1]^2/(x[2] - x[1]),  mu = x[1])
            }

            return(replis)
        }), nrow = simulation@numberReps))

        flatProfiles <- which(profileInfo == 'flat')

        # Disabled: adjusting to specified depth
        # temp.depth <- apply(temp, 2, function(x) x * object@depth / sum(x))
        temp.depth <- temp

        # Generate random means
        noised.means.rand <- stats::rnorm(n = length(counts),
                             mean = counts,
                             sd = 0.025 * counts)

        # For induction profiles, select the minimum value
        noised.means <- ifelse(grepl("induction", profileInfo),
               pmin(counts, noised.means.rand),
               pmax(counts, noised.means.rand)
               )

        # Override for flat profiles with a less variable noise
        if (length(flatProfiles)) {
            noised.means[flatProfiles] <- stats::rnorm(n = length(flatProfiles),
                                                       mean = counts[flatProfiles],
                                                       sd = 0.01 * counts[flatProfiles])
        }

        # Note: assign as [] to maintain structure, as apply will return a matrix or a vector depending
        # on the number of replicates of the simulation.
        temp.depth[] <- t(apply(cbind(noised.means, temp.depth), 1,
                                              function(geneRow) {
                                                  # Original replicate values
                                                  mean.adjusted <- geneRow[1]
                                                  gene.values <- geneRow[-1]

                                                  # Mean
                                                  diff.mean <- mean(gene.values) - mean.adjusted

                                                  gene.adjusted <- gene.values - diff.mean

                                                  return(gene.adjusted)
                                              }))

        return(temp.depth)
    }

    # Profile table of simulator
    simProfiles <- simulation@simSettings$geneProfiles[[class(object)]]

    # Flat DE genes
    flatProfiles <- simulation@simSettings$geneProfiles$FlatGroups[[class(object)]] #$SimRNAseq

    # Select only the active rows of the regulator profile table, discarding
    # the duplicated rows (all active rows for a given regulator will have the
    # same profile among groups by now).
    if (object@regulator) {
        # If the object is pregenerated (e.g. methylation) then keep the profiles
        # as they are.
        if (object@pregenerated) {
            simProfiles <- dplyr::select(simProfiles, .data$ID, dplyr::starts_with("Group"), dplyr::starts_with("Tmax")) %>%
                dplyr::distinct_()
        } else {
            simProfiles <- dplyr::filter(simProfiles, ! is.na(.data$Effect)) %>%
                dplyr::select(.data$ID, dplyr::starts_with("Group"), dplyr::starts_with("Tmax")) %>%
                dplyr::distinct_()

            # Add remaining IDs present on data with a flat profile
            simProfiles <- rbind(simProfiles, do.call(cbind, setNames(append(
                # ID
                list(rownames(object@data)[! rownames(object@data) %in% as.character(simProfiles$ID)]),
                # Groups
                append(rep('flat', simulation@numberGroups), rep(NA, simulation@numberGroups)) # Group & Tmax columns
            ), colnames(simProfiles))))
        }
    }

    # Filtering step
    simProfiles <- dplyr::filter(simProfiles, .data$ID %in% rownames(object@data))

    if (object@pregenerated) {
        # TODO: fix this in other place. When working with blocks, sometimes there are duplicated rows
        # because they do not have the tmax columns properly filled.
        simProfiles <- dplyr::group_by(simProfiles, .data$ID) %>%
            dplyr::summarise_all(dplyr::funs(dplyr::first(stats::na.omit(.)))) %>%
            dplyr::ungroup()
    }

    # If data is already generated (i.e. methylation simulator) skip some
    # steps like generating random counts or replicates, adjusting also the
    # parameters for mapply.

    # Pass only "group" columns
    columnsParam <- dplyr::select(simProfiles, dplyr::starts_with("Group"))
    # Pass only "tmax" columns
    tmaxParam <- dplyr::select(simProfiles, dplyr::starts_with("Tmax."))
    # Keep track of iteration
    iterationParam <- seq(simulation@numberGroups)
    # Pass IDs
    idsParam <- simProfiles[, rep('ID', simulation@numberGroups), drop = FALSE]

    if (object@pregenerated) {
        # Repeat for every replicate
        columnsParam <- columnsParam[, rep(colnames(columnsParam), each = simulation@numberReps), drop = FALSE]
        iterationParam <- rep(iterationParam, each = simulation@numberReps)
        idsParam <- simProfiles[, rep('ID', simulation@numberGroups * simulation@numberReps), drop = FALSE]

        tmaxParam <- tmaxParam[, rep(colnames(tmaxParam), each = simulation@numberReps), drop = FALSE]
    }

    simulatedParams <- setNames(mapply(function(counts, profiles, group, ids) {
        simulateParams(object, simulation, counts, profiles, group, ids)
    },
    data.frame(object@data[as.character(simProfiles$ID), , drop = FALSE], stringsAsFactors = FALSE),
    columnsParam,
    iterationParam,
    idsParam,
    USE.NAMES = TRUE,
    SIMPLIFY = FALSE
    ), paste0("Group", seq_len(simulation@numberGroups)))

    features2Groups <- c()

    # Genes flat in one group but not the other
    if (simulation@numberGroups > 1) {
        # Check if at least one gene has a flat profile in one group
        features2Groups <- dplyr::select(simProfiles, dplyr::starts_with("Group")) %>%
            purrr::map(`==`, 'flat') %>%
            purrr::pmap_int(any)

        # Exclude those with 'flat' in all profiles
        flatIndexes <- simProfiles$ID %in% flatProfiles$ID

        features2Groups <- simProfiles[features2Groups & ! flatIndexes, 'ID', drop = TRUE]

        # dplyr::select(profilesDE, dplyr::starts_with("Group")) %>% purrr::pmap(., function(...) {
        #     geneProfiles <- list(...)
        #
        #     return(any(geneProfiles == 'flat'))
        # })
    }

    if (simulation@numberReps > 1) {
        nonDE.selected <- simulation@simSettings$featureSamples[[class(object)]]$noiseNonDE


        # Disabled: chose totally opposite profiles for noise genes.
        #
        # sampleProfiles <- function(featureID) {
        #     profileOptions <- names(simulation@profileProbs)
        #     profileTypes <- c("transitory", "continuous")
        #
        #     sampleTypes <- sample(profileTypes, size = floor(simulation@numberReps/2), replace = TRUE)
        #
        #     chosenOptions <- unlist(lapply(sampleTypes, function(typeName) profileOptions[grep(typeName, profileOptions)]))
        #
        #     # For an odd number of replicates, add a final flat profile
        #     if (simulation@numberReps %% 2 != 0) {
        #         chosenOptions <- c(chosenOptions, "flat")
        #     }
        #
        #     return(chosenOptions)
        # }
        profileOptions <- names(simulation@profileProbs)[-grep("flat", names(simulation@profileProbs))]

        replacement <- (simulation@numberReps > length(profileOptions))

        repProfiles <- do.call(rbind,
                               replicate(n = length(nonDE.selected),
                                         sample(profileOptions, size = simulation@numberReps, replace = replacement),
                                         simplify = FALSE))

        # repProfiles <- do.call(rbind, lapply(nonDE.selected, sampleProfiles))

        tmax.T <- length(simulation@times) - 1

        repTmax <- matrix(rep(stats::runif(length(nonDE.selected), min = 0.25 * tmax.T, max = 0.75 * tmax.T),
                            times = simulation@numberReps,
                            simplify = FALSE), ncol = simulation@numberReps)
    }

    # Generate time series (if any) with replicates
    object@simData <- data.frame(
        mapply(
            function(simulateParams, profiles, tmaxValues, group, ids) {
                if (! object@pregenerated) {
                    message(sprintf("- Simulating count values for group %d.", group))
                } else {
                    message(sprintf("[Pregenerated replicate] Simulating count values for group %d.", group))
                }

                # Create the coefficients based on the "tmax" parameter
                nt <- length(simulation@times)

                x <- c(0:(nt - 1))

                varT <- x[nt]

                calculateTimeVectors <- function(localProfiles, localTmaxvalues) {
                    t(apply(cbind(localProfiles, localTmaxvalues), 1, function(rowInfo) {
                    # Set in the current scope the proper values
                    profileValue <- rowInfo[1]
                    tmaxValue <- as.numeric(rowInfo[2])

                    data_values <- list(
                        "0" = 0,
                        a1 = 0,
                        a2 = 0,
                        a2.neg = 0,
                        b1 = 1 / x[nt],
                        b1.neg = -1 / x[nt],
                        b2 = 4 * (x[nt]) / (x[nt] * x[nt]),
                        b2.neg = -4 * (x[nt]) / (x[nt] * x[nt]),
                        c2 = -4 / (x[nt] * x[nt]),
                        c2.neg = 4 / (x[nt] * x[nt])
                    )

                    # If tmaxValue is set, the profile is transitory: override variables
                    if (! is.na(tmaxValue)) {

                        if (tmaxValue >= varT/2) {
                            data_values$b2 <- 2 / tmaxValue
                            data_values$c2 <- - 1 /tmaxValue^2
                        } else {
                            data_values$c2 <- - 1 /(varT - tmaxValue)^2
                            data_values$a2 <- 1 + (data_values$c2) * tmaxValue^2
                            data_values$b2 <- - 2 * (data_values$c2) * tmaxValue
                        }

                        data_values$a2.neg <- - data_values$a2
                        data_values$b2.neg <- - data_values$b2
                        data_values$c2.neg <- - data_values$c2
                    }

                    return(unlist(data_values[as.character(simulation@profiles[[profileValue]])]))

                })) %*% rbind(c(rep(1, nt)), x, x * x)
                }

                generateCountsProfiles <- function(timeVectors, profiles, indexID, applyCorrection=TRUE, verbose=TRUE) {

                    # Generate base counts using the formulas:
                    # X = initial_counts + noise [for flat]
                    # X = m + lambda(M - m) + noise [for lambda in [0,1]]
                    # X = M + lambda(M - m) + noise [for lambda in [-1,0]]

                    # Delegate creation of specific simulation parameters
                    counts <- simulateParams$counts[indexID]
                    randomCounts <- simulateParams$randomCounts[indexID]
                    noiseValues <- simulateParams$noiseValues[indexID]
                    M <- simulateParams$M[indexID]
                    m <- simulateParams$m[indexID]
                    repressionMean <- simulateParams$repression.Mean[indexID]
                    inductionMean <- simulateParams$induction.Mean[indexID]

                    # Rows with flat profile
                    indexFlat <- which(profiles == "flat")

                    diffRepression <- repressionMean - m
                    diffInduction <- M - inductionMean

                    if (! object@pregenerated) {
                        simData <-
                            ifelse(grepl('repression', profiles, fixed = TRUE),
                                   repressionMean,
                                   inductionMean) +
                            # ((M - m) * timeVectors) +
                            ifelse(grepl('repression', profiles, fixed = TRUE),
                                   diffRepression,
                                   diffInduction) * timeVectors +
                            noiseValues
                    } else {
                        # For methylation, avoid using means so that we can replicate the plot showing
                        # a maximum in 0 and 1 values
                        simData <-
                            ifelse(grepl('repression', profiles, fixed = TRUE), M, m) +
                            ((M - m) * timeVectors) +
                            noiseValues
                    }

                    # Overwrite flat rows with the correct value
                    otherGroup.mean <- purrr::pmap_dbl(lapply(grep(paste0("Group", group), names(simulatedParams), invert = TRUE), function(otherGroup) {
                        simulatedParams[[otherGroup]]$mean.counts[indexID]
                    }), mean)

                    if (length(indexFlat)) {
                        simData[indexFlat, ] <- ifelse(ids[indexID][indexFlat] %in% flatProfiles$ID,
                                                       counts[indexFlat],
                                                       # Reduce differentes between groups
                                                       ifelse(ids[indexID][indexFlat] %in% features2Groups,
                                                              otherGroup.mean[indexFlat],
                                                              randomCounts[indexFlat]))
                        + noiseValues[indexFlat]
                    }

                    # Some genes might have negative values due to the way the time profile formula is applied.
                    # For DEG, we correct these as to manually increase the amount in all time profiles to avoid this, increasing it
                    # to a minimum value.
                    if (applyCorrection && ! object@pregenerated) {
                        minimumDEGvalue <- object@minMaxDist$minDE

                        # Apply the correction to transitory profiles only
                        rowsToCheck <- grep("transitory", profiles)

                        if (length(rowsToCheck)) {
                            simData[rowsToCheck, ] <- t(apply(simData[rowsToCheck, , drop = FALSE], 1, function(rowValues) {

                                if (min(rowValues) < 0 ) {
                                    increaseAmount <- minimumDEGvalue - min(rowValues)

                                    correctedValues <- rowValues + increaseAmount

                                    return(correctedValues)
                                }

                                return(rowValues)
                            }))
                        }
                    }

                    # If the object is flagged as pregenerated it should already have the
                    # replicates (i.e. methyl seek), so skip this step.
                    if (! object@pregenerated) {
                        simData <-
                            mapply(
                                makeReplicates,
                                counts = as.data.frame(simData, stringsAsFactors = FALSE),
                                groupInfo = group,
                                timeInfo = simulation@times,
                                profileInfo = as.data.frame(profiles)[, rep(1, length(simulation@times)),drop=FALSE],
                                verbose = verbose,
                                SIMPLIFY = FALSE,
                                USE.NAMES = FALSE
                            )
                    }

                    return(simData)
                }

                globalTimeVectors <- calculateTimeVectors(profiles, tmaxValues)

                simData <- generateCountsProfiles(globalTimeVectors, profiles, seq_along(ids))

                if (simulation@numberReps > 1 && ! object@pregenerated) {
                    ids.positions <- match(nonDE.selected, ids)

                    for (repNum in seq(from = 1, to = simulation@numberReps)) {
                        customTimeVectors <- calculateTimeVectors(repProfiles[, repNum],
                                                                  repTmax[, repNum])

                        customSimdata <- generateCountsProfiles(customTimeVectors,
                                                                repProfiles[, repNum],
                                                                ids.positions,
                                                                applyCorrection = FALSE,
                                                                verbose = FALSE)

                        # Replace replicates in already simulated data
                        for (timePoint in seq_len(ncol(customTimeVectors))) {
                            simData[[timePoint]][ids.positions, repNum] <- customSimdata[[timePoint]][, repNum]
                        }
                    }
                }

                return(simData)
            },
            # Count column for each group in the correct order
            simulatedParams,
            # Pass only "group" columns
            columnsParam,
            # Coefficients parameters
            tmaxParam,
            # Keep track of iteration
            iterationParam,
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
#' @keywords internal
#' @rdname Generics
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
                seq_len(simulation@numberGroups),
                each = length(simulation@times) * simulation@numberReps
            ),
            ".Time",
            rep(simulation@times,
                each = simulation@numberReps),
            ".Rep",
            seq_len(simulation@numberReps)
        )

    if (! object@pregenerated) {
        object@simData[object@simData < object@min] <- object@min
        object@simData[object@simData > object@max] <- object@max
    }

    # Before adjusting to depth, we include a "dummy" feature that would adjust
    # the sum of each column to be the same, so as to maintain the profiles on
    # the genes that matter.
    diffToMax <- max(colSums(object@simData)) - colSums(object@simData)

    object@simData <- data.frame(mapply(function(simuCol, diffCol) {

        # First rescale the column
        rescaledColDifferences <- simuCol * diffCol / sum(simuCol)

        rescaledCol <- simuCol + rescaledColDifferences

        return(rescaledCol)
    }, as.data.frame(object@simData), diffToMax, SIMPLIFY = FALSE),
    row.names = rownames(object@simData),
    stringsAsFactors = FALSE)

    # Adjust to simulated depth
    object@simData <- apply(object@simData, 2, function(x) x * object@depth / sum(x))

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
#' @keywords internal
#' @rdname Generics
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
#' @keywords internal
#' @rdname Generics
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
#' @keywords internal
#' @rdname Generics
#'
setMethod("simulateParams", signature="Simulator", function(object, simulation, counts, profiles, group, ids) {
    # Reorder random counts
    randomCounts <- object@randData[ids]

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

    repressionMean <- stats::runif(length(M), min = (M+m)/2, max = M)
    inductionMean <- stats::runif(length(m), min = m, max = (M+m)/2)

    return(list(
        'randomCounts' = randomCounts,
        'noiseValues' = noiseValues,
        'm' = m,
        'M' = M,
        'repression.Mean' = repressionMean,
        'induction.Mean' = inductionMean,
        'mean.counts' = ifelse(grepl('repression', profiles),(repressionMean + m)/2, (inductionMean + M)/2),
        'counts' = counts
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
#' @keywords internal
#' @rdname Generics
#'
setMethod("adjustProfiles", signature="Simulator", function(object, simulation, profiles, step) {
    # By default return the same
    return(profiles)
})

#' @rdname Generics
setMethod("show", signature="Simulator", function(object) {
    print(object@simData)
})

#' @rdname Generics
setMethod("simSettings", signature="Simulator", function(object) {
    cat(sprintf("Simulation settings of class %s:\n", class(object)))
    cat(sprintf("- Depth: %d\n", object@depth))
})


setValidity("Simulator", function(object) {
    errors <- c()

    if (! is.declared(object@data) && ! object@pregenerated)
        errors <- c(errors, "Every omic needs to have the initial data set.")

    if (object@regulator && ! is.declared(object@idToGene))
        errors <- c(errors, "Regulators must provide the association gene list.")

    return(if(length(errors)) errors else TRUE)
})
