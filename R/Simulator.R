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

    # Remove those rownames not present on idToGene association list
    # if (.Object@regulator) {
    #     .Object@data <- .Object@data[rownames(.Object@data) %in% .Object@idToGene$ID, ]
    # }

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
    if (! is.declared(object@data)) {
        stop(sprintf("Data not declared on simulator %s.", object@name))
    }

    # Adjust counts to depth
    object@data$Counts <- object@depth * object@data$Counts / sum(object@data$Counts)

    object@min <- min(object@data)
    object@max <- max(object@data)

    # if (any(object@minMaxQuantile > 0)) {
    #     dataPerc <- quantile(object@data$Counts, object@minQuantile)
    #
    #     # Minimum value to change when simulating data
    #     # TODO: this does not use indexes anymore but numeric positions instead
    #     object@increment <- (dataPerc[2] - dataPerc[1]) * 0.1
    # }
    dataPerc <- quantile(object@data$Counts, object@minMaxQuantile)

    # TODO: for flat genes, override
    dataPerc <- quantile(object@data[object@data > 0], object@minMaxQuantile)

    # TODO: this might fail when data does have more than one column (pregenerated?)
    object@minMaxDist <- list(
        "min" = dataPerc[1],
        "max" = dataPerc[2]
    )

    # Generate a random counts series


    # Set to 0 the same amount as in the original data
    # excluding the DEG
    # nonExprNumber <- object@data[object@data < 5]
    #
    # indexSample <- setdiff(
    #     seq_along(object@data),
    #     match(simulation@simSettings$geneSamples$DE, rownames(object@data))
    # )
    omicSettings <- simulation@simSettings$geneProfiles[[class(object)]]

    if (object@regulator) {
        featuresDE <- unique(dplyr::filter(omicSettings, ! is.na(Effect))$ID)
    } else {
        featuresDE <- unique(dplyr::filter(omicSettings, DE == TRUE)$ID)
    }

    features.Flat <- simulation@simSettings$geneProfiles$FlatGroups[[class(object)]]$ID

    featuresDE.notFlat <- setdiff(featuresDE, features.Flat)

    featuresNONDE <- setdiff(rownames(object@data), featuresDE)

    # geneDE.index <- match(geneDE.notFlat, rownames(object@data))
    # geneNONDE.index <- match(geneNONDE, rownames(object@data))

    # TODO: if we keep the "min.diff" (see below) this code might be redundant.
    # Remove it or leave it (as it swaps positions with DE so it may keep the "density" more
    # similar to the original data)
    minDE.value <- quantile(object@data[object@data > 0], c(0.3))

    #
    # object@randData[sample(indexSample, min(nonExprNumber, length(indexSample)))] <- 0

    # # TODO: remove this
    # object@randData <- object@data[sample(rownames(object@data)),]

    kdeOrgData <- density(object@data$Counts, n = length(unique(object@data$Counts)))

    # Avoid negative probs
    # kdeProbs <- replace(kdeOrgData$y, kdeOrgData$y < 0, 0)

    # TODO: when nrow(data) > range of values in x (should be impossible... but who knows) sample will presumably fail
    object@randData <- sample(kdeOrgData$x, nrow(object@data), replace = TRUE, prob = kdeOrgData$y)
    # object@randData <- runif(nrow(object@data),
    #                          min = object@minMaxDist$min,
    #                          max = object@minMaxDist$max)
    names(object@randData) <- rownames(object@data)

    # Replace negative values with 0
    object@randData[object@randData < 0] <- 0

    # Adjust randomCounts to depth
    object@randData <- object@depth * object@randData / sum(object@randData)

    # Determine which DE non-flat genes have been assigned a non-expressed status AND have an original
    # value lower than that, then swap their values with non-DE genes
    minDEassignments <- names(object@randData)[object@randData < minDE.value] # Keep IDS from randomCounts separate
    minDEassignments.org <- rownames(object@data)[object@data$Counts < minDE.value]

    badAssignmentDE <- featuresDE.notFlat[featuresDE.notFlat %in% intersect(minDEassignments, minDEassignments.org)]
    goodAssignmentNONDE <- featuresNONDE[! featuresNONDE %in% minDEassignments]

    # TODO: what to do when there are not enough IDS to swap
    if (length(badAssignmentDE) && length(goodAssignmentNONDE)) {

        # Select non-DE IDs to swap. Select the maximum amount available
        swapID <- sample(as.character(goodAssignmentNONDE),
                       size = min(length(badAssignmentDE), length(goodAssignmentNONDE)),
                       replace = FALSE)

        # Take a sample to replace
        swapID.bad <- sample(as.character(badAssignmentDE),
                             size = length(swapID),
                             replace = FALSE)

        swapValues <- object@randData[swapID]

        object@randData[swapID] <- object@randData[swapID.bad]
        object@randData[swapID.bad] <- swapValues

        # Update bad assignments
        badAssignmentDE <- setdiff(badAssignmentDE, swapID.bad)
    }

    # In case there are still bad "DE" cases with too low amount of counts,
    # replace the random ones with a proper number
    if (length(badAssignmentDE)) {
        object@randData[badAssignmentDE] <- rnorm(n = length(badAssignmentDE),
                                                  mean = minDE.value,
                                                  sd = 1)
    }

    # TODO: remove this temp mod
    random.M <- pmax(object@randData, object@data$Counts)
    random.m <- pmin(object@randData, object@data$Counts)

    nonzero.org.counts <- object@data[object@data > 0]

    max.diff <- quantile(nonzero.org.counts, c(0.9))
    min.diff <- quantile(nonzero.org.counts, c(0.3))

    # For flat DEG genes in all groups
    min.diff.flat <- quantile(nonzero.org.counts, c(0.8))
    max.diff.flat <- quantile(nonzero.org.counts, c(0.99))

    diff.M <- random.M - random.m

    index.DE <- names(object@randData) %in% featuresDE.notFlat

    # Maximum difference between M and m
    index.match <- (diff.M > max.diff & index.DE)
    # Minimum difference between M and m
    index.min.match <- (diff.M < min.diff & index.DE)

    diff.mask <- (object@randData > object@data$Counts)


    # Limit max diffs
    object@randData[index.match] <- ifelse(diff.mask[index.match],
                                           # rnorm(length(index.match), mean = max.diff, sd = 0.5),
                                           object@randData[index.match] - (diff.M[index.match] - max.diff),
                                           object@randData[index.match] + (diff.M[index.match] - max.diff))

    # Limit min diffs
    # TODO: when substracting the diff it could be under 0 so it will not assure a "minDE" diff.
    # Try with this to see how it behaves.
    object@randData[index.min.match] <- ifelse(diff.mask[index.min.match],
                                           # rnorm(length(index.match), mean = max.diff, sd = 0.5),
                                           object@randData[index.min.match] + (min.diff - diff.M[index.min.match]),
                                           object@randData[index.min.match] - (min.diff - diff.M[index.min.match]))

    # TODO: remove this?
    # Re-adjust to deptp after changes?
    # object@randData <- object@depth * object@randData / sum(object@randData)


    # filt1 <- sum(object@data[dePos,] == 0)
    # filt2 <- sum(object@randData[dePos] == 0)
    # filt3 <- sum(object@data[dePos,] == 0 & object@randData[dePos] == 0)
    # message(sprintf("Hay %d en originales y %d en random [%d]", filt1, filt2, filt3))

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
        # TODO: remove this, randData now has names.
        # TODO: add min/max.diff.de as in regulators
        flatIndexes <- match(flatProfiles$ID, rownames(object@data))

        # Different treatments between RNA-seq and the rest of simulators
        if ( ! object@regulator) {
            # Select genes
            # dataSubset <- object@data[flatProfiles$ID, ]

            # Iterate over the rows on a unique matrix with the structure:
            # Gene_ID | Profile_Group_1 | ... | Counts_Group_1 | ...
            # TODO: change to data instead of randData?
            # object@data[flatIndexes, ] <- t(apply(cbind(flatProfiles, object@randData[flatIndexes]), 1, function(idRow) {


            # Profiles + Counts
            profileCounts <- cbind(flatProfiles, object@data[flatIndexes, ])

            profileColumns <- dplyr::starts_with("Group", vars = colnames(profileCounts))
            countColumns <- dplyr::starts_with("Counts.Group", vars = colnames(profileCounts))

            object@data[flatIndexes, ] <- t(apply(profileCounts, 1, function(idRow) {
                # Types of profiles (exclude first column)
                idProfiles <- as.character(idRow[profileColumns])
                # Count values (columns after profiles)
                idCounts <- as.numeric(idRow[countColumns])

                # Profiles
                return(ifelse(idProfiles == 'flat',
                       # If the profile is flat, keep the current value
                       idCounts,
                       ifelse(idProfiles == 'enhancer',
                              # If the profile is enhancer, generate a new count value in the range
                              # (actual_value + excludeMargin) to maxValue
                              rep.int(runif(1,
                                            min = min(idCounts[1] + excludeMargin + min.diff.flat, object@max), #object@minMaxDist$max),
                                            max = min(idCounts[1] + excludeMargin + max.diff.flat, object@max)),
                                            # max = object@minMaxDist$max),
                                      simulation@numberGroups),
                              # If the profile is repressor, generate a new count value in the range
                              # minValue to (actual_value - excludeMargin)
                              rep.int(runif(1,
                                            # min = object@minMaxDist$min,
                                            # max = max(idCounts[1] - excludeMargin, object@minMaxDist$min)),
                                            min = max(idCounts[1] - excludeMargin - max.diff.flat, object@min),
                                            max = max(idCounts[1] - excludeMargin - min.diff.flat, object@min)),
                                      simulation@numberGroups)
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
                # object@data[profileSubset$ID, ] <- t(apply(cbind(profileSubset, object@randData[profileSubset$ID]), 1, function(idRow) {
                    # Regulator effect
                    regEffect <- idRow[effectColumn]
                    # Types of profiles (skip first 2 columns: ID and Effect)
                    idProfiles <- as.character(idRow[profileColumns])
                    # Count values (following columns)
                    idCounts <- as.numeric(idRow[countColumns])

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
                                         # rep.int(runif(1, min = min(idCounts[1] + excludeMargin, object@max), max = object@max), simulation@numberGroups),
                                         rep.int(runif(1,
                                                       min = min(idCounts[1] + excludeMargin + min.diff.flat, object@max),
                                                       max = min(idCounts[1] + excludeMargin + max.diff.flat, object@max)),
                                                 simulation@numberGroups),
                                         # If the profile is repressor, generate a new count value in the range
                                         # minValue to (actual_value - excludeMargin)
                                         rep.int(runif(1,
                                                       min = max(idCounts[1] - excludeMargin - max.diff.flat, object@min),
                                                       max = max(idCounts[1] - excludeMargin - min.diff.flat, object@min)),
                                                 simulation@numberGroups)
                                  )
                    ))
                }))
            }
        }

        # Adjust again to depth each column
        object@data <- apply(object@data, 2, function(x) x * object@depth / sum(x))

        # Make sure that no-DE features have the same mean around groups
        # TODO: remove this now?
        if (object@regulator) {

            idDE <- unique(dplyr::filter(simulation@simSettings$geneProfiles[[class(object)]],
                                  ! is.na(Effect))$ID)

            nonDE <- unique(dplyr::filter(simulation@simSettings$geneProfiles[[class(object)]],
                                   ! ID %in% idDE)$ID)


        } else {
            nonDE <- unique(dplyr::filter(simulation@simSettings$geneProfiles$SimRNAseq, DE == FALSE)$ID)
        }

        # Replace with a mean value between the groups
        object@data[nonDE, ] <- apply(object@data[nonDE, ], 1, mean)
    }

    return(object)
})

#' @rdname simulate-methods
#' @aliases simulate,Simulation
setMethod("simulate", signature="Simulator", function(object, simulation) {

    message(sprintf("Starting simulation of %s.", object@name))

    object <- initializeData(object, simulation)

    makeReplicates <- function(counts, groupInfo, timeInfo, profileInfo) {

        message(sprintf("\t- Making replicates for group %d on time %s.", groupInfo, timeInfo))

        counts[is.na(counts) | counts < 0.1] <- 0.1

        varest <- 10^object@replicateParams$a * (counts + 1)^object@replicateParams$b - 1

        # Flat variance estimation lower to avoid false "time profile changes"
        # varest.flat <- 10^object@replicateParams$a * (counts + 1)^(object@replicateParams$b/2) - 1

        # Force a minimum variance
        varest[varest < 0.03] <- 0.03


        # varest[flatProfiles] <- varest.flat[flatProfiles]

        # Apply to every row:
        #   Param x: [mu.noise.cpm, stdev]
        temp <- t(matrix(apply(cbind(counts, varest), 1, function (x) {
            if (x[1] >= x[2]) {
                replis <- rpois(n = simulation@numberReps, lambda = x[1])
            } else {
                replis <- rnbinom(n = simulation@numberReps, size = x[1]^2/(x[2] - x[1]),  mu = x[1])
            }

            return(replis)
        }), nrow = simulation@numberReps))

        flatProfiles <- (profileInfo == 'flat')

        # temp[flatProfiles, ] <- t(apply(cbind(counts[flatProfiles], temp[flatProfiles, ]), 1,
        #                                       function(geneRow) {
        #                                           # Original replicate values
        #                                           mean.original <- geneRow[1]
        #                                           gene.values <- geneRow[-1]
        #
        #                                           diff.mean <- mean(gene.values) - mean.original
        #
        #                                           gene.adjusted <- gene.values - diff.mean
        #
        #                                           return(gene.adjusted)
        #                                       }))


        temp.depth <- apply(temp, 2, function(x) x * object@depth / sum(x))

        # For flat profiles row, re-adjust the values so they are centered around the original mean
        adjusted.means <-  (object@depth * counts) / sum(counts)

        # Generate random means
        noised.means <- rnorm(n = length(adjusted.means),
                             mean = adjusted.means,
                             sd = 0.025 * adjusted.means)

# temp.bak <- temp.depth
        # temp.depth[flatProfiles, ] <- t(apply(cbind(adjusted.means[flatProfiles], temp.depth[flatProfiles, ]), 1,
        temp.depth[flatProfiles, ] <- t(apply(cbind(noised.means[flatProfiles], temp.depth[flatProfiles, ]), 1,
                                              function(geneRow) {
                                                  # Original replicate values
                                                  mean.adjusted <- geneRow[1]
                                                  gene.values <- geneRow[-1]

                                                  # Mean
                                                  diff.mean <- mean(gene.values) - mean.adjusted

                                                  gene.adjusted <- gene.values - diff.mean

                                                  return(gene.adjusted)
                                                  }))

        object@debugInfo <<- simDebug(
            object,
            simulation,
            method = sprintf("makeReplicates[group=%d,time=%d]", groupInfo, timeInfo),
            temp,
            # temp.bak,
            temp.depth,
            adjusted.means
        )

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
            simProfiles <- dplyr::select(simProfiles, ID, dplyr::starts_with("Group"), dplyr::starts_with("Tmax")) %>%
                dplyr::distinct_()
        } else {
            simProfiles <- dplyr::filter(simProfiles, ! is.na(Effect)) %>%
                dplyr::select(ID, dplyr::starts_with("Group"), dplyr::starts_with("Tmax")) %>%
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
    simProfiles <- dplyr::filter(simProfiles, ID %in% rownames(object@data))

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
        # columnsParam <- do.call(cbind, replicate(simulation@numberReps, columnsParam, simplify = FALSE))
        # columnsParam <- do.call(cbind, lapply(columnsParam, function(x) replicate(simulation@numberReps, x, simplify=TRUE)))
        columnsParam <- columnsParam[, rep(colnames(columnsParam), each = simulation@numberReps), drop = FALSE]
        # iterationParam <- rep(iterationParam, times = simulation@numberReps)
        iterationParam <- rep(iterationParam, each = simulation@numberReps)
        idsParam <- simProfiles[, rep('ID', simulation@numberGroups * simulation@numberReps), drop = FALSE]

        tmaxParam <- tmaxParam[, rep(colnames(tmaxParam), each = simulation@numberReps), drop = FALSE]

        # TODO: fix this properly. When working with blocks, sometimes there are duplicated rows
        # because they do not have the tmax columns properly filled.
        simProfiles <- dplyr::group_by(simProfiles, ID) %>%
            dplyr::summarise_all(dplyr::funs(dplyr::first(na.omit(.)))) %>%
            dplyr::ungroup()
    }

    # quant.info <- quantile(object@data[[1]], probs = c(0.25, 0.75))
    # simuRCounts <- runif(nrow(simProfiles), min = quant.info[1], max = quant.info[2])
    # TODO: remove this
    # simuRCounts <- jitter(object@data[[1]], factor = 1)
    # simuRCounts.adj <- sum(object@data[[1]]) * simuRCounts / sum(simuRCounts)

    # Generate time series (if any) with replicates
    object@simData <- data.frame(
        mapply(
            function(counts, profiles, tmaxValues, group, ids) {
                if (! object@pregenerated) {
                    message(sprintf("- Simulating count values for group %d.", group))
                } else {
                    message(sprintf("[Pregenerated replicate] Simulating count values for group %d.", group))
                }


                # Create the coefficients based on the "tmax" parameter
                nt <- length(simulation@times)

                x <- c(0:(nt - 1))

                varT <- x[nt]

                timeVectors <- t(apply(cbind(profiles, tmaxValues), 1, function(rowInfo) {
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

                    # return(as.numeric(unlist(rlang::eval_tidy(simulation@profiles[profileValue],
                                                              # data = data_values))))
                })) %*% rbind(c(rep(1, nt)), x, x * x)


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

                repressionMean <- runif(length(M), min = (M+m)/2, max = M)
                inductionMean <- runif(length(m), min = m, max = (M+m)/2)


                if (! object@pregenerated) {
                    simData <-
                        ifelse(grepl('repression', profiles, fixed = TRUE),
                               repressionMean,
                               inductionMean) +
                        # ((M - m) * timeVectors) +
                        ifelse(grepl('repression', profiles, fixed = TRUE),
                               (repressionMean - m),
                               (M - inductionMean)) * timeVectors +
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
                simData[indexFlat, ] <- ifelse(ids[indexFlat] %in% flatProfiles$ID,
                                               counts[indexFlat],
                                               (counts[indexFlat] + randomCounts[indexFlat])/2) + noiseValues[indexFlat, ]

                # TODO: remove this
                mmnoflat <- (M-m)[! indexFlat]

                debugSimData <- simData
                rownames(debugSimData) <- ids

                object@debugInfo <<- simDebug(
                    object,
                    simulation,
                    method = sprintf("simData[group=%d]", group),
                    debugSimData,
                    mmnoflat
                )

                # Calculate a base SD to keep the same between different makeReplicate calls.
                # That way the flat genes will be more stable during time.
                # groupStdev <- makeReplicates.generateSD(simData[, 1])

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
                            # flatPos = as.data.frame(indexFlat)[, rep(1, length(simulation@times))],
                            # baseStdev = as.data.frame(groupStdev)[, rep(1, length(simulation@times))],
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
            # Coefficients parameters
            tmaxParam,
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

    if (! object@pregenerated) {
        object@simData[object@simData < object@min] <- object@min
        object@simData[object@simData > object@max] <- object@max
    }

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

    # # TODO: REMOVE THIS? TRYING WITH THE SAME randomCounts between groups
    # if (! is.null(rcounts)) {
    #     randomCounts <- rcounts
    # } else {
    #     randomCounts <- runif(length(counts), min = object@min, max = object@max)
    # }

    # Reorder
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
