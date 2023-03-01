#' @include Simulator.R functions.R
#' @import lazyeval
#' @importFrom rlang .data
NULL

setMethod("initialize", signature="MOSimulation", function(.Object, ...) {
    # Set slots
    .Object <- callNextMethod()

    # Rename the simulators names to match the classes names
    names(.Object@simulators) <- ifelse(startsWith(names(.Object@simulators), "Sim"),
                                        names(.Object@simulators),
                                        paste0("Sim", gsub("-", "", names(.Object@simulators)))
                                        )

    # RNA-seq simulator must always be present
    if (! exists('SimRNAseq', where = .Object@simulators)) {
        .Object@simulators[['SimRNAseq']] <- list()
    }

    # Load default data if none is provided
    dataProvided <- sapply(.Object@simulators, is.declared, key = 'data')

    if (! any(dataProvided)) {
        data(sampleData, envir = environment())

        .Object@simulators <- mapply(function(sim, data) {
            # Assign data only on list option, not already instantiated objects or those
            # that already have it on the options.
            if (is.declared(data) && ! inherits(sim, "MOSimulator")) {

                if (! is.declared(sim, "data")) {
                    sim$data <- if (exists('data', where = data)) data$data else NULL
                }

                if (! is.declared(sim, "idToGene")) {
                    sim$idToGene <- if (exists('idToGene', where = data)) data$idToGene else NULL
                }
            }

            return(sim)
        }, sim = .Object@simulators, data = sampleData[names(.Object@simulators)], SIMPLIFY = FALSE)

        if (! is.null(.Object@TFtoGene) && .Object@TFtoGene) {
            .Object@TFtoGene <- sampleData$SimRNAseq$TFtoGene
        }
    }

    # Inherited params from simulation
    inheritParams <- list(
        "noiseFunction",
        "noiseParams",
        "depth",
        "minMaxQuantile",
        "minMaxFC"
        )

    # Initialize simulators
    .Object@simulators <- mapply(function(sim, class) {
        # Directly return the simulator if it has already been initialized
        if (! inherits(sim, "MOSimulator")) {
            # Add default params if not provided in the options list
            for (param in inheritParams) {
                if (! exists(param, sim)) {
                    sim[[param]] <- slot(.Object, param)
                }
            }

            # Initialize the simulator
            sim <- do.call(new, c("Class" = class, sim))
        } else {
            # Set the proper slots in the initialized object
            for (param in inheritParams) {
                slot(sim, param) <- slot(.Object, param)
            }
        }

        # If required, adjust the total number of features to simulate on every omic
        if (nrow(sim@data) > sim@totalFeatures) {

            if (sim@regulator) {
                # TODO: limit the sample to those regulators affecting the possibly
                # limited gene list?
            }

            sim@data <- sim@data[sample(rownames(sim@data), sim@totalFeatures), ,drop=FALSE]
        }

        return(sim)

    }, .Object@simulators, names(.Object@simulators), SIMPLIFY = FALSE)

    # Assign gene identifiers and total number from the data
    .Object@geneNames <- rownames(.Object@simulators[['SimRNAseq']]@data)
    .Object@totalGenes <- nrow(.Object@simulators[['SimRNAseq']]@data)

    # Check gene number
    .Object@diffGenes <- round(.Object@totalGenes * min(.Object@diffGenes, 0.99))

    # Print settings
    simSettings(.Object)

    # Create the simulation settings only if they haven't been passed as an option.
    if (! is.declared(.Object@simSettings)) {

        message("Generating simulation settings for RNA-seq.")

        # Available profiles with probabilities
        profileOptions <- names(.Object@profiles)
        profileProbs <- .Object@profileProbs[profileOptions]

        # If there is no time series, allow only the flat profile
        if (length(.Object@times) < 2) {
            profileProbs <- replace(profileProbs, profileOptions != 'flat', 0)
        }

        # If we consider the TF, force that at least 85% of them to be DE
        if (! is.null(.Object@TFtoGene)) {
            # Rename columns
            colnames(.Object@TFtoGene) <- c("Symbol", "TFgene", "LinkedGene")

            # Exclude those not present on gene expression data
            .Object@TFtoGene <- dplyr::filter(.Object@TFtoGene,
                                              TFgene %in% .Object@geneNames)

            # TF to be DE
            allGeneTF <- unique(.Object@TFtoGene$TFgene)

            sampleDE <- sample(as.character(allGeneTF),
                               size = min(round(length(allGeneTF) * 0.85), .Object@diffGenes),
                               replace = FALSE)

            numDiffGenes <- .Object@diffGenes - length(sampleDE)

            if (numDiffGenes) {
                sampleDE <- c(sampleDE, sample(
                    setdiff(.Object@geneNames, sampleDE),
                    size = numDiffGenes,
                    replace = FALSE
                ))
            }
        } else {
            # Assign random gene names to DE
            sampleDE <- sample(
                .Object@geneNames,
                size = .Object@diffGenes,
                replace = FALSE
            )
        }

        # Assign gene names to NonDE
        sampleNonDE <- setdiff(.Object@geneNames, sampleDE)

        # Randomly assign a profile to every DE gene
        # If we have only ONE group, remove flat profile from the
        # DE options because otherwise could not be considered DE genes.
        profileProbsDE <- profileProbs

        if (.Object@numberGroups < 2) {
            profileProbsDE <- replace(profileProbs, profileOptions == 'flat', 0)
        }

        # If number of times is smaller than 3, remove transitory profiles
        if (length(.Object@times) < 3) {
            profileProbsDE <- replace(profileProbsDE, grep("transitory", profileOptions), 0)
        }

        # Construct the data frame with DE genes profiles
        profilesDE <-
            if (.Object@diffGenes)
                data.frame(
                    ID = sampleDE,
                    DE = TRUE,
                    replicate(
                        .Object@numberGroups,
                        sample(
                            profileOptions,
                            .Object@diffGenes,
                            replace = TRUE,
                            prob = profileProbsDE
                        )
                    ), stringsAsFactors = FALSE)
            else
                stop("There are no DEG profiles available, please increase the amount of DE genes.")


        # Assign a random "Tmax" for transitory profiles

        # Declare helper variable
        tmax.T <- length(.Object@times) - 1

        for (g in seq(.Object@numberGroups)) {
            profilesDE[, paste0("Tmax.Group", g)] <- ifelse(grepl("transitory", profilesDE[, 2 + g]), # Third column = Group 1
                                                            stats::runif(nrow(profilesDE), min = 0.25 * tmax.T, max = 0.75 * tmax.T),
                                                            NA)
        }

        # Randomly assign a profile to every NonDE
        profilesNonDE <-
            if (length(sampleNonDE))
                data.frame(
                    ID = sampleNonDE,
                    DE = FALSE,
                    replicate(.Object@numberGroups,
                              rep('flat', times = length(sampleNonDE))),
                    replicate(.Object@numberGroups,
                              rep(NA, times = length(sampleNonDE))),
                    stringsAsFactors = FALSE)
            else
                stop("There are no non-DEG profiles available, please reduce the amount of DE genes.")

        # Set correct colnames
        # Important: pattern "Group[n]" is used to subset with dplyr later
        colNames <- paste0("Group", seq_len(.Object@numberGroups))

        # Symbols needed to split the DF by dplyr
        dplyrGroup <- lapply(colNames, as.symbol)

        colnames(profilesDE) <- colnames(profilesNonDE) <- c("ID", "DE",
                                         colNames, utils::tail(colnames(profilesDE),
                                                     .Object@numberGroups))

        # Merge DE + nonDE for easier access
        # profilesExpr <- rbind(profilesDE, profilesNonDE)
        profilesAll <- rbind(profilesDE, profilesNonDE)

        # Create a profile table for each regulator
        # Note: prevent harmless warning about converting factors to character vector
        suppressWarnings(profilesReg <- sapply(.Object@simulators[names(.Object@simulators) != 'SimRNAseq'], function(sim) {
            message(sprintf("Generating simulation settings for %s.", sim@name))

            # First assign an effect to each regulator, including NE
            allRegulators <- rownames(sim@data)

            # Remove NE
            availableEffects <- sim@regulatorEffect[grep("NE", names(sim@regulatorEffect), invert = TRUE)]

            # Number of regulators with effect
            numberWithEffect <- length(allRegulators) * sum(as.numeric(availableEffects))

            # Select genes ID regulator ID from all genes considered
            regTable <- IDtoGenes(sim, rownames(sim@data), simplify = FALSE)

            # Choose from regulators associated to some gene
            availableRegulators <- unique(regTable$ID)

            # If there are not enough, chose the available ones
            chosenRegulators <- sample(x = availableRegulators,
                                       size = min(length(availableRegulators), numberWithEffect),
                                       replace = FALSE)
            effectByID <- data.frame(
                ID = allRegulators,
                Effect = ifelse(allRegulators %in% chosenRegulators,
                                TRUE, # Do not assign the effect as it will be assigned later
                                "NE"
                                )
                )

            # Skip the process
            if (! nrow(regTable))
                stop(sprintf("No genes were retrieved from the association list for %s, be sure to provide the correct table.", sim@name))

            # Randomly assign an effect to every regulator from the available
            # options defined on every simulator, excluding NE effect.
            regTable$Effect <- sample(
                names(availableEffects),
                size = nrow(regTable),
                replace = TRUE,
                prob = as.numeric(availableEffects)
            )

            # Set non-DEG or non-selected regulators to <NA>
            discardedReg <- dplyr::filter(effectByID, Effect == "NE") %>% dplyr::pull(ID)

            regTable <- dplyr::mutate(regTable, Effect = replace(Effect, (! Gene %in% sampleDE) | (ID %in% discardedReg), NA))

            # Copy the effect to each group
            regTable[, paste0("Effect.Group", seq(.Object@numberGroups))] <- regTable[, rep("Effect", times =.Object@numberGroups)]

            # Initialize tmax.Group columns
            # Copy the effect to each group
            regTable[, paste0("Tmax.Group", seq(.Object@numberGroups))] <- NA

            # Allow simulator class to adjust the generated profiles
            # (i.e. methyation using blocks)
            regTable <- adjustProfiles(sim, .Object, regTable, step = "Effect")

            # Select all (hence the double check) rows duplicated, only for DE genes.
            # Duplicated ID means that the regulator is linked to more than one gene.
            regDupsDE <-  (duplicated(regTable$ID) |
                            duplicated(regTable$ID, fromLast=TRUE)) &
                            (regTable$Gene %in% sampleDE) &
                            (! regTable$ID %in% discardedReg)

            # For every regulator with more than one regulated gene, decide which
            # class (combination of profiles among conditions) to regulate.
            message("- Linking to gene classes...")

            # Legacy code to enable Tmax
            # df: regulator df
            # genesTmax: matrix with gene Tmax values in the proper order
            # assignTmax <- function(df, genesTmax) {
            #     for (group in seq(.Object@numberGroups)) {
            #         colName <- paste0("Tmax.Group", group)
            #         df <- dplyr::mutate(df, !!colName := (if(is.na(genesTmax[[group]]) || is.infinite(genesTmax[[group]])) NA else stats::runif(n = 1, 0.25 * tmax.T, as.numeric(genesTmax[group]))))
            #     }
            #
            #     return(df)
            # }

            # Populate tmax values with the ones found in gene settings table
            regTable <- dplyr::select(regTable, -dplyr::starts_with("Tmax")) %>%
                dplyr::left_join(dplyr::select(profilesDE, ID, dplyr::starts_with("Tmax")), by = c("Gene" = "ID"))

            if (any(regDupsDE)) {
                # Reset effect on the duplicated DE
                # regTable[regDupsDE, c("Effect")] <- NA

                # Choose a new effect
                # Function to use inside classifyDups
                # .selectEffect <- function(groupRow) {
                #     # Relies on the row having these columns (check function classifyDups)
                #     percEnhancer <- max(0, groupRow$Percentage.Enhancer, na.rm = T)
                #     percRepressor <- max(0, groupRow$Percentage.Repressor, na.rm = T)
                #
                #     # In case of tie, select a new effect excluding NE option
                #     if (percEnhancer == percRepressor) {
                #         # If both % are zero, this means the regulator had only ONE class
                #         # that was chosen, and that class was assigned initially as NE,
                #         # thus we keep that.
                #         if (sum(percEnhancer, percRepressor) == 0) {
                #             selEffect <- "NE"
                #         } else {
                #             sampleOptions <- grep("NE", names(sim@regulatorEffect), invert = T)
                #
                #             selEffect <- sample(names(sim@regulatorEffect)[sampleOptions],
                #                                 size = 1,
                #                                 prob = as.numeric(sim@regulatorEffect[sampleOptions]))
                #         }
                #     } else {
                #         selEffect <- if (percEnhancer > percRepressor) 'enhancer' else 'repressor'
                #     }
                #
                #     return(selEffect)
                # }

                classifyDups <- function(profileGroup) {
                    # Select classes of the genes associated with the regulator
                    # Note: dplyr renames the unneeded ID column of "." to ID.y
                    # TODO: room for optimization, remove this join from the loop.
                    regGenes <- dplyr::inner_join(profilesDE, profileGroup, by = c("ID" = "Gene")) %>%
                        dplyr::select(-DE, -ID.y, -dplyr::starts_with("Effect.Group"), -dplyr::ends_with(".y"), dplyr::starts_with("Keep."))

                    # Initialize Tmax.GroupX considering all, for those cases in which
                    # there is only one row.
                    tmaxGroup <- dplyr::select(regGenes, dplyr::starts_with("Tmax.")) %>%
                        dplyr::summarise_all(min, na.rm = TRUE)

                    # Skip the process when there is only one row (regulators affects more than 1 gene
                    # but the others have been classified has non-DE) or if there is no effect at all
                    # in the genes.
                    if (nrow(profileGroup) > 1 && ! all(is.na(regGenes$Effect) | regGenes$Effect == "NE")) {
                        # Split by class and count number of genes on each one.
                        #
                        # I.e. 3 conditions:
                        # induction - induction - repression: X rows
                        # induction - induction - induction: Y rows
                        classSplit <-
                            dplyr::group_by_(regGenes, .dots = dplyrGroup)

                        # Select genes of the majoritary class.
                        #
                        # The majoritary class is selected based on the % of the
                        # most assigned Effect in the group, then, in case of tie, the number
                        # of affected genes is taken into account. If there is another tie,
                        # one is chosen at random or two classes if there are totally
                        # opposite classes.
                        #
                        # Force sampling options to character to prevent a seq when there
                        # is only a result.
                        groupSizes <- dplyr::group_size(classSplit)

                        # If there is only 1 groupSize, skip the process.
                        # If the max % effect is 0, it means that there were multiple groups but all of them
                        # had NE assigned by default.
                        if (length(groupSizes) > 1){# || max(effectSizes$MaxPercentage) == 0) {
                            # Remove groups that have no effect at all
                            # I.e. in methylation we keep <NA> effect up to this point
                            # so there could be two majoritary groups and one of them
                            # without effect.
                            excludedGroups <- dplyr::summarize(classSplit, Exclude = all(is.na(Effect))) %>%
                                dplyr::pull(Exclude)

                            # Remove those groups that have no effect (in methylation)
                            groupSample <- setdiff(which(groupSizes == max(groupSizes[! excludedGroups])), which(excludedGroups))

                            groupSelected <- sample(as.character(groupSample), size = 1)

                            selGroup <- dplyr::slice(regGenes,
                                                     which(dplyr::group_indices(classSplit) == groupSelected))

                            # Retrieve the most used effect
                            propEffect <- prop.table(table(selGroup$Effect))
                            selEffect <- names(propEffect)[propEffect == max(propEffect)]

                            # In case of tie, select at random
                            if (length(selEffect) > 1) {
                                selEffect <- sample(
                                    names(availableEffects),
                                    size = 1,
                                    replace = TRUE,
                                    prob = as.numeric(availableEffects)
                                )
                            }

                            # Select the minimum "Tmax" value for the genes in the group, and generate one
                            # for the regulator.
                            #
                            # Note: there should not be any NA values (in a group) when there is a transitory effect (na.rm not needed)
                            tmaxGroup <- dplyr::select(selGroup, dplyr::starts_with("Tmax.")) %>%
                                dplyr::summarise_all(min, na.rm = TRUE)

                            selGenes <- selGroup$ID

                            # Save profile info (remove ID column)
                            groupCols <- regGenes[, - 1]

                            # Append new columns with associated effect
                            # Add column of effect in Group <i>
                            matchTable <- match(regGenes$ID, profileGroup$Gene)

                            # Selected profile
                            selectedProfile <- dplyr::select(selGroup, dplyr::starts_with("Group")) %>%
                                dplyr::distinct() %>% unlist()

                            # Opposite effect
                            revEffect <- if (selEffect == "NE") NA else setdiff(c('enhancer', 'repressor'), selEffect)

                            for (i in seq(.Object@numberGroups)) {
                                # Opposite class for this group
                                revClass <- stri_replace_all_fixed(
                                    selGroup[1, i + 1], # Skip ID column
                                    c(".induction", ".repression", "*"),
                                    c("*repression", "*induction", "."),
                                    vectorize_all = FALSE
                                )

				## CM_Fix
                                ## When simulating methylation for more than 1 
                                # timepoint, there are repeated matches between 
                                # profile group and groupCols, so there are 
                                # repeated indexes in matchTable. When this 
                                # happens, get the unique overlap, as the 
                                # methylation block identifyer keeps the info
                                # of which CpGs are involved

                                if (sim@name == 'Methyl-seq'){
                                    if (any(duplicated(matchTable))){
                                        matchTable <- unique(matchTable)
                                        profileGroup <- profileGroup[c(matchTable), ]
                                        rownames(profileGroup) <- c(1:length(matchTable))
                                        groupCols <- groupCols[matchTable, ]
                                        rownames(groupCols) <- c(1:length(matchTable))
                                        matchTable <- c(1:length(matchTable))
                                    }
                                }

                                profileGroup[matchTable, paste0("Effect.Group", i)] <- ifelse(
                                    # Check if group is the same as the selected class in this group
                                    groupCols[, i] == as.character(selectedProfile[i]), #selGroup[1, i + 1]),
                                    # If that is the case, choose selected effect
                                    selEffect,
                                    # If not, check if it has the opposite profile
                                    ifelse(
                                        groupCols[, i] == revClass,
                                        # If it has the opposite profile, assign the contrary effect
                                        revEffect,
                                        # If not, assign an NA
                                        NA
                                    )
                                )
                            }

                            profileGroup <- dplyr::mutate(profileGroup, Effect = ifelse(Gene %in% selGenes, selEffect, NA))

                        } else {
                            # Select the group effect (either only one group or multiple but all of them with NE effect)
                              selEffect <- sample(
                                names(availableEffects),
                                size = 1,
                                replace = TRUE,
                                prob = as.numeric(availableEffects)
                            )

                            # Return the original data frame with the same effect for all the genes in the group
                            # Note: use paste as function because newer versions of dplyr cannot work with "~ selEffect"
                            # anymore.
                            profileGroup <- dplyr::mutate(profileGroup, Effect = selEffect) %>%
                                dplyr::mutate_at(dplyr::vars(dplyr::starts_with("Effect.Group")), ~ paste0(selEffect))
                        }
                    }

                    # Assign tmax
                    for (group in seq(.Object@numberGroups)) {
                        colName <- paste0("Tmax.Group", group)

                        profileGroup <- dplyr::mutate(profileGroup,
                                                      !!colName := (if(is.na(tmaxGroup[[group]]) || is.infinite(tmaxGroup[[group]])) NA
                                                                    else as.numeric(tmaxGroup[group])))
                    }

                    # Replace "NE" effect with NA
                    # TODO: check if this is really required, no NE at this point?
                    profileGroup <- dplyr::mutate_at(profileGroup,
                                                     dplyr::vars(dplyr::contains("Effect")),
                                                     list(~ifelse(. == "NE", NA, .)))

                    # Return the processed (or not) profileGroup
                    return(profileGroup)
                }

                regTable[regDupsDE, ] <- dplyr::group_by(regTable[regDupsDE, ], ID) %>% dplyr::do(classifyDups(.))
            }

            message("- Adjusting regulator effect")

            # Expression profiles of all expressed genes
            condTable <- dplyr::left_join(regTable, dplyr::select(profilesAll, -dplyr::starts_with("Tmax.")),
                                          by = c("Gene" = "ID")) %>% dplyr::ungroup()

            subCondTable <- dplyr::select(condTable, ID, Gene, Effect,
                                          dplyr::starts_with("Effect.Group"),
                                          dplyr::starts_with("Tmax.Group"),
                                          dplyr::starts_with("Keep."))

            # Add those regulators not associated with any gene
            allIds <- rownames(sim@data)

            if (! sim@pregenerated) {
                notInTable <- allIds[!allIds %in% condTable$ID]

                subCondTable <- dplyr::bind_rows(subCondTable, data.frame(ID = notInTable))
                condTable <- dplyr::bind_rows(condTable, data.frame(ID = notInTable))
            } else {
                # We assume that Keep.* is the one including the original ID like methyl-seq simulator
                colID <- colnames(regTable)[grep("Keep.", colnames(regTable))]

                notInTable <- dplyr::select(dplyr::ungroup(regTable), "ID", colID)

                notInTable <- notInTable[! as.character(notInTable[, colID]) %in% as.character(condTable[, colID]), ]

                subCondTable <- dplyr::bind_rows(subCondTable, notInTable)
                condTable <- dplyr::bind_rows(condTable, notInTable)
            }

            # Modify the profile of every condition based on the regulator effect
            # on the gene. If the regulator has an "enhancer" effect, the profile
            # will be the same as in RNA-seq, whereas with "repressor" effect, the
            # profile will be the opposite. For no effect, it will be "flat".
            #
            # Examples:
            #   Enhancer:
            #       RNA-seq: continuous.induction
            #       Regulator: continuous.induction
            #   Repressor:
            #       RNA-seq: continous.induction
            #       Regulator: continuous.repressor
            #
            # Note: the fixed pattern/replacement of stri_replace_* overwrites previous
            # results, hence the need to substitute the dot with asterisk first.
            settingsTable <- data.frame(
                    # Variable containing required columns and also those marked as "keep.Name"
                    # in some particular simulators (i.e. methylation)
                    subCondTable,
                    # dplyr::select(condTable, ID, Gene, Effect, dplyr::starts_with("Keep.")),
                    apply(dplyr::select(condTable, dplyr::starts_with("Group")), 2, function(col)
                        ifelse(
                            is.na(condTable$Effect),
                            # Flat profile for no effect regulator
                            'flat',
                            ifelse(
                                condTable$Effect == "enhancer",
                                # Same profile as RNA-seq for enhancer regulators
                                col,
                                # Opposite profile for repressor regulators
                                stri_replace_all_fixed(
                                    col,
                                    c(".induction", ".repression", "*"),
                                    c("*repression", "*induction", "."),
                                    vectorize_all = FALSE
                                )
                            )
                        ))
                , stringsAsFactors = FALSE)

            # Let individual simulator final adjustments after the processing step
            settingsTable <- adjustProfiles(sim, .Object, settingsTable, step = "Groups")

            return(settingsTable)
        }, simplify = FALSE))

        featuresAll <- sapply(profilesReg, function(simulatorProfiles) {

            test <- dplyr::select(simulatorProfiles, ID, dplyr::starts_with(("Effect.Group")))

            # nonDE.pos <- dplyr::select(simulatorProfiles, dplyr::starts_with(("Effect.Group"))) %>%
            #     pmap(~ all(is.na(.))) %>%
            #     unlist()
            #
            # nonDE.ids <- unique(simulatorProfiles$ID[nonDE.pos])

            nonDE.ids <- unique(dplyr::group_by(simulatorProfiles, ID) %>%
                                    dplyr::filter(all(is.na(Effect))) %>%
                                    dplyr::pull(ID))

            DE.ids <- setdiff(simulatorProfiles$ID, nonDE.ids)

            output <- list(
                "DE" = DE.ids,
                "nonDE" = nonDE.ids,
                "noiseNonDE" = sample(nonDE.ids, size = ceiling(length(nonDE.ids) * 0.05))
            )

            return(output)
        }, simplify = FALSE)

        featuresAll$SimRNAseq <- list(
            "DE" = sampleDE,
            "nonDE" = sampleNonDE,
            "noiseNonDE" = sample(sampleNonDE, size = ceiling(length(sampleNonDE) * 0.05))
        )

        # Expressed genes profiles (DE + NonDE)
        profilesReg$SimRNAseq <- profilesAll#profilesExpr

        # In case there are no times in exp design, change the initial count values between
        # samples on DE genes. The associated regulators must change as well,
        # depending on the linked effect.
        if (.Object@numberGroups > 1) {
            # Determine the rows with the same condition on all columns
            # Change only base count values when the profiles are flat on all conditions,
            # otherwise they are already DE genes.
            # sameCond <- apply(apply(dplyr::select(profilesDE, dplyr::starts_with("Group")), 2, '==', 'flat'), 1, all)
            sameCond <- dplyr::select(profilesDE, dplyr::starts_with("Group")) %>%
                purrr::map(`==`, 'flat') %>%
                purrr::pmap_lgl(all)

            # NULL by default
            FlatRNAseq <- NULL

            if ((nSameCount <- sum(sameCond))) {
                message(sprintf("Creating settings to change count values on %d DE genes with the same flat profile on all groups.", nSameCount))

                # Select ID, Group1-GroupN (remove DE -2nd- column)
                sameCondDE <- t(apply(profilesDE[sameCond, -c(2), drop = FALSE], 1, function(geneRow) {
                    # TODO: consider DE as one gene that changes in one condition compared to any of the others,
                    # or just consider the fist group as a reference?
                    #
                    # At this time exclude group 1 from being modified.

                    # Convert to character vector to avoid unexpected sample behaviour
                    # in case of only two groups
                    colChange <- sample(as.character(2:.Object@numberGroups), size = sample.int(.Object@numberGroups - 1, 1))

                    # Associate a trend to every group to change. If enhancer is set,
                    # the new generated count value will be higher than the original,
                    # and the opposite with repressor effect. This has nothing to do with
                    # regulator effects and are used only as indicators.
                    #
                    # Note: add 1 to the indexes because the first column is the gene ID
                    geneRow[as.numeric(colChange) + 1] <- sample(c('enhancer', 'repressor'), size = 1)

                    return(geneRow)
                }))

                FlatRNAseq <- data.frame(sameCondDE, stringsAsFactors = FALSE)
            }

            # Generate "FlatGroup" list for each regulator
            profilesReg$FlatGroups <- lapply(.Object@simulators[names(.Object@simulators) != 'SimRNAseq'], function(regulator) {
                # Skip if there aren't any flat DEG genes in gene expression
                if (is.null(FlatRNAseq))
                    return(NULL)

                profileSubsetID <- dplyr::rename(FlatRNAseq, Gene = ID) %>%
                    dplyr::inner_join(profilesReg[[class(regulator)]][, c('ID', 'Effect', 'Gene')], by = c("Gene" = "Gene")) %>%
                    dplyr::filter(! is.na(Effect)) %>% dplyr::select(ID)

                return(profileSubsetID)
            })

            profilesReg$FlatGroups$SimRNAseq <- FlatRNAseq
        }

        # Resolve references of profiles to proper values
        message("Finishing generation of configuration settings.")

        .Object@simSettings <- list(
            geneProfiles = profilesReg,

            geneSamples = list(
                DE = sampleDE,
                nonDE = sampleNonDE
            ),

            featureSamples = featuresAll,

            expDesign = list(
                times = .Object@times,
                numberReps = .Object@numberReps,
                numberGroups = .Object@numberGroups
            )
        )

        message("Configuration generated.\n")
    }

    return(.Object)
})

setMethod("simulate", signature="MOSimulation", function(object) {
    object@simulators <- sapply(object@simulators, simulate, object)

    # TODO: as a special-case omic, keep it here or move to a proper location?
    # Return TF <-> gene data
    if (! is.null(object@TFtoGene) && nrow(object@TFtoGene)) {

        tableTF <- object@TFtoGene

        # TF data (symbol -> geneTF)
        # Select only genes present in data
        expression.data <- object@simulators$SimRNAseq@simData
        expression.settings <- object@simSettings$geneProfiles$SimRNAseq
        expression.settings.flat <- object@simSettings$geneProfiles$FlatGroups$SimRNAseq

        tableTF.data <- unique(tableTF[, c("Symbol", "TFgene")])

        # Select the proper rows and change the rownames
        simDataTF <- expression.data[tableTF.data$TFgene, ]
        rownames(simDataTF) <- tableTF.data$Symbol

        # Generate the association

        # Symbol | LinkedGene | Effect | Effect.Group1 | ... | Effect.GroupX | Group1 | ... | GroupX
        message("Generating settings for TF based on gene expression settings...")

        settingsTF <- dplyr::inner_join(tableTF, expression.settings, by = c("TFgene" = "ID")) %>%
            dplyr::inner_join(expression.settings, by = c("LinkedGene" = "ID"), suffix = c(".TF", ".Gene")) %>%
            dplyr::left_join(expression.settings.flat, by = c("TFgene" = "ID")) %>%
            dplyr::left_join(expression.settings.flat, by = c("LinkedGene" = "ID"), suffix = c(".FlatTF", ".FlatGene")) %>%
            dplyr::select(ID = Symbol, Gene = LinkedGene, DE.TF, DE.Gene, dplyr::starts_with("Group")) %>%
            dplyr::mutate(Effect = NA)

        # Check effect for every group
        for (i in seq(object@numberGroups)) {

            # TF profile in this group
            TFprofile <- as.character(settingsTF[, sprintf("Group%d.TF", i)])

            # Linked genes profiles in this group
            GeneProfile <- as.character(settingsTF[, sprintf("Group%d.Gene", i)])

            # FLAT TF profile in this group
            flatTFprofile <- as.character(settingsTF[, sprintf("Group%d.FlatTF", i)])

            # FLAT Linked genes profiles in this group
            flatGeneProfile <- as.character(settingsTF[, sprintf("Group%d.FlatGene", i)])

            # Opposite effect of linked genes
            GeneProfile.rev <- stri_replace_all_fixed(
                GeneProfile,
                c(".induction", ".repression", "*"),
                c("*repression", "*induction", "."),
                vectorize_all = FALSE
            )

            # Assign the effect in the group
            settingsTF[, sprintf("Effect.Group%d", i)] <- ifelse(
                # For a proper effect, both genes (regulator and regulated)
                # must be DE
                ! (settingsTF$DE.TF & settingsTF$DE.Gene),
                NA,
                ifelse(
                    # If none of the genes are DE flat, check normally
                    (is.na(flatTFprofile) & is.na(flatGeneProfile)),
                    ifelse(
                        TFprofile == GeneProfile,
                        "enhancer",
                        ifelse(
                            TFprofile == GeneProfile.rev,
                            "repressor",
                            NA
                        )
                    ),
                    # If both genes are flat, check if the effect in the
                    # group is the same or not. Otherwise if one is flat
                    # and the other not, assign NA
                    ifelse(
                        # Both genes DE flat and the effect in the group is different than flat (that is,
                        # it does not change the original counts in that group, that is kept as reference)
                        ((! is.na(flatTFprofile)) & (! is.na(flatGeneProfile))) & flatTFprofile != "flat",
                        ifelse(
                            flatTFprofile == flatGeneProfile,
                            "enhancer",
                            "repressor"
                        ),
                        NA
                    )
                )
            )

            # For flat-flat profiles we also need to check the effect
            # in each group

            # Override the global "Effect" column if it has any effect
            settingsTF[, "Effect"] <- ifelse(
                is.na(settingsTF[, "Effect"]),
                      settingsTF[, sprintf("Effect.Group%d", i)],
                      settingsTF[, "Effect"]
            )
        }

        # Assign an effect
        settingsTF <- dplyr::select(settingsTF,
                                    ID, Gene, Effect,
                                    dplyr::starts_with("Effect"),
                                    dplyr::ends_with(".TF"))

        # Rename columns
        colnames(settingsTF) <- gsub(".TF", "", colnames(settingsTF))

        # Append empty tmax columns just for compatibility with the others
        settingsTF[, paste0("Tmax.Group", seq(object@numberGroups))] <- NA

        # Assign to the simulation object
        object@simSettings$geneProfiles$SimTF <- settingsTF

        # Clone the RNA-seq simulator then override some params
        TFobject <- object@simulators$SimRNAseq

        TFobject@idToGene <- tableTF[, c("Symbol", "LinkedGene")]
        colnames(TFobject@idToGene) <- c("ID", "Gene")

        TFobject@simData <- simDataTF
        TFobject@name <- "TF"
        TFobject@regulator <- TRUE

        class(TFobject) <- "SimTF"

        object@simulators$SimTF <- TFobject
    }

    return(object)
})


setMethod("show", signature="MOSimulation", function(object) {
    # TO DO: temp...
    simSettings(object)
})

setMethod("simSettings", signature="MOSimulation", function(object) {
    cat(sprintf("Simulation settings of class %s:\n", class(object)))
    cat(sprintf("- Default depth: %d\n", object@depth))
    cat(sprintf("- Total genes: %d\n", object@totalGenes))
    cat(sprintf("- Dif. expressed genes: %d\n", object@diffGenes))
    cat(sprintf("- Replicates: %d\n", object@numberReps))
    cat(sprintf("- Factor levels (groups): %d\n", object@numberGroups))
    cat(sprintf("- Time vector length: %d\n", length(object@times)))

    if (is.declared(object@simSettings)) {
        # Print anything else
    }
})

setValidity("MOSimulation", function(object) {
    # Initialize list of errors
    errors <- c()

    # Make sure the design allows DE genes
    if (length(object@times) < 2 && object@numberGroups < 2)
        errors <- c(errors, "The design must have a minimum of 2 times or 2 groups.")

    # Check that noiseParams has names
    if (is.declared(object@noiseParams) & is.null(names(object@noiseParams)))
        errors <- c(errors, "The 'noiseParams' parameter must be a named list.")

    # Convert the names of the simulators to the real classnames
    names(object@simulators) <- ifelse(startsWith(names(object@simulators), "Sim"),
                                       names(object@simulators),
                                       paste0("Sim", gsub("-", "", names(object@simulators)))
                                       )

    # Exclude those simulators already initialized
    simuInitialized <- sapply(object@simulators, inherits, what = "MOSimulator")

    # Simulators excluded from association checking (add RNA-seq)
    assocExcluded <- c(names(simuInitialized), "SimRNAseq")

    # If we provide real data it needs to be done for every loaded
    # simulator and also with the association list.
    dataProvided <- sapply(object@simulators, is.declared, key = "data")

    if (length(object@simulators) > 1) {
        # Check the provided association list. RNAseq is the only one with "NA" by default because
        # it doesn't need one.
        assocProvided <-
            sapply(object@simulators[! names(object@simulators) %in% assocExcluded], is.declared, key = 'idToGene')

        # Currently methylation does not accept custom data
        # TODO: rewrite this check to make it more general instead of relying on the simulator name.
        if ('SimMethylseq' %in% names(dataProvided) && dataProvided['SimMethylseq'] == TRUE)
            message("Currently there is no support for including custom methylation data. You need to provide the 'idToGene' slot, which will be used to retrieve the number and location of the CpGs.")

        # Remove Methylation simulator from the list and avoid checking the already
        # initialized simulators.
        dataProvided <- dataProvided[names(dataProvided) != 'SimMethylseq' &
                                         ! names(dataProvided) %in% names(simuInitialized)]

        # If at least one simulator has provided data, all of them must have
        # both data and association list.
        # Consider also the methylation simulation which can create its own data
        # using the association list.
        if (any(dataProvided, assocProvided) && ! all(assocProvided, dataProvided))
            errors <- c(errors, "If you provide an initial sample for one simulator, you need to do it for all the other ones, and include the association list as well.")
    }

    # Load default data if none is provided
    if (any(dataProvided)) {
        # Check that initial sample has only 1 column.
        # Skip initialized simulator.
        dataCheck <- sapply(object@simulators[! names(object@simulators) %in% names(simuInitialized)], function(sim) {
            # If provided, data should have only 1 column
            return(if (is.declared(sim, 'data')) ncol(sim$data) == 1 else TRUE)
        })

        # Check that every association list contains required columns
        assocCheck <-
            sapply(object@simulators[! names(object@simulators) %in% assocExcluded], function(sim, genes) {
                return(all(colnames(sim$idToGene) %in% c('ID', 'Gene')))
            }, genes = rownames(object@simulators[['SimRNAseq']]$data))

        if (! all(dataCheck, assocCheck))
            errors <- c(errors, "The initial sample data must have only 1 column with counts, and all association lists must contain genes present on RNA-seq data.")
    }

    return(if(length(errors)) errors else TRUE)
})
