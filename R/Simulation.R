#' @include Simulator.R functions.R
NULL

#' @rdname initialize-methods
#' @aliases initialize,Simulation
setMethod("initialize", signature="Simulation", function(.Object, ...) {
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
        .Object@simulators <- mapply(function(sim, data) {
            # Assign data only on list option, not already instantiated objects or those
            # that already have it on the options.
            if (is.declared(data) && ! inherits(sim, "Simulator")) {

                if (! is.declared(sim, "data")) {
                    sim$data <- if (exists('data', where = data)) data$data else NULL
                }

                if (! is.declared(sim, "idToGene")) {
                    sim$idToGene <- if (exists('idToGene', where = data)) data$idToGene else NULL
                }
            }

            return(sim)
        }, sim = .Object@simulators, data = sampleData[names(.Object@simulators)], SIMPLIFY = FALSE)
    }

    # Set random number generator seed
    set.seed(.Object@randomSeed)

    # Inherited params from simulation
    inheritParams <- list(
        "noiseFunction",
        "noiseParams",
        "depth",
        "minQuantile"
        )

    # Initialize simulators
    .Object@simulators <- mapply(function(sim, class) {
        # Directly return the simulator if it has already been initialized
        if (! inherits(sim, "Simulator")) {
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
        # sim@totalFeatures <- min(nrow(sim@data), sim@totalFeatures)
        if (nrow(sim@data) > sim@totalFeatures) {

            if (sim@regulator) {
                # TODO: limit the sample to those regulators affecting the possibly
                # limited gene list?
            }

            # TODO: keep copy of the original data?
            sim@data <- dplyr::sample_n(sim@data, sim@totalFeatures)
        }

        return(sim)

    }, .Object@simulators, names(.Object@simulators), SIMPLIFY = FALSE)

    # Assign gene identifiers and total number from the data
    .Object@geneNames <- rownames(.Object@simulators[['SimRNAseq']]@data)
    .Object@totalGenes <- nrow(.Object@simulators[['SimRNAseq']]@data)

    # Check gene number
    # .Object@exprGenes <- sum(.Object@simulators[['SimRNAseq']]@data$Counts > 0)
    # .Object@exprGenes <- min(.Object@totalGenes, .Object@exprGenes)
    # .Object@diffGenes <- round(.Object@exprGenes * min(.Object@diffGenes, 1))
    .Object@diffGenes <- round(.Object@totalGenes * min(.Object@diffGenes, 1))

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

        # Assign random gene names to DE
        sampleDE <- sample(
            .Object@geneNames,
            size = .Object@diffGenes,
            replace = FALSE
            )

        # Assign random gene names to NonDE
        # sampleNonDE <- sample(
        #     setdiff(.Object@geneNames, sampleDE),
        #     size = .Object@exprGenes - .Object@diffGenes,
        #     replace = FALSE
        #     )

        # Expressed genes (DE + NonDE) identifiers
        # sampleExpr <- union(sampleNonDE, sampleDE)

        # Non-expressed genes identifiers
        # sampleNonExpr <- setdiff(.Object@geneNames, sampleExpr)
        sampleNonDE <- setdiff(.Object@geneNames, sampleDE)

        # Randomly assign a profile to every DE gene
        # If we have only ONE group, remove flat profile from the
        # DE options because otherwise could not be considered DE genes.
        profileProbsDE <- profileProbs

        if (.Object@numberGroups < 2) {
            profileProbsDE <- replace(profileProbs, profileOptions == 'flat', 0)
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
                data.frame()

        # Randomly assign a profile to every NonDE
        # profilesNonDE <-
        #     if (.Object@exprGenes > .Object@diffGenes)
        #         data.frame(
        #             ID = sampleNonDE,
        #             DE = FALSE,
        #             replicate(.Object@numberGroups,
        #                       rep('flat', times = .Object@exprGenes - .Object@diffGenes)),
        #             stringsAsFactors = FALSE)
        #     else
        #         data.frame()

        profilesNonDE <-
            if (.Object@totalGenes > .Object@diffGenes)
                data.frame(
                    ID = sampleNonDE,
                    DE = FALSE,
                    replicate(.Object@numberGroups,
                              rep('flat', times = .Object@totalGenes - .Object@diffGenes)),
                    stringsAsFactors = FALSE)
            else
                data.frame()

        # Set correct colnames
        # Important: pattern "Group[n]" is used to subset with dplyr later
        colNames <- paste0("Group", 1:.Object@numberGroups)

        # Symbols needed to split the DF by dplyr
        dplyrGroup <- c(as.symbol("Effect"), lapply(colNames, as.symbol))

        colnames(profilesDE) <- colnames(profilesNonDE) <- c("ID", "DE", colNames)

        # Merge DE + nonDE for easier access
        # profilesExpr <- rbind(profilesDE, profilesNonDE)
        profilesAll <- rbind(profilesDE, profilesNonDE)

        # Create a profile table for each regulator
        # Note: prevent harmless warning about converting factors to character vector
        suppressWarnings(profilesReg <- sapply(.Object@simulators[names(.Object@simulators) != 'SimRNAseq'], function(sim) {
            message(sprintf("Generating simulation settings for %s.", sim@name))

            # Select genes ID regulator ID from all genes considered
            # TODO: use ALL genes available on the association list (current) or only those PRESENT on RNA-seq?
            # regTable <- IDfromGenes(sim, sampleExpr, simplify = FALSE) %>% filter(ID %in% rownames(sim@data))
            regTable <- IDtoGenes(sim, rownames(sim@data), simplify = FALSE)# %>% filter(ID %in% rownames(sim@data))

            # Skip the process
            if (! nrow(regTable))
                stop(sprintf("No genes were retrieved from the association list for %s, be sure to provide the correct list."))

            # Randomly assign an effect to every regulator from the available
            # options defined on every simulator.
            regTable$Effect <- sample(
                sim@regulatorEffect,
                size = nrow(regTable),
                replace = TRUE
            )

            # Set non-DE expressed or non-expressed genes to <NA>
            # TODO: same treatment for expressed/non-expressed genes?
            regTable <- dplyr::mutate(regTable, Effect = replace(Effect, ! Gene %in% sampleDE, NA))

            # Allow simulator class to adjust the generated profiles
            # (i.e. methyation using blocks)
            regTable <- adjustProfiles(sim, .Object, regTable, step = "Effect")

            # Select all (hence the double check) rows duplicated, only for DE genes.
            # Duplicated ID means that the regulator is linked to more than one gene.
            # regDups <- duplicated(regTable$Gene) | duplicated(regTable$Gene, fromLast=TRUE)
            regDupsDE <-  (duplicated(regTable$ID) |
                            duplicated(regTable$ID, fromLast=TRUE)) &
                            (regTable$Gene %in% sampleDE)

            # For every regulator with more than one regulated gene, decide which
            # class (combination of profiles among conditions) to regulate.
            message("- Linking to gene classes...")

            if (any(regDupsDE)) {
                regTable[regDupsDE, ] <- dplyr::group_by(regTable[regDupsDE, ], ID) %>% dplyr::do({
                    # Select classes of the genes associated with the regulator
                    # regGenes <- dplyr::semi_join(profilesDE, ., by = c("ID" = "Gene")) %>% dplyr::select(-DE)
                    # Note: dplyr renames the unneeded ID column of "." to ID.y
                    regGenes <- dplyr::inner_join(profilesDE, ., by = c("ID" = "Gene")) %>% dplyr::select(-DE, -ID.y)

                    # Split by class and count number of genes on each one.
                    #
                    # I.e. 3 conditions:
                    # induction - induction - repression: X rows
                    # induction - induction - induction: Y rows
                    classSplit <-
                        dplyr::group_by_(regGenes, .dots = dplyrGroup)

                    # Select genes of the majoritary class, or one at random in case
                    # of tie.
                    #
                    # Force sampling options to character to prevent a seq when there
                    # is only a result.
                    groupSizes <- dplyr::group_size(classSplit)
                    tieIndexes <- as.character(which(groupSizes == max(groupSizes)))

                    selGenes <- dplyr::slice(regGenes,
                                             which(dplyr::group_indices(classSplit) == sample(tieIndexes, 1)))$ID

                    # Set to <NA> those genes not included in the primary class
                    return(dplyr::mutate(., Effect = replace(Effect, ! Gene %in% selGenes, NA)))
                })
            }

            message("- Adjusting regulator effect")

            # Expression profiles of all expressed genes
            # condTable <- inner_join(profilesExpr, regTable, by = "Gene")
            # condTable <- dplyr::left_join(regTable, profilesExpr, by = c("Gene" = "ID")) %>% dplyr::ungroup()

            # TODO: select ALL genes (expr + non-expr) or not?
            condTable <- dplyr::left_join(regTable, profilesAll, by = c("Gene" = "ID")) %>% dplyr::ungroup()

            # TODO: dplyr version 0.5 has a bug, returning 0 columns when no one satifies
            # a 'regex' criteria (starts_with). Prevent this by selecting only those when they
            # are truly present.
            if (any(startsWith(colnames(condTable), "Keep."))) {
                subCondTable <- dplyr::select(condTable, ID, Gene, Effect, dplyr::starts_with("Keep."))
            } else {
                subCondTable <- dplyr::select(condTable, ID, Gene, Effect)
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
                            ifelse( # TODO: change this ifelse for replace?
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

        # Expressed genes profiles (DE + NonDE)
        profilesReg$SimRNAseq <- profilesAll#profilesExpr

        # In case there are no times in exp design, change the initial count values between
        # samples on DE genes. The associated regulators must change as well,
        # depending on the linked effect.
        if (.Object@numberGroups > 1) {
            # sampleDE
            # profileSubset <- dplyr::filter(profilesExpr, ID %in% sampleDE)
            # dataSubset <- object@data[as.character(profileSubset$ID), ]

            # Determine the rows with the same condition on all columns
            # Change only base count values when the profiles are flat on all conditions,
            # otherwise they are already DE genes.
            sameCond <- apply(apply(dplyr::select(profilesDE, dplyr::starts_with("Group")), 2, '==', 'flat'), 1, all)

            # NULL by default
            profilesReg$FlatGroups <- NULL

            if ((nSameCount <- sum(sameCond))) {
                message(sprintf("Creating settings to change count values on %d DE genes with the same flat profile on all groups.", nSameCount))

                # Select ID, Group1-GroupN
                sameCondDE <- t(apply(profilesDE[sameCond, -c(2)], 1, function(geneRow) {
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

                profilesReg$FlatGroups <- data.frame(sameCondDE, stringsAsFactors = FALSE)
            }
        }

        # Resolve references of profiles to proper values
        message("Finishing generation of configuration settings.")

        nt <- length(.Object@times)

        x <- c(0:(nt - 1))

        a1 <- 0
        b1 <- 1 / x[nt]

        a2 <- 0
        b2 <- 4 * (x[nt]) / (x[nt] * x[nt])
        c2 <- -4 / (x[nt] * x[nt])

        # Create a list containing all possible options as tags: 0, a1, -a1, b1, -b1, ...
        # TODO: add a new slot to provide a list of parameters instead of a predefined one?
        replace.values <- c("0" = 0, unlist(lapply(c('a1', 'b1', 'a2', 'b2', 'c2'), function(tag)
            setNames(c(get(tag), - get(tag)), c(tag, paste0('-', tag)))
        )))

        # Convert references to real values
        simProfiles <-
            sapply(.Object@profiles, function(p)
                replace.values[as.character(p)] %*% rbind(c(rep(1, nt)), x, x * x), simplify = FALSE)


        .Object@simSettings <- list(
            geneProfiles = profilesReg,

            geneSamples = list(
                DE = sampleDE,
                nonDE = sampleNonDE#,
                # expr = sampleExpr,
                # nonExpr = sampleNonExpr
            ),

            profiles = simProfiles,

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

#' @rdname simulate-methods
#' @aliases simulate,Simulation
setMethod("simulate", signature="Simulation", function(object) {
    object@simulators <- sapply(object@simulators, simulate, object)

    return(object)
})

setMethod("show", signature="Simulation", function(object) {
    # TO DO: temp...
    simSettings(object)
})

setMethod("simSettings", signature="Simulation", function(object) {
    cat(sprintf("Simulation settings of class %s:\n", class(object)))
    cat(sprintf("- Default depth: %d\n", object@depth))
    cat(sprintf("- Total genes: %d\n", object@totalGenes))
    # cat(sprintf("- Expressed genes: %d\n", object@exprGenes))
    cat(sprintf("- Dif. expressed genes: %d\n", object@diffGenes))
    cat(sprintf("- Replicates: %d\n", object@numberReps))
    cat(sprintf("- Factor levels (groups): %d\n", object@numberGroups))
    cat(sprintf("- Time vector length: %d\n", length(object@times)))

    if (is.declared(object@simSettings)) {
        # Print anything else
    }
})

setValidity("Simulation", function(object) {
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
    simuInitialized <- sapply(object@simulators, inherits, what = "Simulator")

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
        if ('SimMethylseq' %in% names(dataProvided) && dataProvided['SimMethylseq'] == TRUE)
            errors <- c(errors, "Currently there is no support for including custom methylation data. You need to provide the 'idToGene' slot, which will be used to retrieve the number and location of the CpGs.")

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

        # Check that every association list contains all genes inside
        # TODO: skip check of containing all genes, as only those in commmon will be selected.
        assocCheck <-
            sapply(object@simulators[! names(object@simulators) %in% assocExcluded], function(sim, genes) {
                return(all(colnames(sim$idToGene) %in% c('ID', 'Gene')))# &&
                           # all(genes %in% sim$idToGene$Gene))
                           # all(sim$idToGene$Gene %in% genes))
            }, genes = rownames(object@simulators[['SimRNAseq']]$data))

        if (! all(dataCheck, assocCheck))
            errors <- c(errors, "The initial sample data must have only 1 column with counts, and all association lists must contain genes present on RNA-seq data.")
    }

    return(if(length(errors)) errors else TRUE)
})
