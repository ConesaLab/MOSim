#' @include Simulator.R SimulatorRegion.R simulate_WGBS_functions.R
#' @import HiddenMarkov
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors queryHits
#' @importFrom rlang .data
NULL

#
# Methyl-seq simulator
#
# Based on simulate_WGBS.R script from WGBSSuite V 0.3 (owen.rackham@imperial.ac.uk)
#
setMethod("initialize", signature="SimMethylseq", function(.Object, idToGene, totalFeatures = nrow(idToGene), ...) {

    # Due to package size restraints, default data has a loc DF used as idToGene for methylation.
    if (ncol(idToGene) > 2) {
        .Object@locs <- idToGene %>% dplyr::mutate(ID = paste(.data$chr, .data$start, .data$end, sep = "_"))

        # Reconstruct original idToGene
        idToGene <- .Object@locs %>% dplyr::select(.data$ID, .data$Gene)
    }

    # Restrict the number of locations to totalFeatures
    # (Methyl-seq lacks "data")
    if (nrow(idToGene) > totalFeatures) {
        idToGene <- dplyr::sample_n(idToGene, totalFeatures)
    }

    .Object <- callNextMethod(.Object, idToGene = idToGene, totalFeatures = totalFeatures, ...)

    # Filter those CHR having less than 2 observations
    .Object@locs <- dplyr::semi_join(.Object@locs, dplyr::count(.Object@locs, .data$chr)
                                     %>% dplyr::filter(n > 1),
                                     by = "chr") %>%
        dplyr::filter(!is.na(.data$chr))

    # Override WGBS function to consider the CpG coordinates instead of the index
    find_adjusted_blocks <- function(a, positions) {
        l <- length(a)
        coords <- list()
        index <- 1
        state <- a[1]
        on <- positions[1]
        off <- 0
        for (i in 2:l) {
            if (a[i - 1] != a[i]) {
                if (on == 0) {
                    on <- positions[i]
                } else{
                    off <- positions[i]
                    coords[[index]] <- c(on, off, a[i], off - on)
                    index <- index + 1
                    on <- positions[i]
                    off <- 0
                }
            } else if (i == l) {
                # Add last block if last part was not executed
                off <- positions[i] + 1
                coords[[index]] <- c(on, off, a[i], off - on)
            }
        }

        coords_mat <- do.call(rbind, coords)

        return(coords_mat)
    }


    # This simulator needs to previously create the status blocks to modify
    # the simulation settings.
    #
    # Repeat for every chromosome
    .Object@WGBSparams$chrBlocks <- sapply(setNames(nm = as.character(sort(unique(.Object@locs$chr)))), function(chr) {
        message("Creating methylation state blocks for chr ", chr)

        # Order by ascendent position
        regions <- dplyr::arrange(.Object@locs[.Object@locs$chr == chr, ,drop=FALSE], .data$start)

        #simulate the state transition based on the location of the CpGs
        a <- simulate_state_transition((nrow(regions) * .Object@WGBSparams$m),
                                       c(0.01, 0.99, 0.08, 0.99),
                                       as.numeric(regions$start), 0.5,
                                       .Object@WGBSparams$transition_size)

        #extract the positions of the blocks
        state_blocks <- find_adjusted_blocks(a, regions$start)

        return(list(
            'a' = a,
            'state' = state_blocks,
            'chr' = chr
            ))
    }, simplify = 'array')

    # TODO: keep the provided data instead of forcing the default?
    .Object@data <- data.frame(Counts=rep(0, length(unique(.Object@idToGene$ID))),
                               row.names = unique(.Object@idToGene$ID), stringsAsFactors = FALSE)

    return(.Object)
})

setMethod("initializeData", signature="SimMethylseq", function(object, simulation) {
    #############################################################
    # Code from function generate_sim_set [simulate_WGBS.R]
    # Version 0.3 (01 June 2015)
    #############################################################
    object@WGBSparams$number_of_replicas <- simulation@numberReps
    object@WGBSparams$number_of_samples <- simulation@numberGroups

    # TODO: check if M values will be used?
    object@min <- 0
    object@max <- 1

    make_differential_fixed <- function(state_blocks,percentage,a){
        l<-length(a)
        states<-unique(state_blocks[,3])
        n <-length(states)
        diff_meth <- matrix(data=0,nrow=1,ncol=l)
        for(i in seq_len(n)){
            blocks <- state_blocks[state_blocks[,3]==i,]
            # MOD: sum only if it's matrix
            block_length <- if(is.null(dim(blocks))) blocks[4] else sum(blocks[,4])
            cutoff_length <- percentage[i]*block_length
            so_far <- 0
            while((so_far <= cutoff_length)&&!is.null(dim(blocks)[1])){
                picked<- sample(seq_len(dim(blocks)[1]), 1)
                if(so_far+blocks[picked,4] < cutoff_length){
                    so_far <- so_far + blocks[picked,4]
                    diff_meth[blocks[picked,1]:blocks[picked,2]]<-1
                }
                blocks<- blocks[-(picked),]
            }
        }
        return(diff_meth)
    }

    # Fixed function adding the last block
    find_the_blocks <- function(a){
        l<-length(a)
        coords<-list()
        index<-1
        state<-a[1]
        on<-1
        off<-0
        for(i in 2:l){
            if(a[i-1] != a[i]){
                if(on == 0){
                    on<-i
                }else{
                    off<-i
                    coords[[index]]<-c(on,off,a[i],off-on)
                    index<-index+1
                    on<-i
                    off<-0
                }
            } else if (i == l) {
                # Add last block if last part was not executed
                off <- i
                coords[[index]] <- c(on, off, a[i], off - on)
            }
        }

        coords_mat <- do.call(rbind,coords)

        return(coords_mat)
    }

    # Load methylation simulation settings
    methProfiles <- simulation@simSettings$geneProfiles[[class(object)]]

    # DE genes identifiers
    # geneSamples
    genesDE <- simulation@simSettings$featureSamples$SimRNAseq$DE

    # Settings for flat on all groups
    flatProfiles <- simulation@simSettings$geneProfiles$FlatGroups$SimRNAseq

    return(with(object@WGBSparams, {

        methData <- vector("list", simulation@numberGroups)
        regTable <- list()

        # Select a random number of groups to be non-modified
        non_mod_groups <- stats::runif(runif(1, 1, simulation@numberGroups - 1), 1, simulation@numberGroups)

        # Repeat process for every chromosome
        for (chr in as.character(unique(object@locs$chr))) {
            message("Simulating methylation data from CpG locations of chr ", chr, " (", distType, ")")

            regions <- object@locs[object@locs$chr == chr, ]

            a <- chrBlocks[['a', chr]]
            state_blocks <- find_the_blocks(a)
            n <- length(a)
            locs <- regions$start
            # state_blocks <- chrBlocks[['state', chr]]

            # Select blocks of the current chromosome that contain at least
            # 1 DE gene and is affected by methylation
            # diffBlocks <- subset(methProfiles, grepl(paste0("Block", chr), Block) &
            #            ! is.na(Effect) & Gene %in% genesDE, 'Block')
            chrMethProfiles <- methProfiles[methProfiles$ID %in% rownames(regions), ]

            # Select the regions that affect at least one DE gene
            diffRegions <- subset(chrMethProfiles, subset = ! is.na(Effect) & Gene %in% genesDE)

            # Select all the regions contained on the previously selected blocks
            diffBlocks <- subset(chrMethProfiles, subset = Block %in% diffRegions$Block)

            # indexBlock <- dplyr::filter(chrMethProfiles, Block %in% diffBlocks) %>%
            #     select(ID, Block) %>% distinct()
            #
            # chrMethProfiles <- dplyr::filter(methProfiles, ID %in% rownames(regions), ! is.na(Effect), Gene %in% genesDE)

            # Select the regions in the profile matrix in the same order as they
            # are here, and check if the block is inside the previously selected
            # list of blocks, transforming the boolean to 0/1 values.
            diff_methed <- matrix(data = as.numeric(rownames(regions) %in% diffBlocks$ID), nrow = 1, ncol = n)

            # Select all

            # The number of the block "Block<Chr>.number" equals to the row
            # number of the state block matrix
            # diffBlocks <- as.numeric(stri_extract_last_regex(diffBlocks, '\\d+'))

            # Emulate the ouput of the function "make_differential", a matrix with
            # 1 row and N (number of CpG) columns, with 0 for non-modified and 1
            # for modified.
            # diff_methed <- matrix(data=0, nrow=1, ncol=n)

            #set the percentage of DM in each block type
            # percs_for_diff <- c(0,0,0,0.5)

            #update the blocks to be differentially methylated
            # diff_methed <- make_differential(state_blocks, percs_for_diff, a)

            #create the simulated reads methylated/unmethylated at each CpG, at the moment this is hard coded to be 3 replicates of each type
            #The phase diff param control how different the methylation is in the differentially methylated regions.
            errors <-
                list(error_m, ((error_d + error_m) / 2), error_d, ((error_d + error_m) / 2))
            means <-
                list(mean_m, ((mean_d + mean_m) / 2), mean_d, ((mean_d + mean_m) / 2))


            # Generate data for every sample and every rep
            repData <- list()
            # Use format to prevent scientific notation on rownames
            # repDataNames <- regionspaste(chr, format(regions$start, scientific = FALSE),
                                  # format(regions$end, scientific = FALSE), sep="_") #regionNames(object, chrNumber=i, start=locs)

            # Make all groups identical and later choose random values like in the
            # other simulators.
            diffPhase <- rep(0, simulation@numberGroups)

            prob_d_factor <- lapply(seq_len(simulation@numberGroups), function(f) {
                return(list(
                    d=hypo_hyper_diffs(prob_d, diffPhase[f], diff_methed, balance),
                    m=hypo_hyper_diffs(prob_m, diffPhase[f], diff_methed, balance)
                ))
            })

            if (! is.null(flatProfiles)) {
                flatMethProfiles <- dplyr::select(chrMethProfiles, ID, Gene, Block, Effect) %>%
                    dplyr::inner_join(flatProfiles, by = c("Gene" = "ID")) %>%
                    dplyr::select(Block, Effect, dplyr::starts_with("Group")) %>%
                    dplyr::filter(!is.na(Effect)) %>%
                    dplyr::distinct()

                if (nrow(flatMethProfiles)) {
                    # By default the methylation prob is high, so for reflecting
                    # a change we lower it on the flat groups instead of increasing
                    # it on the proper group (read below).
                    prob_m_mod <- prob_m - max(phase_diff)
                    prob_d_mod <- prob_d + max(phase_diff)

                    # Repeat for every block
                    for (i in seq_len(nrow(flatMethProfiles))) {
                        blockProfiles <- unlist(flatMethProfiles[i, ])

                        # ID of the block
                        blockId <- blockProfiles[1]

                        # Effect of the block
                        blockEffect <- blockProfiles[2]

                        # Expression profiles
                        blockExpr <- blockProfiles[seq(3, length.out = simulation@numberGroups)]

                        # Retrieve the groups used as reference
                        blockFlat <- which(blockExpr == 'flat')

                        # Modified groups
                        blockNonFlat <- which(blockExpr != 'flat')

                        # Associated effect: enhancer or repressor
                        exprEffect <- blockExpr[blockNonFlat][1]

                        # Positions of the probs vector to update
                        # blockIndexes <- match(blockId, chrMethProfiles$Block)
                        blockIndexes <- match(unique(subset(chrMethProfiles, Block == blockId)$ID), rownames(regions))

                        # Check every RNA-seq group for increased (enhancer) or
                        # decreased (repressor) expression, altering the probs of
                        # methylated/non-methylated reads depending on the effect
                        # associated to the block:
                        #   - enhancer: will increase methylation probs with enhancer,
                        # and decrease them with repressor.
                        #   - repressor: decrease methylation probs with enhancer, and
                        # increase them with repressor.
                        #
                        # If the effect of the regulator is complementary to the
                        # change on the expression data, increase methylation probs
                        # and reduce non-methylated.
                        if (blockEffect == exprEffect) {
                            # Modify flat (reference) groups
                            changeBlocks <- blockFlat
                        } else {
                            # Modify non-flat groups
                            changeBlocks <- blockNonFlat
                        }

                        # Repeat for each affected group
                        for (f in changeBlocks) {
                            prob_d_factor[[f]]$d[blockIndexes] <- prob_d_mod
                            prob_d_factor[[f]]$m[blockIndexes] <- prob_m_mod
                        }
                    }
                }
            }

            # Get one "non-mod" data
            # prob_nonmod <- prob_d_factor[[which(diffPhase == 0)[1]]]
            # prob_nonmod <- prob_d_factor[[which.max(diffPhase == 0)]]
            prob_nonmod <- prob_d_factor[[1]]

            probs <- vector("list", simulation@numberGroups)

            for (f in seq_len(simulation@numberGroups)) {
                # Diff. sites for every factor
                prob_f <- prob_d_factor[[f]]

                # Non-mod
                # if (diffPhase[f] == 0) {
                #     prob_f <- list(prob_f$m, ((prob_nonmod$m + prob_nonmod$d) / 2), prob_f$d, ((prob_f$m + prob_f$d) / 2))
                # } else {
                #     prob_f <- list(prob_f$m, ((prob_f$m + prob_f$d) / 2), prob_f$d, ((prob_f$m + prob_f$d) / 2))
                # }

                prob_f <- list(prob_f$m, ((prob_f$m + prob_f$d) / 2), prob_f$d, ((prob_f$m + prob_f$d) / 2))

                probs[[f]] <- prob_f

                # New matrix for every sample
                repData <- matrix(data=0, nrow=simulation@numberReps, ncol=n)

                for (r in seq_len(simulation@numberReps)) {
                    if (distType == 'binomial') {
                        sRepData <-
                            generate_replicat_methyl_bin_data(a, probs, means, 0, errors, locs, f, output_path)
                    } else if (distType == 'truncated') {
                        sRepData <-
                            generate_replicat_methyl_truncated_nbin_data(a, probs, means, 0, errors, locs, f, 20, output_path)
                    } else{
                        sRepData <-
                            generate_replicat_methyl_nbin_data_model_3(a, probs, means, 0, errors, locs, f, 30)
                    }
                    # Columns:
                    #         V1: location in base pairs
                    #         V2: differentially methylated flag
                    #
                    # Blocks of 4 columns for each replica as follows:
                    #
                    #         Vn: Number of methylated reads
                    #         Vn+1: Number of de-methylated reads
                    #         Vn+2: Total number of reads
                    #         Vn+3: Proportion of methylated vs de-methylated reads
                    sRepData <- sRepData[returnColumn,]

                    repData[r,] <- sRepData
                }

                repData <- t(repData)

                # Add rows to the factor data
                rownames(repData) <- regions$ID #repDataNames

                # Repeating for every group
                methData[[f]] <- rbind(methData[[f]], repData)
            }
        }

        # Save random values for every group to be used later on simulate function.
        # This needs to be done on a "block basis".
        uniqueBlocks <- unique(methProfiles$Block)

        object@WGBSparams$randomCounts <- lapply(seq_len(simulation@numberGroups), function(group) {
            # Min/max values
            simRange <- range(methData[[group]])

            # Create random values for every block
            randomValues <- stats::runif(length(uniqueBlocks), min = simRange[1], max = simRange[2])

            return(dplyr::inner_join(methProfiles, data.frame(Block=uniqueBlocks,
                                                              M=randomValues,
                                                              stringsAsFactors = FALSE),
                                     by = c("Block" = "Block")) %>%
                       dplyr::distinct(ID, M)) %>%
                       dplyr::mutate(M = jitter(M))
        })

        # Merge simulated data from every group
        methData <- do.call(cbind, methData)

        object@WGBSparams$blocks <- chrBlocks
        object@data <- methData

        # Pattern "Counts.Group" required for post-simulation. Later it will
        # be correctly renamed.
        colnames(object@data) <- paste('Counts.Group', seq_len(ncol(object@data)))

        return(object)
    }))
})


setMethod("simulateParams", signature="SimMethylseq", function(object, simulation, counts, profiles, group, ids) {
    # Blocks
    randomValues <- object@WGBSparams$randomCounts[[group]]
    randomValues <- randomValues[match(ids, randomValues$ID), ]

    # Order [m, M]
    m <- pmin(counts, randomValues$M)
    M <- pmax(counts, randomValues$M)

    # Do not add extra noise
    # noiseValues <- matrix(0, ncol = length(simulation@times), nrow = length(counts))
    # Noise with fixed parameters and function (the value can only go from 0 to 1)
    # TODO: remove this extra noise?
    noiseValues <-
        replicate(length(simulation@times),
                  do.call(
                      object@noiseFunction,
                      list(n = length(counts), sd = 0.03)
                  ))

    repressionMean <- stats::runif(length(M), min = (M+m)/2, max = M)
    inductionMean <- stats::runif(length(m), min = m, max = (M+m)/2)

    return(list(
        'randomCounts' = randomValues$M,
        'noiseValues' = noiseValues,
        'm' = m,
        'M' = M,
        'repression.Mean' = repressionMean,
        'induction.Mean' = inductionMean,
        'mean.counts' = ifelse(grepl('repression', profiles),(repressionMean + m)/2, (inductionMean + M)/2),
        'counts' = counts
    ))
})

setMethod("postSimulation", signature="SimMethylseq", function(object, simulation) {
    # Limit the range of proportions to 0-1
    object@simData[object@simData < 0] <- 0
    object@simData[object@simData > 1] <- 1

    # Add M-values
    if (object@Mvalues) {
        # methLevel & rowData comes from BiSeq package?
        betaThreshold <- object@betaThreshold

        # Convert to numeric with data.matrix
        beta <- pmin(pmax(data.matrix(object@simData), betaThreshold), 1 - betaThreshold)

        colnames(beta) <- colnames(object@simData)

        rownames(beta) <- rownames(object@simData)

        # beta <- beta[rowVars(beta, na.rm=T)!=0,]
        M <- log2(beta/(1-beta))

        # Return M values
        # TODO: decide which data to return
        # object@simData <- list(
        #     raw=object@simData,
        #     beta=beta,
        #     M=M
        # )

        object@simData <- as.data.frame(M)
    }

    # Reorder the columns to respect the order timeX.<reps>
    numberTimes <- length(simulation@times)

    columnReorder <- unlist(
        lapply(seq_len(numberTimes),
           seq,
           by = numberTimes,
           length.out = simulation@numberReps)
    )

    # Add extra groups
    if (simulation@numberGroups > 1) {
        for (groupNumber in seq(from = 2, to = simulation@numberGroups)) {
            columnReorder <- c(columnReorder, columnReorder +
                                   (numberTimes * simulation@numberReps)*(groupNumber - 1))
        }
    }

    object@simData <- object@simData[, columnReorder]

    object <- callNextMethod(object, simulation)

    return(object)
})

setMethod("adjustProfiles", signature="SimMethylseq", function(object, simulation, profiles, step) {

    message("Adjusting methylation profiles")

    availableEffects <- object@regulatorEffect[grep("NE", names(object@regulatorEffect), invert = TRUE)]

    switch(step,
           Effect = return(dplyr::left_join(profiles, unique(dplyr::select(object@locs, .data$chr, .data$start, .data$end, .data$ID)), by = c("ID" = "ID")) %>%
                          dplyr::group_by(.data$chr) %>% dplyr::do({
                              # Chromosome name
                              chrName <- .$chr[[1]]

                              message(sprintf("Generating block methylation data for chromosome %s.", chrName))

                              # State blocks for the current chromosome
                              # Note: almost half of the blocks generated by WGBS consist of a single CpG island repeated on another block
                              # so we filter them.
                              state_blocks <- data.frame(object@WGBSparams$chrBlocks[['state', as.character(chrName)]],
                                                         stringsAsFactors = FALSE)

                              stateIRanges <- IRanges::IRanges(start = state_blocks[,1], end = state_blocks[,2] - 1)
                              chrIRanges <- IRanges::IRanges(start = .$start, end = .$end)

                              # The query column in overlaps will correspond to the block number.
                              blockOverlaps <- IRanges::findOverlaps(stateIRanges, chrIRanges)

                              blockNames <- S4Vectors::queryHits(blockOverlaps)

                              if (length(blockNames) < nrow(.)) {
                                  blockNames <- c(blockNames, seq(from = max(blockNames) + 1,
                                                                  length.out = nrow(.) - length(blockNames)))
                              }

                              within(., {
                                  # Rename ID column to CpG
                                  Keep.CpG <- .$ID
                                  # Different name for every block
                                  # Overwrite ID to allow grouping
                                  ID <- paste('Block', chrName, blockNames, sep = ".")
                                  # Assign a random effect per block
                                  # Keep the already set NA (non DE genes)
                                  Effect <- ifelse(is.na(.$Effect), NA, sample(names(availableEffects), size = 1,
                                                                                   prob = as.numeric(availableEffects)))
                              })
                          }) %>% dplyr::ungroup()),
            Groups = return({
                    # Retrieve true group effect for each block
                    dplyr::select(profiles, .data$ID, .data$Effect, dplyr::starts_with("Group")) %>%
                        dplyr::filter(!is.na(.data$Effect)) %>%
                        dplyr::select(-.data$Effect) %>%
                        dplyr::distinct_() %>%
                        dplyr::bind_rows(dplyr::select(profiles, .data$ID, .data$Effect, dplyr::starts_with("Group")) %>%
                                             dplyr::filter(is.na(.data$Effect)) %>%
                                             dplyr::select(-.data$Effect) %>%
                                             dplyr::distinct_()
                        )%>%
                    # Join with the full profile table
                    dplyr::left_join(dplyr::select(profiles, - dplyr::starts_with("Group")), by = c("ID" = "ID")) %>%
                    dplyr::mutate(Block = .data$ID, ID = .data$Keep.CpG) %>%
                    dplyr::select(.data$ID, .data$Gene, .data$Block, .data$Effect,
                                  dplyr::starts_with("Effect.Group"), dplyr::starts_with("Group"), dplyr::starts_with("Tmax.")) %>%
                    dplyr::distinct_()
            })
    )

    # By default return the same
    return(profiles)
})
