#' @include Simulation.R
#' @import ggplot2
#' @importFrom dplyr %>%
#' @importFrom rlang := !! .data
#' @importFrom stats setNames
NULL

# Avoid harmless note with R CMD check
if(getRversion() >= "2.15.1")  utils::globalVariables(c(".", "n", "sampleData"))

#' MOSim
#'
#' Multiomics simulation package.
#'
#'
#' @docType package
#' @name MOSim-package
#'
NULL
#> NULL

# The default 'mosim' values should be identical to the prototype values of
# simulation class. To avoid warnings checking the package we use
defaultPrototypeValues <- getClass("Simulation")@prototype

#' mosim
#'
#' Performs a multiomic simulation by chaining two actions: 1) Creating the
#' "Simulation" class with the provided params. 2) Calling "simulate" method on
#' the initialized object.
#'
#' @param omics Character vector containing the names of the omics to simulate,
#'   which can be "RNA-seq", "miRNA-seq", "DNase-seq", "ChIP-seq" or
#'   "Methyl-seq"  (e.g. c("RNA-seq", "miRNA-seq")). It can also be a list with
#'   the omic names as names and their options as values, but we recommend to
#'   use the argument omicSim to provide the options to simulated each omic.
#' @param omicsOptions List containing the options to simulate each omic. We
#'   recommend to apply the helper method \link{omicSim} to create this list in
#'   a friendly way, and the function \link{omicData} to provide custom data
#'   (see the related sections for more information). Each omic may have
#'   different configuration parameters, but the common ones are:
#'   \describe{
#'       \item{simuData/idToGene}{Seed sample and association tables for
#'       regulatory omics. The helper function \link{omicData} should be used to
#'       provide this information (see the following section).}
#'       \item{regulatorEffect}{For regulatory omics. List containing the
#'       percentage of effect types (repressor, activator or no effect) over the
#'       total number of regulators. See vignette for more information.}
#'       \item{totalFeatures}{Number of features to simulate. By default, the
#'       total number of features in the seed dataset.}
#'       \item{depth}{Sequencing depth in millions of reads. If not provided,
#'       it takes the global parameter passed to \link{mosim} function.}
#'       \item{replicateParams}{List with parameters \emph{a} and \emph{b} for
#'       adjusting the variability in the generation of replicates using the
#'       negative binomial. See vignette for more information.}
#'    }
#' @param diffGenes  Number of differentially expressed genes to simulate, given
#'   in percentage (0 - 1) or in absolute number (> 1). By default 0.15
#' @param numberReps Number of replicates per experimetal condition (and time
#'   point, if time series are to be generated). By default 3.
#' @param numberGroups Number of experimental groups or conditions to simulate.
#' @param times Vector of time points to consider in the experimental design.
#' @param TFtoGene A logical value indicating if default transcription factors
#'   data should be used (TRUE) or not (FALSE), or a 3 column data frame
#'   containing custom associations. By default FALSE.
#' @param depth Sequencing depth in millions of reads.
#' @param profileProbs Numeric vector with the probabilities to assign each of
#'   the patterns. Defaults to 0.2 for each.
#' @param minMaxFC Numeric vector of length 2 with minimum and maximum fold-change
#'   for differentially expressed features, respectively.
#'
#' @return Instance of class "Simulation" containing the multiomic simulation
#'   data.
#' @export
#'
#' @examples
#'
#'  moSimulation <- mosim(
#'      omics = c("RNA-seq"),
#'      numberReps = 3,
#'      times = c(0, 2, 6, 12, 24)
#'  )
#'
#'  # Retrieve simulated count matrix for RNA-seq
#'  dataRNAseq <- omicResults(moSimulation, "RNA-seq")
#'
#'
mosim <-
    function(omics,
             omicsOptions,
             diffGenes,
             numberReps,
             numberGroups,
             times,
             depth,
             profileProbs,
             minMaxFC,
             TFtoGene
    ) {
    # Check for mandatory parameters
    if (missing(omics))
        stop("You must provide the list of omics to simulate.")

    # Params to initialize simulation instance
    mosimCall <- as.list(match.call(expand.dots = FALSE))[-1]
    simParams <- mapply(function(argName, argValue) {
        if (is.symbol(argValue))
            return(get(argName))

        return(argValue)
    }, names(mosimCall), mosimCall, SIMPLIFY = FALSE, USE.NAMES = TRUE)

    # Select unique names
    omics <- if (is.list(omics)) omics[unique(names(omics))] else unique(omics)

    # 'omics' parameter alias of 'simulators'
    simParams$simulators <- omics

    # If it is a plain vector, transform it to a list of empty lists
    if (! is.list(simParams$simulators)) {
        simParams$simulators <- setNames(rep(list(list()), length(simParams$simulators)),
                                         simParams$simulators)
    }

    if (! missing(omicsOptions)) {
        if (! all(names(omicsOptions) %in% names(simParams$simulators)))
            stop("There are keys in 'omicOptions' not present on 'omics' list.")

        # Set options (not required)
        for (omicName in names(omicsOptions)) {
            omicParams <- omicsOptions[[omicName]]

            # Transform the "regulatorEffect" parameter to change "activator" to "enhancer"
            # effect.
            if (is.element("regulatorEffect", names(omicParams))) {
                names(omicParams$regulatorEffect) <- gsub("activator", "enhancer", names(omicParams$regulatorEffect))
            }

            # Modify the param list or the slots directly
            if (! inherits(simParams$simulators[[omicName]], "Simulator")) {
                simParams$simulators[[omicName]] <- omicParams
            } else {
                # Override every slot manually
                for (slotName in names(omicParams)) {
                    slot(simParams$simulators[[omicName]], slotName) <- omicParams[[slotName]]
                }
            }
        }
    }

    # Remove keys not present in 'Simulation' class slots
    simParams <- simParams[names(simParams) %in% slotNames("Simulation")]

    oSim <- do.call(new, c("Class"="Simulation", simParams))
    oSim <- simulate(oSim)

    return(oSim)
}

#' Set customized data for an omic.
#'
#' @param omic The name of the omic to provide data.
#' @param data Data frame with the omic identifiers as row names and just one
#'   column named Counts containing numeric values used as initial sample for
#'   the simulation.
#' @param associationList Only for regulatory omics, a data frame with 2
#'   columns, the first called containing the regulator ID and the second called
#'   Gene with the gene identifier.
#'
#' @return Initialized simulation object with the given data.
#' @export
#'
#' @examples
#'
#' # Take a subset of the included dataset for illustration
#' # purposes. We could also load it from a csv file or RData,
#' # as long as we transform it to have 1 column named "Counts"
#' # and the identifiers as row names.
#' custom_rnaseq <- head(sampleData$SimRNAseq$data, 100)
#'
#' # In this case, 'custom_rnaseq' is a data frame with
#' # the structure:
#' head(custom_rnaseq)
#' ##                    Counts
#' ## ENSMUSG00000000001   6572
#' ## ENSMUSG00000000003      0
#' ## ENSMUSG00000000028   4644
#' ## ENSMUSG00000000031      8
#' ## ENSMUSG00000000037      0
#' ## ENSMUSG00000000049      0
#'
#'
#' # The helper 'omicData' returns an object with our custom data.
#' rnaseq_customdata <- omicData("RNA-seq", data = custom_rnaseq)
#'
omicData <- function(omic, data = NULL, associationList = NULL) {
    # Convert the name to the proper class name
    omicClass <- paste0("Sim", gsub("-", "", omic))

    if (is(data, "ExpressionSet")) {
        data <- as.data.frame(Biobase::exprs(data))
    }

    # Create the new instance
    omicSim <- new(omicClass,
                   "data" = data,
                   "idToGene" = associationList)

    return(omicSim)
}

#' Set the simulation settings for an omic.
#'
#' @param omic Name of the omic to set the settings.
#' @param depth Sequencing depth in millions of counts. If not provided will
#'   take the global parameter passed to mosim function.
#' @param totalFeatures Limit the number of features to simulate. By default
#'   include all present in the dataset.
#' @param regulatorEffect only for regulatory omics. Associative list containing
#'   the percentage of effects over the total number of regulator, including
#'   repressor, association and no effect (NE).
#'
#' @return A list with the appropiate structure to be given as options in mosim
#'   function.
#' @export
#'
#' @examples
#'
#' omic_list <- c("RNA-seq")
#'
#' rnaseq_options <- omicSim("RNA-seq", totalFeatures = 2500)
#'
#' # The return value is an associative list compatible with
#' # 'omicsOptions'
#' rnaseq_simulation <- mosim(omics = omic_list,
#'                            omicsOptions = rnaseq_options)
#'
omicSim <- function(omic, depth = NULL, totalFeatures = NULL, regulatorEffect = NULL) {

    paramList <- list(
        "totalFeatures" = totalFeatures,
        "regulatorEffect" = regulatorEffect,
        "depth" = depth
    )

    paramList <- setNames(list(Filter(Negate(is.null), paramList)), omic)

    return(paramList)
}


#' Retrieves the settings used in a simulation
#'
#' @param simulation A Simulation object.
#' @param omics List of omics to retrieve the settings.
#' @param association A boolean indicating if the association must also be
#'   returned for the regulators.
#' @param reverse A boolean, swap the column order in the association list in
#'   case we want to use the output directly and the program requires a
#'   different ordering.
#' @param only.linked Return only the interactions that have an effect.
#' @param prefix Logical indicating if the name of the omic should prefix the
#'   name of the regulator.
#' @param include.lagged Logical indicating if interactions with transitory
#'   profile and different minimum/maximum time point between gene and regulator
#'   should be included or not.
#'
#' @return A list containing a data frame with the settings used to simulate
#'   each of the indicated omics. If association is TRUE, it will be a list with
#'   3 keys: 'associations', 'settings' and 'regulators', with the first two
#'   keys being a list containing the information for the selected omics and the
#'   last one a global data frame giving the merged information.
#' @export
#'
#' @examples
#' \dontrun{
#' omic_list <- c("RNA-seq")
#' rnaseq_simulation <- mosim(omics = omic_list)
#'
#' # This will be a data frame with RNA-seq settings (DE flag, profiles)
#' rnaseq_settings <- omicSettings(rnaseq_simulation, "RNA-seq")
#'
#' # This will be a list containing all the simulated omics (RNA-seq
#' # and DNase-seq in this case)
#' all_settings <- omicSettings(rnaseq_simulation)
#' }
omicSettings <- function(simulation, omics = NULL, association = FALSE, reverse = FALSE, only.linked = FALSE, prefix = FALSE, include.lagged = TRUE) {
    # Select all omics by default
    if (is.null(omics)) {
        omics <- lapply(simulation@simulators, slot, name = "name")
    } else {
        names(omics) <- paste0("Sim", gsub("-", "", omics))
    }

    # Convert the name to the proper class name
    omicsClasses <- setNames(paste0("Sim", gsub("-", "", omics)), omics)

    selectedSettings <- simulation@simSettings$geneProfiles[omicsClasses]

    outputList <- setNames(selectedSettings, omics)

    # Remove Effect from settings and rename Effect.Linked to Effect
    # selectedSettings <- lapply(selectedSettings, function(x) {
    #     if (is.element("Effect", colnames(x))) {
    #         x <- dplyr::select(x, -Effect) %>%
    #                          dplyr::rename(Effect = Effect.Linked)
    #     }
    #
    #     return(x)
    # })

    # Retrieve gene settings always to calculate lagged relationships.
    geneSettings <- simulation@simSettings$geneProfiles[["SimRNAseq"]]


    filter_valid_effect <- function(df) {
        filtered_settings <- df %>%
            dplyr::filter_at(dplyr::vars(dplyr::starts_with("Effect")),
                             dplyr::any_vars(! is.na(.)))

        return(filtered_settings)
    }

    # Remove & rename effects from data frames
    replace_effect_filter <- function(df) {
        df <- dplyr::mutate_at(df, dplyr::vars(dplyr::starts_with("Effect.Group")), dplyr::funs(gsub("enhancer", "activator", .)))

        if (all(c("Group1", "Effect") %in% colnames(df))) {

            propagate_profile <- function(group_values, effects) {
                with_effect <- ! is.na(effects)

                if (any(with_effect)) {
                    # Everything with an effect should be the same
                    return(group_values[with_effect][1])
                }

                return(group_values)
            }

            add_lagged_info <- function(pipedDF) {

                replaceTmax <- function(col, fullDF) {
                    groupName <- gsub("Tmax.", "", deparse(substitute(col)))
                    profilePatterns <- fullDF %>% dplyr::pull(groupName)

                    return(ifelse(grepl("transitory", profilePatterns), col, NA))
                }

                # Replace Tmax with NA when flat profile
                pipedDF <- pipedDF %>% dplyr::mutate_at(dplyr::vars(dplyr::starts_with("Tmax.Group")), replaceTmax, pipedDF)

                mergedDF <- pipedDF %>% dplyr::left_join(geneSettings,
                                                         suffix = c(".x", ".y"), by = c("Gene" = "ID"))

                for (groupNumber in seq_len(simulation@numberGroups)) {
                    newColName <- paste0("Lagged.Group", groupNumber)

                    regTmaxCol <- rlang::sym(paste0("Tmax.Group", groupNumber, ".x"))
                    geneTmaxCol <- rlang::sym(paste0("Tmax.Group", groupNumber, ".y"))

                    laggedValue <- mergedDF %>%
                        dplyr::mutate(Lagged = !!(regTmaxCol) != !!(geneTmaxCol)) %>%
                        dplyr::pull(.data$Lagged)

                    pipedDF <- pipedDF %>% dplyr::mutate(!!newColName :=  laggedValue)
                }

                return(pipedDF)
            }

            # Assign the same effect to regulator rows.
            df <- dplyr::group_by(df, .data$ID) %>%
                dplyr::mutate_at(dplyr::vars(dplyr::starts_with("Group")),
                                 dplyr::funs(propagate_profile(., .data$Effect))) %>%
                dplyr::ungroup() %>%
                add_lagged_info()

            if (! include.lagged) {
                df <- dplyr::mutate(df, IsLagged = purrr::pmap_int(dplyr::select(df, dplyr::starts_with("Lagged.Group")), sum, na.rm=TRUE)) %>%
                    dplyr::filter(.data$IsLagged == FALSE) %>%
                    dplyr::select(-.data$IsLagged, -dplyr::starts_with("Lagged"))
            }

            if (only.linked) {
                df <- filter_valid_effect(df)
            }
        }

        return(df[, colnames(df) != "Effect"])
    }


    # Add the associations
    if (association) {

        regClasses <- omicsClasses[grep("SimRNAseq", omicsClasses, invert = TRUE)]

        associationLists <- setNames(lapply(simulation@simulators[regClasses],
                                   slot, name = "idToGene"), regClasses)

        # Include only those regulators having an effect on genes
        if (only.linked) {
            associationLists <- sapply(names(associationLists), function(x) {
                omic.settings <- filter_valid_effect(selectedSettings[[x]])

                omic.association <- associationLists[[x]] %>%
                    dplyr::filter(.data$ID %in% omic.settings$ID)

                return(omic.association)
            }, simplify = FALSE)
        }

        # Prepend the name of the omic to the regulator ID
        if (prefix) {
            associationLists <- sapply(names(associationLists), function(x) {

                prefix.settings <- associationLists[[x]] %>% dplyr::mutate(ID = paste0(x, .data$ID))

                return(prefix.settings)
            }, simplify = FALSE)
        }

        if (reverse) {
            associationLists <- lapply(associationLists, function(x) if(ncol(x) > 1) x[, c(2,1)] else x)
        }

        # Generate a global table containing the associated regulator per gene and omic
        globalTable <- do.call(rbind, lapply(regClasses, function(x) {
            x.data <- selectedSettings[[x]]
            x.settings <- outputList[[omics[[x]]]]

            # Select only the effect
            x.settings <- dplyr::filter(x.settings, !is.na(.data$Effect)) %>%
                dplyr::select(.data$ID, dplyr::starts_with("Group")) %>%
                dplyr::distinct()

            output.df <- dplyr::select(x.data, .data$Gene, .data$ID, .data$Effect, dplyr::starts_with("Effect.Group"), dplyr::starts_with("Tmax.Group")) %>%
                dplyr::mutate(Omic = x) %>% dplyr::left_join(x.settings, by = c("ID" = "ID")) %>%
                dplyr::mutate_at(dplyr::vars(dplyr::starts_with("Group")), dplyr::funs(ifelse(is.na(.), "flat", .)))

            if (only.linked) {
                linked.IDs <- dplyr::filter_at(output.df, dplyr::vars(dplyr::starts_with("Effect")), dplyr::any_vars(! is.na(.)))$ID

                output.df <- dplyr::filter(output.df, .data$ID %in% linked.IDs)
            }

            return(output.df)
        }))

        # Restore correct names
        associationLists <- setNames(associationLists, names(omicsClasses)[match(names(associationLists), omicsClasses)])

        outputList <- list(
            "association" = associationLists,
            "settings" = sapply(outputList, replace_effect_filter, simplify = FALSE, USE.NAMES = TRUE),
            "regulators" = replace_effect_filter(globalTable)
        )
    } else {
        # Replace effects on settings
        outputList <- sapply(outputList, replace_effect_filter, simplify = FALSE, USE.NAMES = TRUE)
    }

    if (length(omics) > 1 || length(outputList) > 1) {
        return(outputList)
    } else {
        return(outputList[[omics]])
    }
}

#' Retrieves the simulated data.
#'
#' @param simulation A Simulation object.
#' @param omics List of the omics to retrieve the simulated data.
#' @param format Type of object to use for returning the results
#'
#' @return A list containing an element for every omic specifiec, with the
#'   simulation data in the format indicated, or a numeric matrix with simulated
#'   data if the omic name is directly provided.
#' @export
#'
#' @examples
#'
#' omic_list <- c("RNA-seq")
#' rnaseq_simulation <- mosim(omics = omic_list)
#' #' # This will be a data frame with RNA-seq counts
#' rnaseq_simulated <- omicResults(rnaseq_simulation, "RNA-seq")
#'
#' #                    Group1.Time0.Rep1 Group1.Time0.Rep2 Group1.Time0.Rep3 ...
#' # ENSMUSG00000073155              4539              5374              5808 ...
#' # ENSMUSG00000026251                 0                 0                 0 ...
#' # ENSMUSG00000040472              2742              2714              2912 ...
#' # ENSMUSG00000021598              5256              4640              5130 ...
#' # ENSMUSG00000032348               421               348               492 ...
#' # ENSMUSG00000097226                16                14                 9 ...
#' # ENSMUSG00000027857                 0                 0                 0 ...
#' # ENSMUSG00000032081                 1                 0                 0 ...
#' # ENSMUSG00000097164               794               822               965 ...
#' # ENSMUSG00000097871                 0                 0                 0 ...
#'
omicResults <- function(simulation, omics = NULL, format = "data.frame") {
    # Select all omics by default
    if (is.null(omics)) {
        omics <- lapply(simulation@simulators, slot, name = "name")
    }

    # Convert the name to the proper class name
    omicsClasses <- paste0("Sim", gsub("-", "", omics))

    selectedSimulators <- lapply(simulation@simulators[omicsClasses], slot, name = "simData")

    formatedSimulators <- lapply(selectedSimulators, function (simuData) {
        formatedData <- switch(
            format,
            "data.frame" = data.frame(simuData, stringsAsFactors = FALSE),
            "ExpressionSet" = Biobase::ExpressionSet(assayData = simuData)
        )

        return(formatedData)
    })

    outputList <- setNames(formatedSimulators, omics)

    if (length(omics) > 1) {
        return(outputList)
    } else {
        return(outputList[[omics]])
    }
}

#' Retrieves the experimental design
#'
#' @param simulation A Simulation object
#'
#' @return A data frame containing the experimental design used to simulate the
#'   data.
#' @export
#'
#' @examples
#'
#' omic_list <- c("RNA-seq")
#' rnaseq_simulation <- mosim(omics = omic_list)
#' # This will be a data frame with RNA-seq counts
#'
#' design_matrix <- experimentalDesign(rnaseq_simulation)
#'
experimentalDesign <- function(simulation) {
    sampleNames <- colnames(simulation@simulators$SimRNAseq@simData)

    outputDF <- data.frame(
        Group = gsub("Group([0-9]+)\\..*", "\\1", sampleNames),
        Time = as.numeric(gsub(".*Time([0-9]+)\\..*", "\\1", sampleNames)),
        Rep = gsub(".*Rep([0-9]+)", "\\1", sampleNames),
        row.names = sampleNames
    )

    outputDF$Condition <- paste(outputDF$Group, outputDF$Time, sep = ".")

    return(outputDF)
}

#' Generate a plot of a feature's profile for one or two omics.
#'
#' @param simulation A Simulation object
#' @param omics Character vector of the omics to simulate.
#' @param featureIDS List containing the feature to show per omic. Must have the
#'   omics as the list names and the features as values.
#' @param drawReps Logical to enable/disable the representation of the
#'   replicates inside the plot.
#'
#' @return A ggplot2 object.
#' @export
#'
#' @examples
#' omic_list <- c("RNA-seq", "miRNA-seq")
#' rnaseq_simulation <- mosim(omics = omic_list)
#'
#' plotProfile(rnaseq_simulation,
#'     omics = c("RNA-seq", "miRNA-seq"),
#'     featureIDS = list("RNA-seq"="ENSMUSG00000007682", "miRNA-seq"="mmu-miR-320-3p")
#' )
#'
plotProfile <-
    function(simulation,
             omics,
             featureIDS,
             drawReps = FALSE) {

    numberOfReps <- simulation@numberReps

    calculateRowMeans <- function(colPattern, df, sd = FALSE) {

        colName <- gsub(".Rep", "", colPattern)

        if (! sd) {
            output <- as.data.frame(df) %>%
                dplyr::mutate(!!colName := purrr::pmap_dbl(dplyr::select(., dplyr::contains(colPattern)),
                                                           function(...) mean(c(...)))) %>%
                dplyr::select(colName)
        } else {
            meanCol <- rlang::sym(paste0(colName, ".Mean"))
            seCol <- rlang::sym(paste0(colName, ".SE"))
            sdCol <- rlang::sym(paste0(colName, ".SD"))

            std <- function(x) sd(x)/sqrt(length(x))

            output <- as.data.frame(df) %>%
                dplyr::mutate(!!meanCol := purrr::pmap_dbl(dplyr::select(., dplyr::contains(colPattern)),
                                                           function(...) mean(c(...))),
                              !!seCol := purrr::pmap_dbl(dplyr::select(., dplyr::contains(colPattern)),
                                                         function(...) std(c(...))),
                              !!sdCol := purrr::pmap_dbl(dplyr::select(., dplyr::contains(colPattern)),
                                                         function(...) sd(c(...)))) %>%
                dplyr::select(!!meanCol, !!seCol, !!sdCol)
        }
        rownames(output) <- rownames(df)

        return(output)
    }

    calculateRepMeans <- function(omicDF) {
        column_names_simu <- unique(stringr::str_extract(colnames(omicDF), "Group[0-9]+\\.Time[0-9]+\\.Rep"))
        data_means_simu <- do.call(cbind, lapply(column_names_simu, calculateRowMeans, omicDF, sd = TRUE))

        return(data_means_simu)
    }

    getOmicFeaturesDF <- function(omicName, feature) {
        simuData <- omicResults(simulation, omicName)
        simuSettings <- omicSettings(simulation, omicName)

        if (! feature %in% simuSettings$ID) {
            stop(sprintf("Feature %s does not exists for omic %s.", feature, omicName))
        }

        featureData <- simuData[feature, , drop=FALSE]

        timeProfile <- dplyr::filter(simuSettings, .data$ID == feature) %>%
                dplyr::select(dplyr::starts_with("Group")) %>%
                tidyr::unite("Profile") %>%
                dplyr::pull(.data$Profile) %>%
                unique()

        featureMeansSE <- calculateRepMeans(featureData)

        yRange <-  range(featureData)

        timeProfiles <- unique(gsub(".*(Time[0-9]+).*", "\\1", colnames(featureData)))

        outputDF <- rbind(data.frame(), featureMeansSE %>%
                              tibble::rownames_to_column("ID") %>%
                              tidyr::gather(key = "Time", value = "Counts", -c(.data$ID)) %>%
                              tidyr::extract(.data$Time, c("Group", "Point", "Measure"), "(.*)\\.(Time[0-9]+)\\.(.*)") %>%
                              tidyr::spread(.data$Measure, .data$Counts) %>%
                              dplyr::mutate(Point = factor(.data$Point, levels = unique(timeProfiles), ordered = TRUE)) %>%
                              dplyr::mutate(Omic = as.factor(omicName))) %>%
                              dplyr::mutate(Rep = "Mean") %>%
                              dplyr::mutate(Point = factor(.data$Point, levels = unique(timeProfiles), ordered = TRUE))

        if(drawReps) {
            repDF <- tibble::rownames_to_column(data.frame(featureData)) %>%
                tidyr::gather(key = .data$Time, value = .data$Mean, -c(.data$rowname)) %>%
                tidyr::extract(.data$Time, c("Group", "Point", "Rep"), "(.*)\\.(Time[0-9]+)\\.(.*)") %>%
                dplyr::mutate(SD = 0, SE = 0, Profile = .data$timeProfile, Point = factor(.data$Point, levels = unique(.data$timeProfiles), ordered = TRUE), Omic = as.factor(.data$omicName)) %>%
                dplyr::rename(ID = .data$rowname)

            outputDF <- outputDF %>% dplyr::bind_rows(repDF)
        }

        return(outputDF)
    }


    if (length(omics) > 1) {
        # Currently limited to 2 omics
        omicNames <- utils::head(intersect(omics, names(featureIDS)), 2)

        omicsDF <- sapply(omicNames, function(omicName) {
            # Only 1 feature per omic
            omicFeature <- utils::head(featureIDS[[omicName]], 1)

            omicDF <- getOmicFeaturesDF(omicName, omicFeature)

            return(omicDF)
        }, simplify = FALSE, USE.NAMES = TRUE)

        featureIDS <- featureIDS[omicNames]

        # Reorder keeping the bigger in the first position
        omicsDF <- omicsDF[order(unlist(lapply(omicsDF, function(df) max(df$Mean))), decreasing = TRUE)]
        primaryDF <- omicsDF[[1]]

        # Scale range based on the total mean plus/minus SE
        scaleValues <- range(c(primaryDF$Mean + primaryDF$SE,
                               primaryDF$Mean - primaryDF$SE))

        omicsDF <- sapply(omicsDF, function(omicDF) {

            meanSE.scaled <- scales::rescale(c(omicDF$Mean,
                        omicDF$Mean + omicDF$SE,
                        omicDF$Mean - omicDF$SE),
                        scaleValues)

            outDF <- omicDF %>%
                dplyr::mutate(ScaledMean = utils::head(meanSE.scaled, nrow(omicDF)),
                              ScaledSE = (.data$ScaledMean*.data$SE)/.data$Mean)

            return(outDF)
        }, simplify = FALSE, USE.NAMES = TRUE)

    } else {
        featureIDS <- setNames(utils::head(unlist(featureIDS), 1), omics)
        omicNames <- omics
        omicsDF <- list(getOmicFeaturesDF(omics, featureIDS))
    }

    primaryDF <- omicsDF[[1]]
    primaryOmic <- omicNames[[1]]

    if (drawReps) {
        outputGgplot <- ggplot2::ggplot(data=primaryDF,
                                        ggplot2::aes(x=.data$Point, y=.data$Mean, group=.data$Rep, color=.data$Rep))
    } else {
        outputGgplot <- ggplot2::ggplot(data=primaryDF,
                                        ggplot2::aes(x=.data$Point, y=.data$Mean, group=.data$Omic, color=.data$Omic))
    }

    outputGgplot <- outputGgplot +
        ggplot2::geom_errorbar(ggplot2::aes(ymin=.data$Mean-.data$SE, ymax=.data$Mean+.data$SE), width=.1) +
        ggplot2::geom_line() +  ggplot2::geom_point()

    plotTitle <- featureIDS[[primaryOmic]]

    if (length(omics) > 1) {
        secondaryDF <- omicsDF[[2]]
        secondaryOmic <- omicNames[[2]]

        absoluteRange <- range(c(secondaryDF$Mean + secondaryDF$SE,
                                 secondaryDF$Mean - secondaryDF$SE))

        outputGgplot <- outputGgplot +
            ggplot2::geom_line(data = secondaryDF, ggplot2::aes(x = .data$Point, y = .data$ScaledMean, group = .data$Omic)) +
            ggplot2::scale_y_continuous(limits = scaleValues, sec.axis = ggplot2::sec_axis(~scales::rescale(., absoluteRange), name = featureIDS[[secondaryOmic]])) +
            ggplot2::geom_errorbar(data = secondaryDF, ggplot2::aes(ymin=.data$ScaledMean-.data$ScaledSE, ymax=.data$ScaledMean+.data$ScaledSE), width=.1) +
            ggplot2::geom_point(data = secondaryDF, ggplot2::aes(y = .data$ScaledMean)) +
            ggplot2::labs(color = "Omic")

        plotTitle <- paste0(plotTitle, ' - ', featureIDS[[secondaryOmic]])
    }

    outputGgplot <- outputGgplot +
        ggplot2::ggtitle(plotTitle) +
        ggplot2::ylab(featureIDS[[primaryOmic]]) +
        ggplot2::xlab("Time point") +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
        ggplot2::facet_grid(. ~ Group)

    return(outputGgplot)
}

#' Default data
#'
#' Dataset with base counts and id-gene tables.
#'
#' @details List with 6 elements:
#' \describe{
#'   \item{SimRNAseq}{\describe{
#'   \item{data}{Dataframe with base counts with gene id as rownames.}
#'   \item{geneLength}{Length of every gene.}}}
#'   \item{SimChIPseq}{\describe{
#'   \item{data}{Dataframe with base counts with regions as rownames.}
#'   \item{idToGene}{Dataframe with region as "ID" column and gene name on "Gene" column.}}}
#'   \item{SimDNaseseq}{\describe{
#'   \item{data}{Dataframe with base counts with regions as rownames.}
#'   \item{idToGene}{Dataframe with region as "ID" column and gene name on "Gene" column.}}}
#'   \item{SimMiRNAseq}{\describe{
#'   \item{data}{Dataframe with base counts with miRNA id as rownames.}
#'   \item{idToGene}{Dataframe with miRNA as "ID" column and gene name on "Gene" column.}}}
#'   \item{SimMethylseq}{\describe{
#'   \item{idToGene}{Dataframe with region as "ID" column and gene name on "Gene" column.}}}
#'   \item{CpGisland}{Dataframe of CpG to be used as initialization data, located on "Region" column}
#' }
"sampleData"
