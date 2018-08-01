#' @include Simulation.R
#' @importFrom dplyr %>%
#' @importFrom rlang := !!
NULL

#' MOSim
#'
#' Multiomics simulation package.
#'
#'
#' @docType package
#' @name MOSim
NULL
#> NULL

#' mosim
#'
#' Performs a multiomic simulation by chaining two actions:
#'  1) Creating the "Simulation" class with the provided params.
#'  2) Calling "simulate" method on the initialized object.
#'
#' @param omics Character vector containing the list of omics to simulate it can also be an
#'  associative list with the omic names as keys and their options
#'  as values. This is the prefered way of indicating the omics to simulate using the default
#'  parameters.
#' @param omicsOptions Associative list containing the options for each omic to simulate. This is
#' used in conjunction with the helper methods omicSim to create the associative list
#' in a friendly way, and omicData to provide custom data; see the related section for
#' more information.
#' @param diffGenes A number with the total number of differential genes (if value > 1) or % or
#'  total genes (if value < 1).
#' @param numberReps Number of replicates of the experiment.
#' @param numberGroups Number of samples considered on the experiment.
#' @param randomSeed Random seed for random number generator state.
#' @param times Numeric vector containing the measured times. If numberGroups < 2,
#'  the number of times must be at least 2.
#' @param TFtoGene Boolean indicating if default transcription factors data should be used or not, or
#' a 3 column data frame containing custom associations.
#' @param noiseFunction Noise function to apply when simulating counts. Must accept the parameter 'n' and
#'  return a vector of the same length. Defaults to `rnorm`
#' @param noiseParam Named list with the parameters to apply to the noise function.
#' @param profiles Named list containing the patterns with their coefficients.
#' @param profileProbs Numeric vector with the probabilities to assign each of the patterns. Defaults to 0.2 for each.
#' @param ... Additional options to be passed to simulation object
#'
#' @return Instance of class "Simulation" with slot "simData" containing the multiomic simulation data.
#' @export
#'
#' @examples
#' \dontrun{
#'
#'  moSimulation <- mosim(
#'      omics = c("RNA-seq"),
#'      numberReps = 3,
#'      times = c(0, 2, 6, 12, 24),
#'      randomSeed = 1234
#'  )
#'
#'  # Retrieve simulated count matrix for RNA-seq
#'  dataRNAseq <- omicResults(moSimulation, "RNA-seq")
#'
#' }
mosim <- function(omics, omicsOptions = NULL, ...) {
    # Params to initialize simulation instance
    simParams <- list(...)

    # Select unique names
    omics <- if (is.list(omics)) omics[unique(names(omics))] else unique(omics)

    # 'omics' parameter alias of 'simulators'
    simParams$simulators <- omics

    # If it is a plain vector, transform it to a list of empty lists
    if (! is.list(simParams$simulators)) {
        simParams$simulators <- setNames(rep(list(list()), length(simParams$simulators)),
                                         simParams$simulators)
    }

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

    # Remove keys not present in 'Simulation' class slots
    simParams <- simParams[names(simParams) %in% slotNames("Simulation")]

    oSim <- do.call(new, c("Class"="Simulation", simParams))
    oSim <- simulate(oSim)

    return(oSim)
}

#' Set customized data for an omic.
#'
#' @param omic The name of the omic to provide data.
#' @param data Data frame with the omic identifiers as row names and just one column named Counts
#' containing numeric values used as initial sample for the simulation.
#' @param associationList Only for regulatory omics, a data frame with 2 columns, the first called
#' containing the regulator ID and the second called Gene with the gene identifier.
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' }
omicData <- function(omic, data = NULL, associationList = NULL) {
    # Convert the name to the proper class name
    omicClass <- paste0("Sim", gsub("-", "", omic))

    # Create the new instance
    omicSim <- new(omicClass,
                   "data" = data,
                   "idToGene" = associationList)

    return(omicSim)
}

#' Set the simulation settings for an omic.
#'
#' @param omic Name of the omic to set the settings.
#' @param depth Sequencing depth in millions of counts. If not provided will take the global
#' parameter passed to mosim function.
#' @param totalFeatures Limit the number of features to simulate. By default include all present
#' in the dataset.
#' @param regulatorEffect only for regulatory omics. Associative list containing the percentage
#' of effects over the total number of regulator, including repressor, association and
#' no effect (NE).
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' }
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
#' @param association A boolean indicating if the association must also be returned for the regulators.
#' @param reverse A boolean, swap the column order in the association list in case we want to use the
#' output directly and the program requires a different ordering.
#' @param only.linked Return only the interactions that have an effect.
#'
#' @return A list containing a data frame with the settings used to simulate each of the indicated omics.
#' If association is TRUE, it will be a list with 2 keys: 'associations' and 'settings', with each one being
#' a list containing the information for the selected omics.
#' @export
#'
#' @examples
omicSettings <- function(simulation, omics = NULL, association = FALSE, reverse = FALSE, only.linked = FALSE, prefix = FALSE) {
    # TODO: should this be a generic <Simulation> method?

    # Select all omics by default
    if (is.null(omics)) {
        omics <- lapply(simulation@simulators, slot, name = "name")
    } else {
        names(omics) <- paste0("Sim", gsub("-", "", omics))
    }

    # Convert the name to the proper class name
    omicsClasses <- setNames(paste0("Sim", gsub("-", "", omics)), omics)

    selectedSettings <- simulation@simSettings$geneProfiles[omicsClasses]

    # Remove Effect from settings and rename Effect.Linked to Effect
    # selectedSettings <- lapply(selectedSettings, function(x) {
    #     if (is.element("Effect", colnames(x))) {
    #         x <- dplyr::select(x, -Effect) %>%
    #                          dplyr::rename(Effect = Effect.Linked)
    #     }
    #
    #     return(x)
    # })

    outputList <- setNames(selectedSettings, omics)

    # Add the associations
    if (association) {

        regClasses <- omicsClasses[grep("SimRNAseq", omicsClasses, invert = TRUE)]

        associationLists <- setNames(lapply(simulation@simulators[regClasses],
                                   slot, name = "idToGene"), regClasses)

        # Include only those regulators having an effect on genes
        if (only.linked) {
            associationLists <- sapply(names(associationLists), function(x) {
                omic.settings <- selectedSettings[[x]] %>% dplyr::filter_at(dplyr::vars(dplyr::starts_with("Effect")), dplyr::any_vars(! is.na(.)))

                omic.association <- associationLists[[x]] %>% dplyr::filter(ID %in% omic.settings$ID)

                return(omic.association)
            }, simplify = FALSE)
        }

        # Prepend the name of the omic to the regulator ID
        if (prefix) {
            associationLists <- sapply(names(associationLists), function(x) {

                prefix.settings <- associationLists[[x]] %>% dplyr::mutate(ID = paste0(x, ID))

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
            x.settings <- dplyr::filter(x.settings, !is.na(Effect)) %>%
                dplyr::select(ID, dplyr::starts_with("Group")) %>%
                dplyr::distinct()

            output.df <- dplyr::select(x.data, Gene, ID, Effect, dplyr::starts_with("Effect.Group"), dplyr::starts_with("Tmax.Group")) %>%
                dplyr::mutate(Omic = x) %>% dplyr::left_join(x.settings, by = c("ID" = "ID")) %>%
                dplyr::mutate_at(dplyr::vars(dplyr::starts_with("Group")), dplyr::funs(ifelse(is.na(.), "flat", .)))

            if (only.linked) {
                linked.IDs <- dplyr::filter_at(output.df, dplyr::vars(dplyr::starts_with("Effect")), dplyr::any_vars(! is.na(.)))$ID

                output.df <- dplyr::filter(output.df, ID %in% linked.IDs)
            }

            return(output.df)
        }))

        # Remove & rename effects from data frames
        replace_effect <- function(df) {
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

                # Assign the same effect to regulator rows.
                df <- dplyr::group_by(df, ID) %>%
                    dplyr::mutate_at(dplyr::vars(dplyr::starts_with("Group")),
                                     dplyr::funs(propagate_profile(., Effect)))
            }

            return(df[, colnames(df) != "Effect"])
        }

        # Restore correct names
        associationLists <- setNames(associationLists,
                                     names(omicsClasses)[match(names(associationLists), omicsClasses)])

        # Replace effects on settings
        outputList <- sapply(outputList, replace_effect, simplify = FALSE, USE.NAMES = TRUE)

        outputList <- list(
            "association" = associationLists,
            "settings" = outputList,
            "regulators" = replace_effect(globalTable)
        )
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
#' @return A list containing an element for every omic specifiec, with
#' the simulation data in the format indicated.
#' @export
#'
#' @examples
omicResults <- function(simulation, omics = NULL, format = data.frame) {
    # Select all omics by default
    if (is.null(omics)) {
        omics <- lapply(simulation@simulators, slot, name = "name")
    }

    # Convert the name to the proper class name
    omicsClasses <- paste0("Sim", gsub("-", "", omics))

    selectedSimulators <- lapply(simulation@simulators[omicsClasses], slot, name = "simData")

    outputList <- setNames(selectedSimulators, omics)

    if (length(omics) > 1) {
        return(outputList)
    } else {
        return(outputList[[omics]])
    }
}

#' Retrieves the experimental design
#'
#' @param simulation
#'
#' @return A Simulation object.
#' @export
#'
#' @examples
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

#' Default data
#'
#' Dataset with base counts and id-gene tables.
#'
#' @format List with 6 elements:
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
