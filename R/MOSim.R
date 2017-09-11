#' @include Simulation.R
#' @importFrom dplyr %>%
NULL

#' MOSim
#'
#' Description
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
#' @param omics
#' @param omicsOptions
#' @param totalGenes A number with the total number of genes including not expressed. Overwrited
#'   if a genome reference is provided. Currently not used as we force to provide real data.
#' @param diffGenes A number with the total number of differential genes (if value > 1) or % or
#'  total genes (if value < 1).
#' @param numberReps Number of replicates of the experiment.
#' @param numberGroups Number of samples considered on the experiment.
#' @param randomSeed Random seed for random number generator state.
#' @param times Numeric vector containing the measured times. If numberGroups < 2,
#'  the number of times must be at least 2.
#' @param exprGenes Total number of expressed genes (if value > 1) or %
#'  of the total genes (if value < 1).
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
#'  # Start empty simulation with default params:
#'  moSimulation <- mosim()
#'
#'  # Retrieve simulated data, it will only contain a "SimRNAseq" key.
#'  moSimulation@simData
#'
#'
#'  # Simulation with every omic and some custom options:
#'  simulatorOptions <- list(
#'      'SimRNAseq'=list(
#'          'depth'=25
#'      ),
#'      'SimMethylseq'=list(
#'      ),
#'      'SimmiRNAseq'=list(),
#'      'SimDNaseseq'=list(),
#'      'SimChIPseq'=list()
#'  )
#'
#'  moSimulation <- mosim(
#'      numberReps = 3,
#'      times = c(0, 2, 6, 12, 24),
#'      randomSeed = 1234,
#'      omics = simulatorOptions
#'  )
#'
#'  # simData will contain a key for every simulator
#'  dataRNAseq <- moSimulation@simData$SimRNAseq
#'
#'  # Methylseq is one exception, having both "beta" and "M" keys
#'  dataMethylSeqBeta <- moSimulation@simData$SimMethylseq$beta
#' }
mosim <- function(omics, omicsOptions = NULL, ...) {
    # Params to initialize simulation instance
    simParams <- list(...)

    # Select unique names
    omics <- unique(omics)

    # 'omics' parameter alias of 'simulators'
    simParams$simulators <- omics

    # Start logging
    # if (simParams$debug) {
    #     # logging::basicConfig(level = "FINEST")
    #     # logging::addHandler(logging::writeToFile, file = "mosim_debug.log", level = "DEBUG")
    # }

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

#' omicData
#'
#' @param omic
#' @param data
#' @param associationList
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' }
omicData <- function(omic, data, associationList = NULL) {
    # Convert the name to the proper class name
    omicClass <- paste0("Sim", gsub("-", "", omic))

    # Create the new instance
    omicSim <- new(omicClass,
                   "data" = data,
                   "idToGene" = associationList)

    return(omicSim)
}

#' omicSim
#'
#' @param omic
#' @param depth
#' @param totalFeatures
#' @param reguEffect
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' }
omicSim <- function(omic, depth = NULL, totalFeatures = NULL, reguEffect = NULL) {

    paramList <- list(
        "totalFeatures" = totalFeatures,
        "regulatorEffect" = reguEffect,
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

        regClasses <- omicsClasses[- grep("SimRNAseq", omicsClasses)]

        associationLists <- setNames(lapply(simulation@simulators[regClasses],
                                   slot, name = "idToGene"), regClasses)

        # Include only those regulators having an effect on genes
        if (only.linked) {
            associationLists <- sapply(names(associationLists), function(x) {

                omic.settings <- selectedSettings[[x]] %>% dplyr::filter(! is.na(Effect)) # TODO: change this to include Effect.GroupX

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

            output.df <- dplyr::select(x.data, Gene, ID, Effect, dplyr::starts_with("Effect.Group")) %>%
                dplyr::mutate(Omic = x)

            if (only.linked) {
                linked.IDs <- dplyr::filter(output.df, ! is.na(Effect))$ID # TODO: change this to include Effect.GroupX

                output.df <- dplyr::filter(output.df, ID %in% linked.IDs)
            }

            return(output.df)
        }))

        # Restore correct names
        associationLists <- setNames(associationLists,
                                     names(omicsClasses)[match(names(associationLists), omicsClasses)])

        outputList <- list(
            "association" = associationLists,
            "settings" = outputList,
            "regulators" = globalTable
        )
    }

    return(outputList)
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

    return(outputList)
}

#' Title
#'
#' @param simulation
#'
#' @return
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
