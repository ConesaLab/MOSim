#' @include AllClass.R
NULL

################################################################################
# General
################################################################################
#' Simulate [internal use]
#'
#' @rdname Generics
#' @section Shared methods:
#'
#' This generic can dispatch methods for two types of classes: \describe{
#' \item{\linkS4class{MOSimulation}}{in which case triggers the simulation on
#' every simulator loaded.}
#' \item{\linkS4class{MOSimulator}}{performs the simulation based on already
#' initialized simulation settings.} }
#'
#' @param object Object of class \linkS4class{MOSimulator} or \linkS4class{MOSimulation}.
#' @param simulation Only in \linkS4class{MOSimulator} class. Initialized instance of Simulation class
#' @param ... Extra parameters for extensibility.
#'
#' @return A \linkS4class{MOSimulator} object containing the simulation data inside @simData
#'  if called on a \linkS4class{MOSimulator} class, or a \linkS4class{MOSimulation} object with all @simulators
#'  containing simulated data.
#' @keywords internal
#' @noRd
#'
setGeneric("simulate", function(object, ...) standardGeneric("simulate"))

#' Simulate [internal use]
#'
#' @rdname Generics
#' @section Shared methods:
#'
#' This generic can dispatch methods from two types of classes:
#' \describe{
#'      \item{\linkS4class{MOSimulation}}{show the general simulation settings (experimental design, default depth...)}
#'      \item{\linkS4class{MOSimulator}}{individual omic options (depth, regulator effect...)}
#' }
#'
#' @param object Object of class \linkS4class{MOSimulator} or \linkS4class{MOSimulation}.
#' @return The simulation settings used for the general process (for MOSimulation class) or the individual omic (MOSimulator class).
#' @keywords internal
#' @noRd
#'
setGeneric("simSettings", function(object) standardGeneric("simSettings"))

################################################################################
# Simulator
################################################################################



#' Initialize data
#'
#' @rdname Generics
#' @section Simulator methods:
#'
#' This method will clone the seed data once per group to simulate and
#' establish differences using the simulation settings.
#'
#' @param object Instance of MOSimulator class.
#' @param simulation Instance of MOSimulation class.
#'
#' @return An instance of MOSimulator class with modified options.
#' @noRd
#'
setGeneric("initializeData", function(object, simulation) standardGeneric("initializeData"))

#' postSimulation
#'
#' @rdname Generics
#' @section Simulator methods:
#'
#' Method called after performing the simulation, allowing different actions
#' like adjusting depth, rounding values and so on, on a simulator basis.
#'
#' @param object Instance of MOSimulator class.
#' @param simulation Instance of MOSimulation class.
#'
#' @return An instance of MOSimulator class with modified options.
#' @keywords internal
#' @noRd

#'
setGeneric("postSimulation", function(object, simulation) standardGeneric("postSimulation"))

#' IDfromGenes
#'
#' @rdname Generics
#' @section Simulator methods:
#'
#' Method for transforming genes to IDs (regions, miRNA...)
#'
#' @param object Instance of a MOSimulator class
#' @param geneNames Names of genes to look up in the association table.
#'
#' @return IDs corresponding to the genes.
#' @keywords internal
#' @noRd
#'
setGeneric("IDfromGenes", function(object, geneNames, simplify = TRUE) standardGeneric("IDfromGenes"))


#' IDtoGenes
#'
#' @rdname Generics
#' @section Simulator methods:
#'
#' Method for transforming IDs (regions, miRNAs, ...) to genes
#'
#' @param object Instance of a MOSimulator class.
#' @param idNames IDs to look up in the association table.
#'
#' @return Genes corresponding to the IDs.
#' @keywords internal
#' @noRd
#'
setGeneric("IDtoGenes", function(object, idNames, simplify = TRUE) standardGeneric("IDtoGenes"))

#' simulateParams
#'
#' @rdname Generics
#' @section Simulator methods:
#'
#' Method to allow individual simulator classes to change the default behaviour
#' of generating noise and time coefficients.
#'
#' @param object Instance of a MOSimulator class.
#' @param simulation  Instance of a MOSimulation class
#' @param counts Initial count values to simulate.
#' @param profiles Generated simulation settings for the group being simulated.
#' @param group Group being simulated.
#' @param ids ID of every row.
#'
#' @return List with the following elements:
#' \describe{
#'      \item{randomCounts}{numeric vector containing random count values.}
#'      \item{noiseValues}{numeric vector containing noise values generated with the noise function and parameters specified.}
#'      \item{m}{numeric vector of lower values comparing original and random counts.}
#'      \item{M}{numeric vector of maximum values comparing original and random counts.}
#' }
#' @keywords internal
#' @noRd
#'
setGeneric("simulateParams", function(object, simulation, counts, profiles, group, ids, ...) standardGeneric("simulateParams"))

#' adjustProfiles
#'
#' @rdname Generics
#' @section Simulator methods:
#'
#' Modify the configuration profile generated when initializing simulation settings
#' to allow for certain constrainsts if necessary.
#'
#' @param object Instance of MOSimulator class.
#' @param simulation Instance of MOSimulation class.
#' @param profiles Data frame containing the generated profile configuration for the
#'  simulator.
#' @param step Two possible options "Effect" and "Groups" indicating in which point
#' of the settings simulation is called.
#'
#' @return The data.frame profile with the necessary modifications depending on
#' each child class.
#' @keywords internal
#' @noRd
#'
setGeneric("adjustProfiles", function(object, simulation, profiles, step) standardGeneric("adjustProfiles"))


################################################################################
# SimulatorRegion
################################################################################
#' locGRanges
#'
#' @rdname Generics
#' @section Simulator methods:
#'
#' Creates a very basic GRanges object based on locations
#'
#' @param object Instance of a SimulatorRegion class
#' @param locs Positions of locations
#'
#' @return GRanges object with 1 as seqname and locs as starting positions.
#' @keywords internal
#' @noRd
#'
setGeneric("locGRanges", function(object, locs) standardGeneric("locGRanges"))

#' regionNames
#'
#' @rdname Generics
#' @section Simulator methods:
#'
#' Creates a valid string of a region
#'
#' @param object Instace of a SimulatorRegion class
#' @param chrNumber IDs of the chromosomes
#' @param start Starts positions of the regions
#' @param end End positions of the regions
#'
#' @return A character vector with the format "<chr>_<start>_<end>"
#' @keywords internal
#' @noRd
#'
setGeneric("regionNames", function(object, chrNumber, start, end=NULL) standardGeneric("regionNames"))
