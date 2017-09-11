#' @include AllClass.R
NULL

################################################################################
# General
################################################################################
#' Simulate [internal use]
#'
#' This generic can dispatch methods for two types of classes:
#' \describe{
#'      \item{\linkS4class{Simulation}}{in which case triggers the simulation on every simulator
#'  loaded.}
#'      \item{\linkS4class{Simulator}}{performs the simulation based on already initialized simulation
#'  settings.}
#' }
#'
#' @rdname initialize-methods
#' @param object Object of class \linkS4class{Simulator} or \linkS4class{Simulation}.
#' @param simulation Only in \linkS4class{Simulator} class. Initialized instance of Simulation class
#' @param ... Extra parameters for extensibility.
#'
#' @return A \linkS4class{Simulator} object containing the simulation data inside @simData
#'  if called on a \linkS4class{Simulator} class, or a \linkS4class{Simulation} object with all @simulators
#'  containing simulated data.
#'
setGeneric("simulate", function(object, ...) standardGeneric("simulate"))

#' Simulate [internal use]
#'
#' This generic can dispatch methods from two types of classes:
#' \describe{
#'      \item{\linkS4class{Simulation}}{show the general simulation settings (experimental design, default depth...)}
#'      \item{\linkS4class{Simulator}}{individual omic options (depth, regulator effect...)}
#' }
#'
#' @rdname simulate-methods
#' @param object Object of class \linkS4class{Simulator} or \linkS4class{Simulation}.
#'
setGeneric("simSettings", function(object) standardGeneric("simSettings"))

################################################################################
# Simulator
################################################################################

setGeneric("initializeData", function(object, simulation) standardGeneric("initializeData"))

#' postSimulation
#'
#' Method called after performing the simulation, allowing different actions
#' like adjusting depth, rounding values and so on, on a simulator basis.
#'
#' @param object Instance of Simulator class.
#' @param simulation Instance of Simulation class.
#'
#' @return An instance of Simulator class with modified options.
#'
setGeneric("postSimulation", function(object, simulation) standardGeneric("postSimulation"))

#' IDfromGenes
#'
#' Method for transforming genes to IDs (regions, miRNA...)
#'
#' @param object Instance of a Simulator class
#' @param geneNames Names of genes to look up in the association table.
#'
#' @return IDs corresponding to the genes.
#'
setGeneric("IDfromGenes", function(object, geneNames, simplify = TRUE) standardGeneric("IDfromGenes"))


#' IDtoGenes
#'
#' Method for transforming IDs (regions, miRNAs, ...) to genes
#'
#' @param object Instance of a simulator class.
#' @param idNames IDs to look up in the association table.
#'
#' @return Genes corresponding to the IDs.
#'
setGeneric("IDtoGenes", function(object, idNames, simplify = TRUE) standardGeneric("IDtoGenes"))

#' simulateParams
#'
#' Method to allow individual simulator classes to change the default behaviour
#' of generating noise and time coefficients.
#'
#' @param object Instance of a simulator class.
#' @param simulation  Instance of a simulation class
#' @param counts Initial count values to simulate.
#' @param profiles Generated simulation settings for the group being simulated.
#' @param group Group being simulated.
#' @param ids ID of every row.
#'
#' @return List with the following elements:
#' \describe{
#'  \item{}{}
#' }
#'
#'TODO: MODIFY THE GENERIC
setGeneric("simulateParams", function(object, simulation, counts, profiles, group, ids, ...) standardGeneric("simulateParams"))

#' adjustProfiles
#'
#' Modify the configuration profile generated when initializing simulation settings
#' to allow for certain constrainsts if necessary.
#'
#' @param object Instance of simulator class.
#' @param simulation Instance of simulation class.
#' @param profiles Data.frame containing the generated profile configuration for the
#'  simulator.
#' @param step Two possible options "Effect" and "Groups" indicating in which point
#' of the settings simulation is called.
#'
#' @return The data.frame profile with the necessary modifications depending on
#' each child class.
#'
setGeneric("adjustProfiles", function(object, simulation, profiles, step) standardGeneric("adjustProfiles"))


################################################################################
# SimulatorRegion
################################################################################
#' locGRanges
#'
#' Creates a very basic GRanges object based on locations
#'
#' @param object Instance of a SimulatorRegion class
#' @param locs Positions of locations
#'
#' @return GRanges object with 1 as seqname and locs as starting positions.
#'
setGeneric("locGRanges", function(object, locs) standardGeneric("locGRanges"))

#' regionNames
#'
#' Creates a valid string of a region
#'
#' @param object Instace of a SimulatorRegion class
#' @param chrNumber IDs of the chromosomes
#' @param start Starts positions of the regions
#' @param end End positions of the regions
#'
#' @return A character vector with the format "<chr>_<start>_<end>"
#'
setGeneric("regionNames", function(object, chrNumber, start, end=NULL) standardGeneric("regionNames"))


setGeneric("simDebug", function(object, simulation, ...) standardGeneric("simDebug"))
