#' @import methods stringi matrixStats
#' @importFrom stats rnorm runif rbeta rbinom rpois quantile
#' @importFrom utils data
NULL

#'
#' This class manages the global simulation process, like associating genes with
#' gene classes, regulatory programs and other settings. Finally it will
#' initialize the simulators with their options that will use the previously
#' generated settings to simulate the data.
#'
#' @slot simulators Vector containing either S4 initialized classes of
#'   simulators or a list with the class name as keys, and its options as value,
#'   see example.
#' @slot totalGenes A number with the total number of genes including not
#'   expressed. Overwritten if a genome reference is provided. Currently not used
#'   as we force to provide real data.
#' @slot diffGenes A number with the total number of differential genes (if
#'   value > 1) or \% or total genes (if value < 1).
#' @slot numberReps Number of replicates of the experiment.
#' @slot numberGroups Number of samples considered on the experiment.
#' @slot times Numeric vector containing the measured times. If numberGroups <
#'   2, the number of times must be at least 2.
#' @slot geneNames Read only. List containing the IDs of the genes. Overwrited
#'   by the genome reference if provided. Currently not used as we force to
#'   provide real data.
#' @slot simSettings List of settings that overrides initializing the
#'   configuration of the simulation by passing a previously generated list.
#'   This could be used to tweak by hand the assigned profiles, genes,
#'   regulatory programs, etc.
#' @slot noiseFunction Noise function to apply when simulating counts. Must
#'   accept the parameter 'n' and return a vector of the same length. Defaults
#'   to `rnorm`
#' @slot profiles Named list containing the patterns with their coefficients.
#' @slot profileProbs Numeric vector with the probabilities to assign each of
#'   the patterns. Defaults to 0.2 for each.
#' @slot noiseParams Default noise parameters to be used with noise function.
#' @slot depth Default depth to simulate.
#' @slot TFtoGene Boolean (for default data) or 3 column data frame containing
#'   Symbol-TFGene-LinkedGene
#' @slot minMaxQuantile Numeric vector of length 2 indicating the quantiles to
#'   use in order to retrieve the absolute minimum and maximum value that a
#'   differentially expressed feature can have.
#' @slot minMaxFC Numeric vector of length 2 indicating the minimum and maximum
#'   fold-change that a differentially expressed feature can have.
#'
#' @export
#' @keywords internal
#'
#'
setClass(
    "MOSimulation",
    slots = c(
        simulators = "list",
        totalGenes = "numeric",
        diffGenes = "numeric",
        numberReps = "numeric",
        numberGroups = "numeric",
        times = "vector",
        geneNames = "character",
        simSettings = "list",
        profiles = "list",
        profileProbs = "list",
        noiseFunction = "function",
        noiseParams = "list",
        depth = "numeric",
        minMaxQuantile = "numeric",
        minMaxFC = "numeric",
        TFtoGene = "ANY"
    ),
    prototype = list(
        diffGenes = .15,
        numberReps = 3,
        numberGroups = 2,
        times = c(0, 2, 4, 12, 24),
        depth = 74,
        noiseFunction = stats::rnorm,
        noiseParams = list("sd" = 0.3),
        profiles = list(
            continuous.induction = c("a1", "b1", 0),
            continuous.repression = c("a1", "b1.neg", 0),
            transitory.induction = c("a2", "b2", "c2"),
            transitory.repression = c("a2.neg", "b2.neg", "c2.neg"),
            flat = c(0, 0, 0)
        ),
        profileProbs = list(
            continuous.induction = .235,
            continuous.repression = .235,
            transitory.induction = .235,
            transitory.repression = .235,
            flat = .06
        ),
        minMaxQuantile = c(0.3, 0.95),
        minMaxFC = c(3, 8),
        TFtoGene = NULL
    )
)

#'
#' Virtual class containing common methods and slots for child classes.
#'
#' @slot name Name of the simulator to be used in messages.
#' @slot data Data frame containing the initial sample to be used, with the
#'   features IDs as rownames and only one column named "Counts".
#' @slot regulator Boolean flag to indicate if the omic is a regulator or not.
#' @slot regulatorEffect Possible regulation effects of the omic (enhancer,
#'   repressor or both).
#' @slot idToGene Data frame with the association table between genes and other
#'   features. The structure must be 2 columns, one named "ID" and the other
#'   "Gene".
#' @slot min Minimum value allowed in the omic.
#' @slot max Maximum value allowed in the omic.
#' @slot depth Sequencing depth to simulate.
#' @slot depthRound Number of decimal places to round when adjusting depth.
#' @slot depthAdjust Boolean indicating whether to adjust by sequencing depth or
#'   not.
#' @slot totalFeatures Number of features to simulate. This will replace the
#'   data with a subset.
#' @slot noiseFunction Noise function to apply when simulating counts. Must
#'   accept the parameter 'n' and return a vector of the same length. Defaults
#'   to `rnorm`
#' @slot increment Read-only. Minimum value to increase when simulating counts.
#' @slot simData Contains the final simulated data.
#' @slot pregenerated Indicates if the child class will generate the simulated
#'   data instead of the general process.
#' @slot randData Auxiliary vector containing the original count data in random
#'   order with other adjustments.
#' @slot noiseParams Noise parameters to be used with noise function.
#' @slot roundDigits Number of digits to round the simulated count values.
#' @slot minMaxQuantile Numeric vector of length 2 indicating the quantiles to
#'   use in order to retrieve the absolute minimum and maximum value that a
#'   differentially expressed feature can have.
#' @slot minMaxFC Numeric vector of length 2 indicating the minimum and maximum
#'   fold-change that a differentially expressed feature can have.
#' @slot minMaxDist Named list containing different minimum and maximum
#'   constraints values calculated at the beginning of the simulation process.
#' @slot replicateParams Named list containing the parameters a and b to be used
#'   in the replicates generation process, see the vignette for more info.
#'
#' @export
#' @keywords internal
#'
setClass(
    "MOSimulator",
    slots = c(
        name = "character",
        data = "ANY",
        simData = "ANY",
        randData = "ANY",
        regulator = "logical",
        regulatorEffect = "vector",
        idToGene = "ANY",
        min = "numeric",
        max = "numeric",
        increment = "numeric",
        depth = "numeric",
        depthRound = "numeric",
        depthAdjust = "logical",
        noiseFunction = "function",
        noiseParams = "list",
        roundDigits = "numeric",
        pregenerated = "logical",
        totalFeatures = "numeric",
        minMaxQuantile = "numeric",
        minMaxFC = "numeric",
        minMaxDist = "list",
        replicateParams = "list"
    ),
    prototype = list(
        regulator = TRUE,
        min = 0,
        increment = 0,
        depthRound = 0,
        depthAdjust = TRUE,
        regulatorEffect = NA,
        roundDigits = 0,
        pregenerated = FALSE,
        totalFeatures = Inf,
        replicateParams = list(
            "a" = 0.01,
            "b" = 1.5
        )
    ),
    contains = c("VIRTUAL")
)

#'
#' Virtual class containing general methods for simulators based on regions of
#' the chromosomes, like DNase-seq, ChIP-seq or Methyl-seq
#'
#' @slot locs Vector containing the list of locations of the sites.
#' @slot locsName Type of the site to simulate, only for debug.
#' @slot splitChar Character symbol used to split identifiers in chr/start/end
#'
#' @export
#' @keywords internal
#' @rdname Simulator-class
#'
setClass(
    "MOSimulatorRegion",
    slots = c(
        locs = "ANY",
        splitChar = "character"
    ),
    prototype = list(
        splitChar = "_"
    ),
    contains = c("MOSimulator", "VIRTUAL")
)



#'
#' Class to simulate RNA-seq data
#'
#' @export
#' @keywords internal
#' @rdname Simulator-class
#'
setClass(
    "SimRNAseq",
    prototype = list(
        name = "RNA-seq",
        # idToGene = matrix(NA),
        regulator = FALSE,
        replicateParams = list(
            "a" = 0.5, #0.01,
            "b" = 1.5 #1.5
        )
    ),
    contains = "MOSimulator"
)

#'
#' Class to simulate transcription factor data
#'
#' @export
#' @keywords internal
#' @rdname Simulator-class
#'
setClass(
    "SimTF",
    prototype = list(
        name = "TF",
        # idToGene = matrix(NA),
        regulator = TRUE
    ),
    contains = "SimRNAseq"
)

#'
#' Class to simulate miRNA-seq
#'
#' @export
#' @keywords internal
#' @rdname Simulator-class
#'
setClass(
    "SimmiRNAseq",
    prototype = list(
        name = "miRNA-seq",
        regulator = TRUE,
        idToGene = matrix(),
        regulatorEffect = list(
            'repressor' = 0.05,
            'NE' = 0.95
        ),
        replicateParams = list(
            "a" = -0.065,
            "b" = 1.566
        )
    ),
    contains = "MOSimulator"
)

#'
#' Class to simulate ChIP-seq data
#'
#' @export
#' @keywords internal
#' @rdname Simulator-class
#'
setClass(
    "SimChIPseq",
    prototype = list(
        name = "ChIP-seq",
        regulator = TRUE,
        idToGene = matrix(),
        regulatorEffect = list(
            'enhancer' = 0.48,
            'repressor' = 0.48,
            'NE' = 0.02
        ),
        replicateParams = list(
            # "a" = -0.042,
            # "b" = 1.265,
            "a" = 0.5,
            "b" = 1.5
        )
    ),
    contains = "MOSimulatorRegion"
)

#'
#' Class to simulate DNase-seq data
#'
#' @export
#' @keywords internal
#' @rdname Simulator-class
#'
setClass(
    "SimDNaseseq",
    prototype = list(
        name = "DNase-seq",
        regulator = TRUE,
        idToGene = matrix(),
        regulatorEffect = list(
            'enhancer' = 0.48,
            'repressor' = 0.48,
            'NE' = 0.02
        ),
        replicateParams = list(
            # "a" = -0.042,
            # "b" = 1.265
            "a" = 0.5,
            "b" = 1.5
        )
    ),
    contains = "MOSimulatorRegion"
)

#'
#' Class to simulate Methyl-seq data.
#'
#' @slot nCpG numeric. Number of CpG sites to simulate.
#' @slot pSuccessMethReg numeric. Probability of success in methylated region.
#' @slot pSuccessDemethReg numeric. Probability of success in non methylated
#'   region
#' @slot errorMethReg numeric. Error rate in methylated region
#' @slot errorDemethReg numeric. Error rate in methylated region
#' @slot nReadsMethReg numeric. Mean number of reads in methylated region.
#' @slot nReadsDemethReg numeric. Mean number of reads in non methylated
#'   regions.
#' @slot phaseDiff numeric. Phase difference in the differentially methylated
#'   regions between two samples
#' @slot balanceHypoHyper numeric. Balance of hypo/hyper methylation
#' @slot ratesHMMMatrix numeric. Matrix of values that describes the exponential
#'   decay functions that define the distances between CpG values.
#' @slot distType character. Distribution used to generate replicates:
#' @slot transitionSize numeric.
#' @slot PhiMeth matrix. Transition matrix for CpG locations.
#' @slot PhiDemeth matrix. <Not used>
#' @slot typesLocation numeric. <Not used>
#' @slot returnValue character. Selected column:
#' @slot betaThreshold numeric. Beta threshold value used to calculate M values.
#'
#' @export
#' @keywords internal
#' @rdname Simulator-class
#'
setClass(
    "SimMethylseq",
    slots = c(
        nCpG = "numeric",
        betaThreshold = "numeric",
        WGBSparams = "list",
        Mvalues = "logical"
    ),
    prototype = list(
        name = "Methyl-seq",
        regulator = TRUE,
        Mvalues = TRUE,
        nCpG = 5000,
        betaThreshold = 0.01,
        # maxValue = 1,
        # nRepeats = 1,
        idToGene = matrix(),
        pregenerated = TRUE,
        roundDigits = 5,
        WGBSparams = list(
            # Phase difference in the differentially methylated regions between two samples
            phase_diff = c(0, 0.1),
            #error_rate_in_differentially_methylated_region = 0.1
            # File to write results
            outfile = '/tmp/',
            #output_path
            probs = c(1, 1, 0.9, 0.8, 0.7, 0.6),
            # Transition size
            transition_size = 0,
            # Number of repeats (NOT REPLICATES!)
            m = 1,
            # Transition matrix for CpG locations
            Pi_m = matrix(c(0.65, 0.35, 0.2, 0.8), byrow = 
                            TRUE, nrow = 2),
            # Transition matrix for probability distributions
            Pi_d = matrix(c(0.9, 0.1, 0.1, 0.9), byrow = 
                            TRUE, nrow = 2),
            # Type of locations
            type_of_locations = 2,
            # Mean number of reads in methylated region
            mean_m = 29,
            # Mean number of reads in non methylated regions
            mean_d = 29,
            # Probability of success in  methylated region
            prob_m = 0.9203,
            # Probability of success in non methylated region
            prob_d = 0.076,
            # Error rate in dif. methylated region
            error_m = 0.1,
            # Error rate in non methylated region
            error_d = 0.1,
            # Balance: this is the balance of hypo/hyper methylation
            balance = 0.5,
            # Rates for HMM for CpG locations: this is the matrix of values that
            # describes the exponential decay functions that define the distances
            # between CpG values.
            rates_for_HMM_for_CpG_locations = c(0.019, 0.002),
            # Type: this can be either binomial or truncated and defines the
            # distribution of methylated reads at each CpG.
            type = "binomial",
            # methReads = 1, demethReads = 2, totalReads = 3, proportionReads = 4
            returnColumn = 4,
            distType = "nbin"#binomial"
            ),
        regulatorEffect = list(
            'enhancer' = 0.48,
            'repressor' = 0.48,
            'NE' = 0.01
        ),
        depthAdjust = FALSE
        ),
    contains = "MOSimulatorRegion"
)
