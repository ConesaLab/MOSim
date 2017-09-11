#' @import methods stringi matrixStats
NULL

#'
#' This class manages the global simulation process, like associating genes
#' with gene classes, regulatory programs and other settings. Finally it will
#' initialize the simulators with their options that will use the previously
#' generated settings to simulate the data.
#'
#' @slot simulators Vector containing either S4 initialized classes of simulators
#'  or a list with the class name as keys, and its options as value, see example.
#' @slot totalGenes A number with the total number of genes including not expressed. Overwrited
#'   if a genome reference is provided. Currently not used as we force to provide real data.
#' @slot diffGenes A number with the total number of differential genes (if value > 1) or % or
#'  total genes (if value < 1).
#' @slot numberReps Number of replicates of the experiment.
#' @slot numberGroups Number of samples considered on the experiment.
#' @slot randomSeed Random seed for random number generator state.
#' @slot times Numeric vector containing the measured times. If numberGroups < 2,
#'  the number of times must be at least 2.
#' @slot geneNames Read only. List containing the IDs of the genes. Overwrited by the
#'  genome reference if provided. Currently not used as we force to provide real data.
#' @slot simSettings List of settings that overrides initializing the configuration of the simulation
#'  by passing a previously generated list. This could be used to tweak by hand the
#'  assigned profiles, genes, regulatory programs, etc.
#' @slot noiseFunction Noise function to apply when simulating counts. Must accept the parameter 'n' and
#'  return a vector of the same length. Defaults to `rnorm`
#' @slot noiseParam Named list with the parameters to apply to the noise function.
#' @slot profiles Named list containing the patterns with their coefficients.
#' @slot profileProbs Numeric vector with the probabilities to assign each of the patterns. Defaults to 0.2 for each.
#'
#' @export
#'
#'
setClass(
    "Simulation",
    slots = c(
        simulators = "list",
        totalGenes = "numeric",
        diffGenes = "numeric",
        numberReps = "numeric",
        numberGroups = "numeric",
        randomSeed = "numeric",
        times = "vector",
        geneNames = "character",
        # defaultPatterns = "list",
        simSettings = "list",
        # exprGenes = "numeric",
        defaultData = "list",
        profiles = "list",
        profileProbs = "list",
        noiseFunction = "function",
        noiseParams = "list",
        depth = "numeric",
        debug = "logical",
        minMaxQuantile = "numeric",
        replicateParams = "list",
        TFtoGene = "ANY"
    ),
    prototype = list(
        diffGenes = .15,
        numberReps = 3,
        numberGroups = 2,
        randomSeed = 12345,
        times = c(0, 12),
        # exprGenes = 0.5,
        depth = 74,
        noiseFunction = rnorm,
        noiseParams = list("sd" = 0.3),
        profiles = list(
            continuous.induction = c("a1", "b1", 0),
            continuous.repression = c("a1", "-b1", 0),
            transitory.induction = c("a2", "b2", "c2"),
            transitory.repression = c("a2", "-b2", "-c2"),
            flat = c(0, 0, 0)
        ),
        profileProbs = list(
            continuous.induction = .43,
            continuous.repression = .43,
            transitory.induction = .06,
            transitory.repression = .06,
            flat = .06
        ),
        minMaxQuantile = c(0.25, 0.75),
        replicateParams = list(
            "a" = 0.01,
            "b" = 1.5
        ),
        TFtoGene = NULL,
        debug = FALSE
    )
)

#'
#' Virtual class containing common methods and slots for child classes.
#'
#' @slot name Name of the simulator to be used in messages.
#' @slot data Data frame containing the initial sample to be used, with the features IDs as rownames
#'  and only one column named "Counts".
#' @slot regulator Boolean flag to indicate if the omic is a regulator or not.
#' @slot regulatorEffect Possible regulation effects of the omic (enhancer, repressor or both).
#' @slot idToGene Data frame with the association table between genes and other features. The structure
#'  must be 2 columns, one named "ID" and the other "Gene".
#' @slot min Minimum value allowed in the omic.
#' @slot max Maximum value allowed in the omic.
#' @slot gammaTable List with the gamma distribution parameters needed.
#' @slot depth Sequencing depth to simulate.
#' @slot noise Noise value for the NB mean.
#' @slot depthRound Number of decimal places to round when adjusting depth.
#' @slot depthAdjust Boolean indicating whether to adjust by sequencing depth or not.
#' @slot totalFeaturesN umber of features to simulate. This will replace the data with a subset.
#' @slot noiseFunction Noise function to apply when simulating counts. Must accept the parameter 'n' and
#'  return a vector of the same length. Defaults to `rnorm`
#' @slot noiseParam Named list with the parameters to apply to the noise function.
#' @slot increment Read-only. Minimum value to increase when simulating counts.
#' @slot simData Contains the final simulated data.
#' @slot pregenerated Indicates if the child class will generate the simulated data instead of the
#' general process.
#'
#' @export
#'
setClass(
    "Simulator",
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
        gammaTable = "ANY",
        depthRound = "numeric",
        depthAdjust = "logical",
        # coeffIncrease = "numeric",
        noiseFunction = "function",
        noiseParams = "list",
        debugInfo = "list",
        roundDigits = "numeric",
        pregenerated = "logical",
        totalFeatures = "numeric",
        minMaxQuantile = "numeric",
        minMaxDist = "list",
        replicateParams = "list"
    ),
    prototype = list(
        regulator = TRUE,
        min = 0,
        increment = 0,
        gammaTable = NULL,
        depthRound = 0,
        depthAdjust = TRUE,
        # coeffIncrease = 10,
        regulatorEffect = NA,
        roundDigits = 0,
        pregenerated = FALSE,
        totalFeatures = Inf,
        gammaTable = list(
            mean=c(1,5,10,15,25,50,75,100,500,5000,15000,1e+06),

            shape=c(5.12069695369583,2.31111819433427,2.79973960715549,2.85912994880808,
                    2.92544440726716,2.71844286130629,1.79038505963967,2.93199534788786,
                    1.58803070351186,1.53988847271216,1.86332612328476,0.580739807001836),

            scale=c(0.0881209571080616,0.804836647784225,1.64942276542406,2.87435701954641,
                    5.21636953504161,8.65961524646501,15.7834953443805,13.6005097645718,
                    53.8061030463463,301.806512019322,622.399272725056,16859.692892429)
        )
    ),
    contains = c("VIRTUAL")
)

#'
#' Virtual class containing general methods for simulators based on
#' regions of the chromosomes, like DNase-seq, ChIP-seq or Methyl-seq
#'
#' @slot locs Vector containing the list of locations of the sites.
#' @slot locsName Type of the site to simulate, only for debug.
#' @slot splitChar Character symbol used to split identifiers in chr/start/end
#'
#' @export
#'
setClass(
    "SimulatorRegion",
    slots = c(
        locs = "ANY",
        splitChar = "character"
    ),
    prototype = list(
        splitChar = "_"
    ),
    contains = c("Simulator", "VIRTUAL")
)



#'
#' Class to simulate RNA-seq data
#'
#' @export
#'
setClass(
    "SimRNAseq",
    prototype = list(
        name = "RNA-seq",
        # idToGene = matrix(NA),
        regulator = FALSE,
        replicateParams = list(
            "a" = 0.01,
            "b" = 1.5
        )
        # gammaTable = list(
        #     mean=c(1,5,10,15,25,50,75,100,500,5000,15000,1e+06),
        #
        #     shape=c(5.12069695369583,2.31111819433427,2.79973960715549,2.85912994880808,
        #             2.92544440726716,2.71844286130629,1.79038505963967,2.93199534788786,
        #             1.58803070351186,1.53988847271216,1.86332612328476,0.580739807001836),
        #
        #     scale=c(0.0881209571080616,0.804836647784225,1.64942276542406,2.87435701954641,
        #             5.21636953504161,8.65961524646501,15.7834953443805,13.6005097645718,
        #             53.8061030463463,301.806512019322,622.399272725056,16859.692892429)
        # )
    ),
    contains = "Simulator"
)

#'
#' Class to simulate RNA-seq data
#'
#' @export
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
#'
setClass(
    "SimmiRNAseq",
    prototype = list(
        name = "miRNA-seq",
        regulator = TRUE,
        idToGene = matrix(),
        regulatorEffect = c('repressor'),
        replicateParams = list(
            "a" = -0.065,
            "b" = 1.566
        )
        # gammaTable = list(
        #     mean = c(0.5, 1, 2, 30, 1e+06),
        #
        #     shape = c(
        #         0.0912798260161212,
        #         2.57146261146076,
        #         2.73546025404724,
        #         1.43199997857507,
        #         0.238274308304151
        #     ),
        #
        #     scale = c(
        #         0.170643135401093,
        #         0.116652504719918,
        #         0.207873707983415,
        #         1.74063989947033,
        #         5842.72780240922
        #     )
        # )
    ),
    contains = "Simulator"
)

#'
#' Class to simulate ChIP-seq data
#'
#' @export
#'
setClass(
    "SimChIPseq",
    #slots=c(bsSize = "numeric"),
    prototype = list(
        name = "ChIP-seq",
        regulator = TRUE,
        idToGene = matrix(),
        regulatorEffect = c('enhancer', 'repressor'),
        replicateParams = list(
            "a" = -0.042,
            "b" = 1.265
        )
        # gammaTable = list(
        #     mean = c(0.5, 1, 3, 5, 10, 20, 1e+06),
        #
        #     shape = c(
        #         0.0769368756183698,
        #         4.87463148740777,
        #         1.77614021569757,
        #         1.27135510976813,
        #         1.10859564566844,
        #         1.11429964760339,
        #         0.724838782801546
        #     ),
        #
        #     scale = c(
        #         0.668787864915357,
        #         0.192535468225583,
        #         1.05202204909108,
        #         2.15309258180485,
        #         2.66610983570308,
        #         3.11451405834024,
        #         10.2903610306696
        #     )
        # )
    ),
    contains = "SimulatorRegion"
)

#'
#' Class to simulate DNase-seq data
#'
#' @export
#'
setClass(
    "SimDNaseseq",
    prototype = list(
        name = "DNase-seq",
        regulator = TRUE,
        idToGene = matrix(),
        regulatorEffect = c('enhancer', 'repressor'),
        replicateParams = list(
            "a" = -0.042,
            "b" = 1.265
        )
        # TODO: revert this
        # gammaTable = list(
        #     mean=c(0.5, 1, 5, 10, 15, 25, 50, 75, 100, 500, 5000,15000,1e+06),
        #
        #     shape=c(0.7, 5.12069695369583,2.31111819433427,2.79973960715549,2.85912994880808,
        #             2.92544440726716,2.71844286130629,1.79038505963967,2.93199534788786,
        #             1.58803070351186,1.53988847271216,1.86332612328476,0.580739807001836),
        #
        #     scale=c(0.668787864915357, 0.0881209571080616,0.804836647784225,1.64942276542406,2.87435701954641,
        #             5.21636953504161,8.65961524646501,15.7834953443805,13.6005097645718,
        #             53.8061030463463,301.806512019322,622.399272725056,16859.692892429)
        # )
        # gammaTable = list(
        #     mean = c(5, 10, 20, 30, 40, 50, 100, 1e+06),
        #
        #     shape = c(
        #         2.65652390488234,
        #         2.61251029572217,
        #         2.32879608523843,
        #         2.26578855649232,
        #         2.24077370888784,
        #         2.22805845900357,
        #         2.27087586897883,
        #         0.955428851299402
        #     ),
        #
        #     scale = c(
        #         0.956687532320765,
        #         1.41065433018424,
        #         2.12865961108064,
        #         2.66203156475529,
        #         3.30442073165964,
        #         3.89639860867894,
        #         4.83437593114683,
        #         28.6319224365287
        #     )
        # )
    ),
    contains = "SimulatorRegion"
)

#'
#' Class to simulate Methyl-seq data.
#'
#' @slot nCpG numeric. Number of CpG sites to simulate.
#' @slot pSuccessMethReg numeric. Probability of success in methylated region.
#' @slot pSuccessDemethReg numeric. Probability of success in non methylated region
#' @slot errorMethReg numeric. Error rate in methylated region
#' @slot errorDemethReg numeric. Error rate in methylated region
#' @slot nReadsMethReg numeric. Mean number of reads in methylated region.
#' @slot nReadsDemethReg numeric. Mean number of reads in non methylated regions.
#' @slot phaseDiff numeric. Phase difference in the differentially methylated regions between two samples
#' @slot balanceHypoHyper numeric. Balance of hypo/hyper methylation
#' @slot ratesHMMMatrix numeric. Matrix of values that describes the exponential
#'  decay functions that define the distances between CpG values.
#' @slot distType character. Distribution used to generate replicates:
#' @slot transitionSize numeric.
#' @slot PhiMeth matrix. Transition matrix for CpG locations.
#' @slot PhiDemeth matrix. <Not used>
#' @slot typesLocation numeric. <Not used>
#' @slot returnValue character. Selected column:
#' @slot betaThreshold numeric. Beta threshold value used to calculate M values.
#'
#' @export
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
        maxValue = 1,
        nRepeats = 1,
        idToGene = matrix(),
        pregenerated = TRUE,
        roundDigits = 5,
        WGBSparams = list(
            # Phase difference in the differentially methylated regions between two samples
            phase_diff = c(0, 0.1),
            #error_rate_in_differentially_methylated_region = 0.1
            # File to write results
            outfile = '/tmp/', #output_path
            probs = c(1,1,0.9,0.8,0.7,0.6),
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
        regulatorEffect = c('enhancer', 'repressor'),
        depthAdjust = FALSE
        ),
    contains = "SimulatorRegion"
)
