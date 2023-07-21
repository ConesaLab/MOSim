#' @import SPARSim
#' @import dplyr
#' @import Seurat
#' @import Signac
#' @import stringr
NULL

# Avoid harmless note with R CMD check
if(getRversion() >= "2.15.1")  utils::globalVariables(c(".", "n", "sc_sampleData"))


#' sc_omicData
#'
#' Checks if the user defined data is in the correct format, or loads the default
#' multiomics pbmc dataset from SeuratData package
#'
#' @param omics_types A list of strings which can be either "scRNA-seq" or "scATAC-seq"
#' @param data A user input matrix with genes (peaks in case of scATAC-seq) as 
#'    rows and cells as columns. If a user input matrix is included, cell columns 
#'    must be sorted by cell type. Alternatively scMOSim allows you to use a default 
#'    dataset (PBMC) by not specifying the argument. 
#' @return a named list with omics type as name and the count matrix as value
#' @export
#'
#' @examples
#' # Simulate from PBMC
#' omicsList <- sc_omicData(list("scRNA-seq", "scATAC-seq"))
#' # Simulate using data from the user
#' omicsList_user <- sc_omicData(omics_type = list("scRNA-seq"), data = list(count_matrix)) 
#'
sc_omicData <- function(omics_types, data = NULL){
  # Check for mandatory parameters
  if (missing(omics_types)){
    stop("You must provide the vector of omics to simulate.")
  }
  
  for (omics in omics_types){
    if (omics != "scRNA-seq" && omics != "scATAC-seq"){
      stop("Omics must be a either 'scRNA-seq' or 'scATAC-seq'")
    }
  }
  
  omics_list <- list()
  # If default data
  if (is.null(data)) {
    ## Check we have the dataset installed
    if (SeuratData::AvailableData()["pbmcMultiome.SeuratData","Installed"] != TRUE){
      SeuratData::InstallData("pbmcMultiome.SeuratData")
    } else {
      suppressPackageStartupMessages(library("SeuratData"))
    }
    
    for (omics in omics_types){
      # Load data from seurat
      if(omics == "scRNA-seq"){
        dat <- pbmcMultiome.SeuratData::pbmc.rna
        dat <- subset(x = dat, subset = seurat_annotations %in% c("CD4 TEM", 
                                              "cDC", "Memory B", "Treg"))
        counts <- as.matrix(dat@assays[["RNA"]]@counts)
        # Tell the user which celltypes are present in the dataset
        message(paste0("Celltypes in loaded Seurat's PBMC dataset: list('CD4_TEM' = ",
                       "c(1:298), 'cDC' = c(299:496), 'Memory_B' = c(497:867),", 
                       " 'Treg' = c(868:1029))"))
      } else if (omics == "scATAC-seq"){
        dat <- pbmcMultiome.SeuratData::pbmc.atac
        dat <- subset(x = dat, subset = seurat_annotations %in% c("CD4 TEM", 
                                              "cDC", "Memory B", "Treg"))
        counts <- as.matrix(dat@assays[["ATAC"]]@counts)
      }

      meta <- dat@meta.data$seurat_annotations
      
      metadf <- data.frame("meta" = meta, "cell" = colnames(counts))
      # Sort the metadata according to cell_type
      metadf <- metadf[order(metadf$meta),]
      # Sort by celltype
      counts <- counts[, metadf$cell]
      omics_list[[omics]] <- counts
    }
    
  # If data inputted by user
  } else {
    # If data was inputted by the user, first check
    if (!is.list(data) || length(data) != 1 && length(data) != 2){
      message(paste0("The length of data is ", length(data)))
      stop("Data must be NULL (default) or a list of 1 or 2 elements")
    }
    
    N_data <- length(data)
    for (i in 1:N_data){
      # Then save in a named list
      if (!is.matrix(data[[i]]) && class(data[[i]]) != "Seurat"){
        stop("Each element of data must be either a matrix or a Seurat object")
      } else if (is.matrix(data[[i]])){
        omics_list[[omics_types[[i]]]] <- data[[i]]
      } else if (class(data[[i]]) == "Seurat" && omics_types[[i]] == "scRNA-seq"){
        counts <- as.matrix(data[[i]]@assays[["RNA"]]@counts)
        omics_list[[omics_types[[i]]]] <- counts
      } else if (class(data[[i]]) == "Seurat" && omics_types[[i]] == "scATAC-seq"){
        counts <- as.matrix(data[[i]]@assays[["ATAC"]]@counts)
        omics_list[[omics_types[[i]]]] <- counts
      }
    }
  }
  return(omics_list)
}
  

#' sc_param_estimation
#'
#' Evaluate the users parameters for single cell simulation and use SPARSim
#' to simulate the main dataset. Internal function
#'
#' @param omics named list containing the omics to simulate as names, which can 
#'    be "scRNA-seq" or "scATAC-seq".
#' @param cellTypes list where the i-th element of the list contains the column 
#'    indices for i-th cell type. List must be a named list.
#' @param diffGenes If number groups > 1, Percentage DE genes to simulate.
#'    List of vectors (one per group to compare to group 1) where the vector
#'    contains absolute number of genes for Up and Down ex: c(250, 500) or a 
#'    percentage for up, down ex: c(0.2, 0.2). The rest will be NE
#' @param minFC Threshold of FC below which are downregulated, by default 0.25
#' @param maxFC Threshold of FC above which are upregulated, by default 4
#' @param numberCells vector of numbers. The numbers correspond to the number 
#'    of cells the user wants to simulate per each cell type. The length of the 
#'       vector must be the same as length of \code{cellTypes}.
#' @param mean vector of numbers of mean depth per each cell type. Must be specified 
#'    just if \code{numberCells} is specified.
#' @param sd vector of numbers of standard deviation per each cell type. Must be 
#'    specified just if \code{numberCells} is specified.
#' @param noiseGroup OPTIONAL. Number indicating the desired standard deviation
#'      between treatment groups
#' @param group Group for which to estimate parameters
#' @param genereggroup List with information of genes, clusters and regulators
#'      that must be related to each other
#' @return a list of Seurat object, one per each omic.
#' @return a named list with simulation parameters for each omics as values.
#' @export
#' @examples
#' omicsList <- sc_omicData(list("scRNA-seq"))
#' cell_types <- list('CD4_TEM' = c(1:60), 'cDC' = c(299:310), 'Memory_B' = c(497:510), 'Treg' = c(868:900))
#' estimated_params <- sc_param_estimation(omicsList, cell_types)
#' 
sc_param_estimation <- function(omics, cellTypes, diffGenes = list(c(0.2, 0.2)), 
                                minFC = 0.25, maxFC = 4, numberCells = NULL, 
                                mean = NULL, sd = NULL, noiseGroup = 0.5, 
                                group = 1, genereggroup){
  # Check for mandatory parameters
  if (missing(omics)){
    stop("You must provide the vector of omics to simulate.")
  }
  
  if (missing(cellTypes)) {
    stop("You must provide the correspondence of cells and celltypes")
  }
  
  all_missing <- is.null(numberCells) && is.null(mean) && is.null(sd)
  all_specified <- !is.null(numberCells) && !is.null(mean) && !is.null(sd)
  only_cellnum <- !is.null(numberCells) && is.null(mean) && is.null(sd)
  
  if( !(all_missing || all_specified || only_cellnum)){
    
    stop("The user must either not provide the optional arguments, provide them 
         all or only provide cell numbers")
    
  }
  
  N_omics <- length(omics)
  
  # Add variability to groups in comparison to group 1
  VARlist <- list()
  if (group > 1) {
    for (e in 1:N_omics){
      VARlist[[paste0("Var_", names(omics)[e])]] <- data.frame(matrix(0, 
          ncol = length(colnames(omics[[e]])), nrow = length(rownames(omics[[e]]))))
      var_group <- stats::rnorm(nrow(as.data.frame(omics[[e]])), 0, noiseGroup)
      # Add variability to each omic compared to group 1
      # transform the matrix into TRUE when > 0 false when 0
      sim_trueFalse <- (omics[[e]] > 0)
      # Multiply the variability vector by a 1/0 to keep the zeros.
      for (c in 1:length(colnames(omics[[e]]))){
        sim_trueFalse[, c] <- as.integer(as.logical(sim_trueFalse[,c]))
        omics[[e]][,c] <- omics[[e]][,c] + (sim_trueFalse[, c] * var_group)
        VARlist[[paste0("Var_", names(omics)[e])]][, c] <- (sim_trueFalse[, c] * var_group)
        
      }
      omics[[e]] <- abs(omics[[e]])
    }
  }
  
  
  # Normalize using scran method, we had to suppress warnings because
  # ATAC has too many zeroes for the normalization to be super comfortable
  norm <- function(om) {
    o <- SingleCellExperiment::SingleCellExperiment(assays=list(counts=as.matrix(om)))
    o <- suppressWarnings(scran::computeSumFactors(o, sizes = seq(20, 100, 5),
                                                     positive = FALSE))
    # Apply normalization factors
    o <- scater::normalizeCounts(o, log = FALSE)
    o[is.na(o)] <- 0
    o[!is.finite(o)] <- 0
    o <- abs(o)

    return(o)
  }
  
  norm_list <- lapply(omics, norm)

  param_est_list <- list()
  FC_used_list <- list()
  
  for(i in 1:N_omics){
    message(paste0("Estimating distribution from original data type: ", i))
    param_est <- SPARSim::SPARSim_estimate_parameter_from_data(raw_data = omics[[i]],
                                                               norm_data = norm_list[[i]],
                                                               conditions = cellTypes)
    param_est_list[[paste0("param_est_", names(omics)[i])]] <- param_est
    
  }
  
  if (group > 1){
    # if its from group 2 upwards, we make the differences by multiplying
    # Times a fold change vector, thus we generate it
    FClist <- list()
    if (N_omics > 1){
      prov <- MOSim::make_association_dataframe(group, genereggroup)
      
      associationMatrix <- prov$associationMatrix
      dfPeakNames <- prov$dfPeakNames
      dfGeneNames <- prov$dfGeneNames
    } else {
      
      # Otherwise define the matrix here with the same format
      dfGeneNames <- colnames(omics[[1]])
      dfPeakNames <- NA
      columns <- c("Gene_ID", "Peak_ID", "Effect", "Gene_cluster", "Peak_cluster", 
                   "Gene_DE", "Peak_DE")
      
      associationMatrix <- data.frame(matrix(nrow = length(dfGeneNames), 
                                             ncol = length(columns)))
      associationMatrix["Gene_ID"] <- dfGeneNames
      associationMatrix["Peak_ID"] <- rep(NA, length(dfGeneNames))
      associationMatrix["Effect"] <- rep("NE", length(dfGeneNames))
      clus <- rep(1:length(genereggroup$`Clusters_scRNA-seq`), 
                  each = length(genereggroup$`Clusters_scRNA-seq`[[1]]))
      associationMatrix["Gene_cluster"] <- c(clus, rep(0, length(dfGeneNames) - length(clus)))
      associationMatrix["Peak_cluster"] <- rep(NA, length(dfGeneNames))
      associationMatrix["Gene_DE"] <- c(rep("Up", length(genereggroup[[paste0("GeneExtraUp_G", group)]])),
                                        rep("Down", length(genereggroup[[paste0("GeneExtraDown_G", group)]])),
                                        rep("NE", length(dfGeneNames) - length(genereggroup[[paste0("GeneExtraUp_G", group)]]) - length(genereggroup[[paste0("GeneExtraDown_G", group)]])))
    }
    
    for(i in 1:N_omics){ ## Loop for differentially expressed genes
      if (identical(names(omics[i]), "scRNA-seq")){
        up <- length(associationMatrix[associationMatrix$Gene_DE == "Up",][[1]])
        down <- length(associationMatrix[associationMatrix$Gene_DE == "Down",][[1]])
      } else {
        up <- length(associationMatrix[associationMatrix$Peak_DE == "Up",][[1]])
        down <- length(associationMatrix[associationMatrix$Peak_DE == "Down",][[1]])
      }
      NE <- length(param_est_list[[i]][[1]][[1]]) - up - down
      message(paste0("Up: ", up, " Down: ", down, " NE: ", NE))
      
      # Now we make the FC vector
      notDE_FCvec <- runif(n = NE, min = minFC + 0.001, max = maxFC - 0.001)
      Up_FCvec <- runif(n = length(up), min = maxFC, max = 100)
      Down_FCvec <- runif(n = length(down), min = 0.0001, max = minFC)
      
      ## Here I have 
      
      # The genes affected by FC will be at the begining (also these were
      # The ones most probably included into co-expression patterns)
      
      if (identical(names(omics)[i], "scRNA-seq")){
        FClist[[paste0("FC_est_", names(omics)[i])]] <- order_FC_forMatrix(associationMatrix$Gene_DE, Up_FCvec, Down_FCvec, notDE_FCvec)
      } else {
        FClist[[paste0("FC_est_", names(omics)[i])]] <- order_FC_forMatrix(associationMatrix$Peak_DE, Up_FCvec, Down_FCvec, notDE_FCvec)
      }
    } 
  } else if (group == 1){
    # If its the first group, we dont need to add FC, so we multiply by one
    # Instead of a fold change vector
    FClist <- list()
    VARlist <- list()
    
    for(i in 1:N_omics){
      FCvec <- rep(1, length(param_est_list[[i]][[1]][[1]]))
      FClist[[paste0("FC_", names(omics)[i])]] <- FCvec
      
      # Same for variability
      ## TESTING ONGOING
      VARvec <- rep(1, length(param_est_list[[i]][[1]][["variability"]]))
      VARlist[[paste0("Var_", names(omics)[i])]] <- VARvec
    }
  }

  N_param_est_list<- length(param_est_list)
  N_cellTypes <- length(cellTypes)
  param_est_list_mod <- list()
  
    
  for(i in 1:N_param_est_list){
    cell_type_list <- list()
    message(paste0("Creating parameters for omic: ", i))
      
    for(j in 1:N_cellTypes){
      message(paste0("Creating parameters for cell type: ", j))
      # Estimate library size  
      if(all_missing){
        libs_param <- param_est_list[[i]][[j]][["lib_size"]]
      } else if (all_specified){
        libs_param <- round(stats::rnorm(n = numberCells[j], mean = mean[j], sd = sd[j]))
      } else if (only_cellnum){
        libs_param <- sample(param_est_list[[i]][[j]][["lib_size"]], 
                             size = numberCells[j], replace = TRUE)
      }
      
      # Since we add the foldchange as a scalar, it's multiplied on top of all
      # celltypes, not only one.
      cond_param <- SPARSim::SPARSim_create_simulation_parameter(
        intensity = param_est_list[[i]][[j]][["intensity"]] * FClist[[i]],
        variability = param_est_list[[i]][[j]][["variability"]],
        library_size = libs_param,
        condition_name = param_est_list[[i]][[j]][["name"]],
        feature_names = names(param_est_list[[i]][[j]][["intensity"]]))
        
      cell_type_list[[names(cellTypes)[j]]] <- cond_param
      
      associationMatrix[[names(FClist[[i]])]] <- FClist[[i]]
        
    }
    param_est_list_mod[[paste0("param_est_", names(omics)[i])]] <- cell_type_list
  }
    
  return(list(param_list = param_est_list_mod,
              VARlist = VARlist, 
              associationMatrix = associationMatrix, dfPeakNames = dfPeakNames, 
              dfGeneNames = dfGeneNames))
}



#' scMOSim
#'
#' Performs multiomic simulation of single cell datasets
#'
#' @param omics named list containing the omic to simulate as names, which can 
#'     be "scRNA-seq" or "scATAC-seq".
#' @param cellTypes list where the i-th element of the list contains the column 
#'     indices for i-th experimental conditions. List must be a named list.
#' @param numberReps OPTIONAL. Number of replicates per group
#' @param numberGroups OPTIONAL. number of different groups
#' @param diffGenes OPTIONAL. If number groups > 1, Percentage DE genes to simulate.
#'    List of vectors (one per group to compare to group 1) where the vector
#'    contains absolute number of genes for Up and Down ex: c(250, 500) or a 
#'    percentage for up, down ex: c(0.2, 0.2). The rest will be NE
#' @param minFC OPTIONAL. Threshold of FC below which are downregulated, by 
#'    default 0.25
#' @param maxFC OPTIONAL. Threshold of FC above which are upregulated, by default 4
#' @param numberCells OPTIONAL. Vector of numbers. The numbers correspond to the number of 
#'     cells the user wants to simulate per each cell type. The length of the 
#'         vector must be the same as length of \code{cellTypes}.
#' @param mean OPTIONAL. Vector of numbers of mean depth per each cell type. Must be specified 
#'     just if \code{numberCells} is specified.The length of the vector must be 
#'         the same as length of \code{cellTypes}.
#' @param sd OPTIONAL. Vector of numbers of standard deviation per each cell type. Must be 
#'     specified just if \code{numberCells} is specified.The length of the vector 
#'         must be the same as length of \code{cellTypes}.
#' @param noiseRep OPTIONAL. Number indicating the desired standard deviation 
#'      between biological replicates.
#' @param noiseGroup OPTIONAL. Number indicating the desired standard deviation
#'      between treatment groups
#' @param regulatorEffect OPTIONAL. To simulate relationship scRNA-scATAC, list 
#'      of vectors (one per group) where the vector contains absolute number of
#'      regulators for Activator and repressor ex: c(150, 200) or a percentage
#'      for Activator and repressor ex: c(0.2, 0.1). The rest will be NE. If not
#'      provided, no table of association between scRNA and scATAC is outputted.
#' @param associationList OPTIONAL. A 2 columns dataframe reporting peak ids 
#'      related to gene names. If not provided, no table of association between
#'      scRNA and scATAC is outputted.
#' @param feature_no OPTIONAL. If only scRNA-seq to simulate or scRNA and scATAC
#'      but no regulatory constraints, total number of features to be distributed 
#'      between the coexpression clusters.
#' @param clusters OPTIONAL. Number of co-expression patterns the user wants
#'      to simulate
#' @param cluster_size OPTIONAL. It may be inputted by the user. Recommended: 
#'      by default, its the number of features divided by the number of patterns 
#'      to generate.
#' @return a list of Seurat object, one per each omic.
#' @export
#'
#' @examples
#'
#' cell_types <- list(cellA = c(1:20), cellB = c(161:191))
#' omicsList <- sc_omicData(list("scRNA-seq", "scATAC-seq"))
#' sim <- scMOSim(omicsList, cell_types)
#' # or
#' sim_with_arg <- scMOSim(omicsList, cell_types, numberReps = 2, 
#'                     numberGroups = 2, diffGenes = list(c(0.2, 0.3)), 
#'                     minFC = 0.25, maxFC = 4, numberCells = c(10,20),
#'                     mean = c(2000000, 100000), sd = c(10^3, 10^3), 
#'                     noiseRep = 0.1, noiseGroup = 0.5)
#'
scMOSim <- function(omics, cellTypes, numberReps = 1, numberGroups = 1, 
                    diffGenes = NULL, minFC = 0.25, maxFC = 4,
                    numberCells = NULL, mean = NULL, sd = NULL, noiseRep = 0.1 , 
                    noiseGroup = 0.5, regulatorEffect = NULL, associationList = NULL, 
                    feature_no = 8000, clusters = 8, cluster_size = NULL){
  
  # Check for mandatory parameters
  if (missing(omics)){
    stop("You must provide the vector of omics to simulate.")
  }
  
  if (missing(cellTypes)) {
    stop("You must provide the correspondence of cells and celltypes")
  }
  
  cellTypesCoex <- cellTypes
  if (is.null(numberCells)){
    numberCells <- lengths(cellTypes)
  }

  ## Check that number of groups and number of differentially expressed
  # probabilities makes sense
  if (numberGroups > 1){
    if (is.null(diffGenes) || length(diffGenes) != (numberGroups - 1)){
      stop(paste0("Number of elements in diffGenes must have a length equal to",
                  " numberGroups -1"))
    }
  }
  
  ## Check that columns of association list are c("Peak_ID", "Gene_ID")
  if (!is.null(associationList) && !identical(colnames(associationList), c("Peak_ID", "Gene_ID"))){
    stop("Column names of the user-inputted association list should be 'Peak_ID' and 'Gene_ID'")
  } else if (is.null(associationList)){
    ## Get the association list loaded in the package
    message("Loading default association list from MOSim package")
    data("associationList")
  }
  
  if (is.null(regulatorEffect) && identical(names(omics[2]), "scATAC-seq")){
    stop("You requested a simulation of scATAC-seq data but did not provide information for variable <regulatorEffect>")
  }
  
  ## format number regulators, if its relative, make absolute
  numfeat <- length(associationList$Peak_ID)
  genereggroup <- list()
  
  if (numberGroups > 1){
    if (!is.null(omicsList[2][[1]])){
      for (i in 2:numberGroups){
        
        if (regulatorEffect[[i -1]][1] < 1) {
          # If relative, make absolute numbers
          numActivator <- round(regulatorEffect[[i -1]][1]*numfeat, digits = 0)
          numRepressor <- round(regulatorEffect[[i -1]][2]*numfeat, digits = 0)
        } else {
          numActivator <- regulatorEffect[[i -1]][1]
          numRepressor <- regulatorEffect[[i -1]][2]
          if (numActivator + numRepressor > numfeat){
            stop(paste0("Number of requested Activators and Repressors is bigger
                      than the possible regulators according to the association
                      dataframe. Activators plus repressors must be below: ",
                        numfeat))
          }
        }
        # Get names of regulators per group
        numActivator <- sample(associationList$Peak_ID, numActivator)
        numRepressor <- sample(setdiff(associationList$Peak_ID, numActivator), numRepressor)
        
        ## Get activated, repressed and other diffexp for genes
        if (diffGenes[[i -1]][1] < 1) {
          numup <- round(diffGenes[[i -1]][1]*nrow(omics[[1]]), digits = 0)
          numdown <- round(diffGenes[[i -1]][2]*nrow(omics[[1]]), digits = 0)
        } else {
          numup <- diffGenes[[i -1]][1]
          numdown <- diffGenes[[i -1]][2]
          if (numup + numdown > nrow(omics[[1]])){
            stop(paste0("Number of requested Upregulated and Downregulated genes
                      is bigger than the number of total genes: ", nrow(omics[[1]])))
          }
        }
        
        # Get the genes corresponding to the activator regulators and repressors
        genesActivated <- associationList[associationList$Peak_ID %in% numActivator, ]
        genesRepressed <- associationList[associationList$Peak_ID %in% numRepressor, ]
        # And add other random genes (not in the association list)
        u <- numup - length(genesActivated[[2]])
        d <- numdown - length(genesRepressed[[2]])
        availGenes <- setdiff(rownames(omics[[1]]), as.vector(associationList[[2]]))
        genesUp <- sample(availGenes, u)
        availGenes <- setdiff(availGenes, genesUp)
        genesDown <- sample(availGenes, d)
        
        remaining <-  setdiff(rownames(omics[[1]]), genesActivated[[2]])
        remaining <- setdiff(remaining, genesRepressed[[2]])
        remaining <- setdiff(remaining, genesUp)
        remaining <- setdiff(remaining, genesDown)
        genereggroup[[paste0("GeneActivated_G", i)]] <- genesActivated
        genereggroup[[paste0("GeneRepressed_G", i)]] <- genesRepressed
        genereggroup[[paste0("GeneExtraUp_G", i)]] <- genesUp
        genereggroup[[paste0("GeneExtraDown_G", i)]] <- genesDown
        genereggroup[[paste0("GeneRemaining_G", i)]] <- remaining
        
        # Here the up and down, have to be as many as possible from the association
        # list, so take rows where its activator and make them up, when its repressor
        # make them down
        if (numup < length(numActivator) || numdown < length(numRepressor)){
          stop("You have asked for many regulators, but there aren't enough 
             differentially expressed genes to be regulated")
        }
        # Get the info for the ATAC
        u <- numup - length(unique(genesActivated[[1]]))
        d <- numdown - length(unique(genesRepressed[[1]]))
        availFeat <- setdiff(rownames(omics[[2]]), unique(associationList[[1]]))
        regAct <- sample(availFeat, u)
        availFeat <- setdiff(availFeat, regAct)
        regRep <- sample(availFeat, d)
        
        remaining <- setdiff(rownames(omics[[2]]), as.vector(genesActivated[[1]]))
        remaining <- setdiff(remaining, as.vector(genesRepressed[[1]]))
        remaining <- setdiff(remaining, regAct)
        remaining <- setdiff(remaining, regRep)
        
        genereggroup[[paste0("FeatExtraUp_G", i)]] <- regAct
        genereggroup[[paste0("FeatExtraDown_G", i)]] <- regRep
        genereggroup[[paste0("FeatRemaining_G", i)]] <- remaining
        
      }
    } else {
      ## Here say what to do if we don't have ATAC-seq
      for (i in 2:numberGroups){
        if (diffGenes[[i -1]][1] < 1) {
          numup <- round(diffGenes[[i -1]][1]*nrow(omics[[1]]), digits = 0)
          numdown <- round(diffGenes[[i -1]][2]*nrow(omics[[1]]), digits = 0)
        } else {
          numup <- diffGenes[[i -1]][1]
          numdown <- diffGenes[[i -1]][2]
          if (numup + numdown > nrow(omics[[1]])){
            stop(paste0("Number of requested Upregulated and Downregulated genes
                        is bigger than the number of total genes: ", nrow(omics[[1]])))
          }
        }
        u <- sample(rownames(omics[[1]]), numup)
        availGenes <- setdiff(rownames(omics[[1]]), u)
        d <- sample(availGenes, numdown)
        genereggroup[[paste0("GeneExtraUp_G", i)]] <- u
        genereggroup[[paste0("GeneExtraDown_G", i)]] <- d
        genereggroup[[paste0("GeneRemaining_G", i)]] <- setdiff(availGenes, d)
      }
    }
  }
  
  ### Start working on the data
  
  N_omics <- length(omics)
  
  ## Subset only columns of interest (our celltypes)
  for (om in 1:N_omics){
    omics[[om]] <- omics[[om]][, unname(unlist(cellTypes))]
  }
  # Reorganize our cellTypes variable according to the subset
  createRangeList <- function(vector, nam) {
    rangeList <- list()
    totalValues <- sum(vector)
    start <- 1
    
    for (i in 1:length(vector)) {
      end <- start + vector[i] - 1
      rangeList[[i]] <- seq(start, end)
      start <- end + 1
    }
    names(rangeList) <- nam
    return(rangeList)
  }
  
  cellTypes <- createRangeList(numberCells, names(cellTypes))
  
  
  # Make the patterns to simulate coexpression
  lpatterns <- make_cluster_patterns(length(cellTypes), clusters = clusters)
  # Get also the indices of cluster patterns that are opposite and will be 
  # important for the regulators
  genereggroup[["opposite_indices"]] <- lpatterns$opposite_indices
  patterns <- lpatterns$patterns
  
  # Simulate coexpression
  for (i in 1:N_omics){
    coexpr_results <- MOSim::simulate_coexpression(omics[[i]],
                                            feature_no = feature_no, cellTypes = cellTypesCoex,
                                            patterns = patterns, cluster_size = cluster_size)
    
    # Get the coexpressed matrix out
    omics[[i]] <- as.data.frame(coexpr_results$sim_matrix)
    rownames(omics[[i]]) <- omics[[i]]$feature
    omics[[i]]$feature <- NULL
    # Get the clusters out
    genereggroup[[paste0("Clusters_", names(omics)[[i]])]] <- coexpr_results$sim_clusters
  }

  # Start the lists we will need to include in the output
  seu_groups <- list()
  Association_list <- list()
  VAR_used_list <- list()
  param_list <- list()
  
  for (g in 1:numberGroups){
    seu_replicates <- list()

    message(paste0("Estimating parameters for experimental group ", g))
    

    param_l <- MOSim::sc_param_estimation(omics, cellTypes, diffGenes, minFC, maxFC, 
                                          numberCells, mean, sd, noiseGroup, g, 
                                          genereggroup)
    
    
    Association_list[[paste0("AssociationMatrix_Group_", g)]] <- param_l$AssociationMatrix
    VAR_used_list[[paste0("VAR_Group_", g)]] <- param_l$VARlist
    param_list <- param_l$param_list
    
    for (r in 1:numberReps){
      message(paste0("Simulating parameters for replicate ", r))
      
      sim_list <- list()
      
      
      for(i in 1:N_omics){
        # Simulate the replicate
        sim <- SPARSim::SPARSim_simulation(dataset_parameter = param_list[[i]])
        sim <- sim[["count_matrix"]]
        # Generate a standard deviation to add to the matrix
        var_rep <- stats::rnorm(nrow(as.data.frame(omics[[i]])), 0, noiseRep)
        # Add the standard deviation of the replicate
        # transform the matrix into TRUE when > 0 false when 0
        sim_trueFalse <- (sim > 0)
        # Multiply the variability vector by a 1/0 to keep the zeros.
        for (e in 1:length(colnames(sim))){
          sim_trueFalse[, e] <- as.integer(as.logical(sim_trueFalse[,e]))
          sim[,e] <- sim[,e] + (sim_trueFalse[,e] * var_rep)
        }
        # Make sure there are no negative numbers
        sim <- abs(sim)
        # Now pass this modified matrix back
        sim_list[[paste0("sim_", names(omics)[i])]] <- sim
        
      }
      
      seu_obj <- list()
      N_sim <- length(sim_list)
      
      for(i in 1:N_sim){
        
        assay_name <- str_split(names(sim_list)[i], "-")[[1]][1]
        assay_name <- sub("sim_sc","",assay_name)
        
        # Rename the 
        if(identical(names(sim_list[[i]]), "sim_scRNA-seq")){
          newNammes <- param_l$dfGeneNames
        } else{
          newNames <- param_l$dfPeakNames
        }
        
        seu <- Seurat::CreateSeuratObject(counts = sim_list[[i]],
                                          assay = assay_name,
                                          rownames = newNames, #explicitly specifying rownames
                                          colnames = colnames(sim_list[[i]])) #and colnames for Seurat obj
        seu_obj[[names(sim_list)[i]]] <- seu
        
      }
      seu_replicates[[paste0("Rep_", r)]] <- seu_obj
    }
    
    seu_groups[[paste0("Group_", g)]] <- seu_replicates
    
  } 
  
  # Bring back final celltypes
  seu_groups[["cellTypes"]] <- cellTypes
  seu_groups[["patterns"]] <- patterns
  # Bring back FC vector
  seu_groups[["AssociationMatrices"]] <- Association_list
  seu_groups[["Variability"]] <- VAR_used_list
  
  return(seu_groups)
  
}


#' sc_omicSim
#'
#' Defines regulatory functions of features for regulatory omics
#'
#' @param sim named list containing the omic simulated as names ("sim_scRNA-seq" and 
#'     "sim_scATAC-seq") , and seurat objects as values.
#' @param cellTypes list where the i-th element of the list contains the column 
#'     indices for i-th experimental conditions. List must be a named list.
#' @param totalFeatures OPTIONAL. Numeric value. Total number of features for 
#'     the regulatory omic.If not provided it uses the same amount of features 
#'         used for the simulation.
#' @param regulatorEffect OPTIONAL. Named list of length 3 where the user can 
#'     pass the percentage of activators, repressors and NE he wants as output. 
#'         If not provided the function outputs the dataframe without sub-setting 
#'             it according to percentages.
#' @param associationList OPTIONAL. A 2 columns dataframe reporting peak ids and 
#'     gene names. If not provided the code uses our own associationlist derived 
#'         from hg19.
#' @return named list containing a 3 columns dataframe (peak_id, activity, 
#'     cell_types), one per each couple of \code{cellTypes}.
#' @export
#'
#' @examples
#'
#'
#' omicsList <- sc_omicData(c("scRNA-seq","scATAC-seq"))
#' cell_types <- list(cellA = c(1:20), cellB = c(161:191))
#' sim <-scMOSim(omicsList, cell_types, numberCells = c(10,20), 
#'             mean = c(2*10^6, 2*10^3), sd = c(10^3, 10^2))
#' cell_types <- list(cellA= c(1:10), cellB = c(11:30))
#' regulatorEffect = list('activator' = 0.8,'repressor' = 0.1,'NE' = 0.1)
#' sc_omicSim(sim, cell_types, totalFeatures = 500, regulatoreEffect = regulatorEffect)
#'
sc_omicSim <- function(sim, cellTypes, totalFeatures = NULL, 
                       regulatorEffect = NULL, associationList = NULL ){
  # Check for mandatory parameters
  if (missing(sim)){
    stop("You must provide the named list of simulated omics.")
  }
  
  if (missing(cellTypes)){
    stop("You must provide the correspondence of cells and celltypes")
  }
  
  if(is.null(totalFeatures)){
    
    atac_counts <- sim[["sim_scATAC-seq"]]@assays[["ATAC"]]@counts
    
  } else if (is.numeric(totalFeatures)){
    
    if(totalFeatures <= nrow(sim[["sim_scATAC-seq"]]@assays[["ATAC"]]@counts)){
      
      # Set the number of features you want to select
      num_features <- totalFeatures
      # Get the data from the ATAC assay of the scATAC-seq object
      atac_data <- sim[["sim_scATAC-seq"]]@assays[["ATAC"]]@counts
      # Subset the ATAC data to retain only the selected number of features
      atac_data <- atac_data[sample(nrow(atac_data), num_features, replace = FALSE), ]
      # Update the ATAC assay with the subsetted data
      sim[["sim_scATAC-seq"]]@assays[["ATAC"]]@counts<- atac_data
      atac_counts <- sim[["sim_scATAC-seq"]]@assays[["ATAC"]]@counts
      
    } else if (totalFeatures > nrow(sim[["sim_scRNA-seq"]]@assays[["RNA"]]@counts)){
      
      message(paste("the number of totalFeatures you have inserted is higher ", 
                    "than what's possible to be generated,", 
                    nrow(sim[["scATAC-seq"]]@assays[["ATAC"]]@counts), 
                    " peaks were generated instead." ))
      atac_counts <- sim[["sim_scATAC-seq"]]@assays[["ATAC"]]@counts
      
    }
  }
  
  rna_counts <- sim[["sim_scRNA-seq"]]@assays[["RNA"]]@counts
  
  #calculate gene expression for each cellTypes
  gene_expression_list <- lapply(names(cellTypes), function(cell_type) {
    
    gene_expression <- rowSums(rna_counts[, cellTypes[[cell_type]]])
    names(gene_expression) <- rownames(rna_counts)
    return(gene_expression)
    
  })
  
  names(gene_expression_list) <- paste0("gene_expr_", names(cellTypes))
  
  if(length(cellTypes) > 2){
    
    message(paste0("the length of cellTypes is ", length(cellTypes)))
    
    da_peaks_atac_list <- lapply(seq_along(cellTypes), function(i) {
      
      j <- i %% length(cellTypes) + 1
      ident1 <- names(cellTypes)[i]
      ident2 <- names(cellTypes)[j]
      Seurat::FindMarkers(object = sim[["sim_scATAC-seq"]], ident.1 = ident1, ident.2 = ident2, min.pct = 0.05)
      
    })
    
    names(da_peaks_atac_list) <- paste0("markers_", names(cellTypes), "_", c(names(cellTypes)[-1], names(cellTypes)[1]))
    
  } else if(length(cellTypes) ==2 ){
    
    da_peaks_atac_list <- list()
    da_peaks_atac <- Seurat::FindMarkers(object = sim[["sim_scATAC-seq"]], ident.1 = names(cellTypes[1]), ident.2 = names(cellTypes[2]) , min.pct = 0.05)
    da_peaks_atac_list[[paste0("markers_", names(cellTypes[1]), "_", names(cellTypes[2]) )]] <- da_peaks_atac
    
  }
  
  #subselecting the dataframe according to upregulated peaks in celltype 1 and 2
  
  subset_list <- lapply(seq_along(da_peaks_atac_list), function(i) {
    
    names_i <- strsplit(names(da_peaks_atac_list), "_")[[i]]
    cell_type_1 <- names_i[2]
    cell_type_2 <- names_i[3]
    
    x <- da_peaks_atac_list[[i]]
    subset_up_cell1 <- x[x$avg_log2FC > 0,]
    subset_up_cell2 <- x[x$avg_log2FC < 0,]
    
    upreg_1 <- rownames(subset_up_cell1)
    upreg_2 <- rownames(subset_up_cell2)
    
    subset_list_upreg <- list(upreg_1, upreg_2)
    names(subset_list_upreg) <- c(paste0("up_reg_", cell_type_1),paste0("up_reg_", cell_type_2))
    subset_list_upreg
    
  })
  names(subset_list) <- paste0(names(da_peaks_atac_list))
  
  #loading human association list
  if (is.null(associationList)){
    
    association_list <- as.data.frame(sc_sampleData$seurat_association_list)
    
  } else if (is.list(associationList)){
    
    if (length(associationList) == 2){
      
      association_list <- associationList
      
    } else {
      
      message("the associationList must be a dataframe having two columns reporting peak ids and gene names")
      return(NA)
      
    }
  }
  
  result_list <- lapply(seq_along(subset_list), function(j) {
    #open empty dataframe
    peak_df <- data.frame(peak_id = character(), 
                          activity = character(), 
                          cell_type = character(), 
                          stringsAsFactors = FALSE)
    
    for(i in 1:nrow(association_list)) {
      
      peak <- association_list[i,1]
      gene <- association_list[i,2]
      
      # Check if the gene is in the rownames of the rna_counts matrix
      if(gene %in% rownames(rna_counts)) {
        activity <- "NE"
        
        # Check if the peak is upregulated in cell type 1 and the gene expression of gene i for cell type 1 is >0
        if(peak %in% subset_list[[j]][[1]] && gene_expression_list[[j]][gene] > 0) {
          
          activity <- "activator"
          peak_df <- rbind(peak_df, data.frame(peak_id = peak, activity = activity, cell_type = strsplit(names(subset_list[[j]]),"_")[[1]][3], stringsAsFactors = FALSE))
          
          #check if the peak is upregulated in cell type 1 and the gene expression of gene i for cell type 1 is ==0
          #or
          # if the peak is downregulated in cell type 1 (up in 2) and the gene expression of gene i for cell type 1 is >0
        } else if((peak %in% subset_list[[j]][[1]] && gene_expression_list[[j]][gene] == 0) || 
                  (peak %in% subset_list[[j]][[2]] && gene_expression_list[[j]][gene] > 0)) {
          
          activity <- "repressor"
          peak_df <- rbind(peak_df, data.frame(peak_id = peak, activity = activity, cell_type = strsplit(names(subset_list[[j]]),"_")[[1]][3], stringsAsFactors = FALSE))
          
        }
        
        next_j<- (j %% length(subset_list))+1
        
        #check if the peak is upregulated in cell type 2 and the gene expression of gene i for cell type 2 is >0
        if(peak %in% subset_list[[j]][[2]] && gene_expression_list[[next_j]][gene] > 0) {
          
          activity <- "activator"
          peak_df <- rbind(peak_df, data.frame(peak_id = peak, activity = activity, cell_type = strsplit(names(subset_list[[j]]),"_")[[2]][3], stringsAsFactors = FALSE))
          
          #check if the peak is upregulated in cell type 2 and the gene expression of gene i for cell type 2 is ==0
          #or
          #check if the peak is downregulated in cell type 2 (up in 1) and the gene expression of gene i for cell type 2 is >0
        } else if((peak %in% subset_list[[j]][[2]] && gene_expression_list[[next_j]][gene] == 0) ||
                  (peak %in% subset_list[[j]][[1]] && gene_expression_list[[next_j]][gene] > 0)) {
          
          activity <- "repressor"
          peak_df <- rbind(peak_df, data.frame(peak_id = peak, activity = activity, cell_type = strsplit(names(subset_list[[j]]),"_")[[2]][3], stringsAsFactors = FALSE))
          
        }
        
        if (activity == "NE") {
          peak_df <- rbind(peak_df, data.frame(peak_id = peak, activity = activity, cell_type = "NA", stringsAsFactors = FALSE))
        } 
      } 
      
    }
    
    peak_df
    
  })
  
  names(result_list) <- paste0(names(da_peaks_atac_list))
  
  #sub-setting the dataframes contained in result list according to regulatorEffect percentages
  
  #if regulatorEffect is null then it returns the dataframe without subsetting
  if (is.null(regulatorEffect)){
    
    return(result_list)
    
    
  } # if regulatorEffect is not null then it return the sub-setted dataframe according to percentages passed
  
  else if(typeof(regulatorEffect) =="list"){
    
    if(length(regulatorEffect) == 3){
      
      subset_df <- function(df) {
        #take dataframe of activators
        activator_df <- df[df$activity == 'activator', ]
        #multiplies the n* of rows for the input percentage
        n_activator <- round(nrow(activator_df) * regulatorEffect[['activator']])
        #create new df randomly containing n_activators activators
        activator_subset <- activator_df[sample(nrow(activator_df), n_activator), ]
        
        #take dataframe of repressors
        repressor_df <- df[df$activity == 'repressor', ]
        #multiplies the n* of rows for the input percentage
        n_repressor <- round(nrow(repressor_df) * regulatorEffect[['repressor']])
        #create new df randomly containing n_repressors repressors
        repressor_subset <- repressor_df[sample(nrow(repressor_df), n_repressor), ]
        
        #take dataframe of NE
        NE_df <- df[df$activity == 'NE', ]
        #multiplies the n* of rows for the input percentage
        n_NE <- round(nrow(NE_df) * regulatorEffect[['NE']])
        #create new df randomly containing n_NE NE
        NE_subset <- NE_df[sample(nrow(NE_df), n_NE), ]
        
        # Combine the subsets and return the df
        subset_df <- rbind(activator_subset, repressor_subset, NE_subset)
        return(subset_df)
      }
      
      # Apply the function to each dataframe in result_list and store the results in a new named list
      sub_result_list <- lapply(result_list, subset_df)
      names(sub_result_list) <- paste0(names(result_list), "_subset")
      return(sub_result_list)
      
    } else {
      
      message("regulatorEffect must be a named list with percentages as values of length 3. The names must be 'activator', 'repressor', 'NE' ")
      return(NA)
      
    }
    
  }
  
}