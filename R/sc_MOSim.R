#' @import SPARSim
#' @import dplyr
#' @import Seurat
#' @import Signac
#' @import stringr
#' @import SeuratData
NULL

# Avoid harmless note with R CMD check
if(getRversion() >= "2.15.1")  utils::globalVariables(c(".", "n", "sc_sampleData",
                                              "feature", "seurat_annotations"))


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
#' count <- omicsList[["scRNA-seq"]]
#' options(Seurat.object.assay.version = "v3")
#' Seurat_obj <- Seurat::CreateAssayObject(counts = count, assay = 'RNA')
#' omic_list_user <- sc_omicData(c("scRNA-seq"), c(Seurat_obj))
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
      message("Installing pbmcMultiome dataset from SeuratData. This may take a while...")
      options(timeout = 3000); SeuratData::InstallData("pbmcMultiome.SeuratData")
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
      if (!is.matrix(data[[i]]) && class(data[[i]]) != "Assay"){
        stop("Each element of data must be either a matrix or a Seurat object")
      } else if (is.matrix(data[[i]])){
        omics_list[[omics_types[[i]]]] <- data[[i]]
      } else if ("Seurat" %in% class(data[[i]]) && omics_types[[i]] == "scRNA-seq"){
        counts <- as.matrix(data[[i]]@assays[["RNA"]]@counts)
        omics_list[[omics_types[[i]]]] <- counts
      } else if ("Seurat" %in% class(data[[i]]) && omics_types[[i]] == "scATAC-seq"){
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
#' cell_types <- list('CD4_TEM' = c(1:60), 'cDC' = c(299:310), 
#'     'Memory_B' = c(497:510), 'Treg' = c(868:900))
#' #estimated_params <- sc_param_estimation(omicsList, cell_types)
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
  
  #We cant move this out of the function because it takes the variability
  # of the group into account
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
      
      # If only scRNA define the matrix here with the same format
      dfGeneNames <- rownames(omics[[1]])
      dfPeakNames <- NA
      columns <- c("Gene_ID", "Peak_ID", "RegulatorEffect", "Gene_cluster", "Peak_cluster", 
                   "Gene_DE", "Peak_DE")
      
      associationMatrix <- data.frame(matrix(nrow = length(dfGeneNames), ncol = length(columns)))
      colnames(associationMatrix) <- columns
      
      associationMatrix["Gene_ID"] <- dfGeneNames
      associationMatrix["Peak_ID"] <- rep(NA, length(dfGeneNames))
      associationMatrix["RegulatorEffect"] <- rep("NE", length(dfGeneNames))
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
        up <- length(associationMatrix[associationMatrix$Gene_DE %in% "Up",][[1]])
        down <- length(associationMatrix[associationMatrix$Gene_DE %in% "Down",][[1]])
      } else {
        up <- length(associationMatrix[associationMatrix$Peak_DE %in% "Up",][[1]])
        down <- length(associationMatrix[associationMatrix$Peak_DE %in% "Down",][[1]])
      }
      NE <- length(param_est_list[[i]][[1]][[1]]) - up - down
      message(paste0("Up: ", up, " Down: ", down, " NE: ", NE))
      
      # Now we make the FC vector
      notDE_FCvec <- runif(n = NE, min = minFC + 0.001, max = maxFC - 0.001)
      Up_FCvec <- runif(n = up, min = maxFC, max = 100)
      Down_FCvec <- runif(n = down, min = 0.0001, max = minFC)
      
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
      VARvec <- rep(1, length(param_est_list[[i]][[1]][["variability"]]))
      VARlist[[paste0("Var_", names(omics)[i])]] <- VARvec
      
      dfPeakNames <- NA
      if (identical(names(omics[i]), "scRNA-seq")){
        dfGeneNames <- rownames(omics[[i]])
      } else {
        dfPeakNames <- rownames(omics[[i]])
      }
    }
    columns <- c("Gene_ID", "Peak_ID", "RegulatorEffect", "Gene_cluster", "Peak_cluster", 
                 "Gene_DE", "Peak_DE")
    
    
    if (length(omics) > 1){
      associationMatrix <- data.frame(matrix(nrow = length(dfPeakNames) + length(dfGeneNames), ncol = length(columns)))
      colnames(associationMatrix) <- columns
      associationMatrix["Gene_ID"] <- c(dfGeneNames, rep(NA, length(dfPeakNames)))
      associationMatrix["Peak_ID"] <- c(rep(NA, length(dfGeneNames)), dfPeakNames)
      associationMatrix["RegulatorEffect"] <- rep("NE", length(associationMatrix["Gene_ID"]))
      clus <- rep(1:length(genereggroup$`Clusters_scRNA-seq`), 
                  each = length(genereggroup$`Clusters_scRNA-seq`[[1]]))
      associationMatrix["Gene_cluster"] <- c(clus, rep(0, length(dfGeneNames) - length(clus)), rep(NA, length(dfPeakNames)))
      associationMatrix["Peak_cluster"] <- c(rep(NA, length(dfGeneNames)), clus, rep(0, length(dfPeakNames)- length(clus)))
    } else {
      associationMatrix <- data.frame(matrix(nrow = length(dfGeneNames), ncol = length(columns)))
      colnames(associationMatrix) <- columns
      associationMatrix["Gene_ID"] <- dfGeneNames
      associationMatrix["Peak_ID"] <- rep(NA, length(dfGeneNames))
      associationMatrix["RegulatorEffect"] <- rep("NE", length(dfGeneNames))
      clus <- rep(1:length(genereggroup$`Clusters_scRNA-seq`), 
                  each = length(genereggroup$`Clusters_scRNA-seq`[[1]]))
      associationMatrix["Gene_cluster"] <- c(clus, rep(0, length(dfGeneNames) - length(clus)))
      associationMatrix["Peak_cluster"] <- rep(NA, length(dfGeneNames))
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
        intensity = param_est_list[[i]][[j]][["intensity"]] * as.numeric(FClist[[i]]),
        variability = param_est_list[[i]][[j]][["variability"]],
        library_size = libs_param,
        condition_name = param_est_list[[i]][[j]][["name"]],
        feature_names = names(param_est_list[[i]][[j]][["intensity"]]))
        
      cell_type_list[[names(cellTypes)[j]]] <- cond_param
        
    }
    param_est_list_mod[[paste0("param_est_", names(omics)[i])]] <- cell_type_list
  }
    
  return(list(param_list = param_est_list_mod, FClist = FClist,
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
#' @param associationList REQUIRED A 2 columns dataframe reporting peak ids 
#'      related to gene names. If user doesnt have one, load from package
#'      data("associationList")
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
#' omicsList <- sc_omicData(list("scRNA-seq"))
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
    associationList <- as.data.frame(MOSim::associationList)
  }
  
  if (is.null(regulatorEffect) && identical(names(omics[2]), "scATAC-seq")){
    stop("You requested a simulation of scATAC-seq data but did not provide information for variable <regulatorEffect>")
  }
  
  ## Message the experimental dessign
  message(paste0("The experimental design includes: 
                 - ", numberReps, " Biological replicates
                 - ", numberGroups, " Experimental groups
                 - ", paste0(diffGenes), " Differentially expressed genes (Up/Down) per group
                 - ", minFC, " FC below which a gene is downregulated
                 - ", maxFC, " FC above which a gene is upregulated
                 - ", paste0(numberCells), " Number of cells per celltype
                 - ", paste0(regulatorEffect), " Regulators (activator/repressor) per group
                 - ", clusters, " Gene co-expression patterns")[1])
  
  ## format number regulators, if its relative, make absolute
  numfeat <- length(associationList$Gene_ID)
  genereggroup <- list()
  
  if (numberGroups > 1){
    if (length(omics) > 1){
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
        u <- abs(numup - length(genesActivated[[2]]))
        d <- abs(numdown - length(genesRepressed[[2]]))
        
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
  FC_used_list <- list()
  VAR_used_list <- list()
  param_list <- list()
  
  for (g in 1:numberGroups){
    seu_replicates <- list()

    message(paste0("Estimating parameters for experimental group ", g))
    
    param_l <- MOSim::sc_param_estimation(omics, cellTypes, diffGenes, minFC, maxFC, 
                                          numberCells, mean, sd, noiseGroup, g, 
                                          genereggroup)
    

    Association_list[[paste0("AssociationMatrix_Group_", g)]] <- param_l$associationMatrix
    FC_used_list[[paste0("FC_Group_", g)]] <- param_l$FClist
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
        
        options(Seurat.object.assay.version = "v3")
        seu <- Seurat::CreateAssayObject(counts = sim_list[[i]],
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
  seu_groups[["FC"]] <- FC_used_list
  seu_groups[["AssociationMatrices"]] <- Association_list
  seu_groups[["Variability"]] <- VAR_used_list
  
  return(seu_groups)
  
}



#' scOmicSettings
#'
#' @param sim a simulated object from scMOSim function
#'
#' @return list of Association matrices explaining the effects of each
#' regulator to each gene
#' @export
#'
#' @examples
#'
#' cell_types <- list(cellA = c(1:20), cellB = c(161:191))
#' omicsList <- sc_omicData(list("scRNA-seq"))
#' sim <- scMOSim(omicsList, cell_types)
#' res <- scOmicSettings(sim)
scOmicSettings <- function(sim){
  asma <- sim$AssociationMatrices
  FC <- sim$FC
  
  # Function to add two new columns to a matrix based on conditions
  add_new_columns <- function(matrix_A, vector1, vector2) {
    new_col1 <- character(nrow(matrix_A))
    new_col2 <- character(nrow(matrix_A))
    vector_index1 <- 1
    vector_index2 <- 1
    for (i in 1:nrow(matrix_A)) {
      if (!is.na(matrix_A[i, 1])) {
        new_col1[i] <- vector1[vector_index1]
        vector_index1 <- vector_index1 + 1
      }
      if (!is.na(matrix_A[i, 2])) {
        new_col2[i] <- vector2[vector_index2]
        vector_index2 <- vector_index2 + 1
      }
    }
    matrix_A <- cbind(matrix_A, new_col1, new_col2)
    return(matrix_A)
  }
  
  if (identical(names(sim$Group_1$Rep_1), c("sim_scRNA-seq", "sim_scATAC-seq"))){
    # Update List A matrices with new columns from List B vectors
    for (i in 1:length(asma)) {
      asma[[i]] <- add_new_columns(asma[[i]], FC[[i]][[1]], FC[[i]][[2]])
      colnames(asma[[i]]) <- c("Gene_ID", "Peak_ID", "RegulatorEffect", "Gene_cluster", 
                               "Peak_cluster", "Gene_DE", "Peak_DE", "Gene_FC", 
                               "Peak_FC")
    } 
  }else{
    for (i in 1:length(asma)){
      asma[[i]]["Gene_FC"] <- FC[[i]][[1]]
    }
  }
  
  return(asma)
}

#' scOmicResults
#'
#' @param sim a simulated object from scMOSim function
#'
#' @return list of seurat objects with simulated data
#' @export
#' @examples
#'
#' cell_types <- list(cellA = c(1:20), cellB = c(161:191))
#' omicsList <- sc_omicData(list("scRNA-seq"))
#' sim <- scMOSim(omicsList, cell_types)
#' res <- scOmicResults(sim)

scOmicResults <- function(sim){
  return(sim[grepl("Group_", names(sim))])
}
