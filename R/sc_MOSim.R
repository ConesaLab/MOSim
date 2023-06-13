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
      library("SeuratData")
    }
    
    for (omics in omics_types){
      # Load data from seurat
      if(omics == "scRNA-seq"){
        dat <- pbmcMultiome.SeuratData::pbmc.rna
        dat <- subset(x = dat, subset = seurat_annotations %in% c("CD4 TEM", 
                                              "cDC", "Memory B", "Treg"))
        counts <- as.matrix(dat@assays[["RNA"]]@counts)
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
    # Tell the user which celltypes are present in the dataset
    message(paste0("Celltypes in loaded Seurat's PBMC dataset: list('CD4_TEM' = ",
    "c(1:298), 'cDC' = c(299:496), 'Memory_B' = c(497:867), 'Treg' = c(868:1029)"))
    
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
#' to simulate the main dataset
#'
#' @param omics named list containing the omics to simulate as names, which can 
#'    be "scRNA-seq" or "scATAC-seq".
#' @param cellTypes list where the i-th element of the list contains the column 
#'    indices for i-th cell type. List must be a named list.
#' @param diffGenes If number groups > 1, Percentage DE genes to simulate. 
#'    Can be a vector with absolute genes for Up and Down ex: c(250, 500) or a 
#'    percentage for up, down by default: c(0.2, 0.2). The rest will be NE
#' @param minFC Threshold of FC below which are downregulated, by default 0.25
#' @param maxFC Threshold of FC above which are upregulated, by default 4
#' @param numberCells vector of numbers. The numbers correspond to the number 
#'    of cells the user wants to simulate per each cell type. The length of the 
#'       vector must be the same as length of \code{cellTypes}.
#' @param mean vector of numbers of mean per each cell type. Must be specified 
#'    just if \code{numberCells} is specified.
#' @param sd vector of numbers of standard deviation per each cell type. Must be 
#'    specified just if \code{numberCells} is specified.
#' @param noiseGroup OPTIONAL. Number indicating the desired standard deviation
#'      between treatment groups
#' @param group Group for which to estimate parameters
#' @return a named list with simulation parameters for each omics as values.
#'
#' @examples
#' omicsList <- sc_omicData(list("scRNA-seq"))
#' cell_types <- list('CD4_TEM' = c(1:60), 'cDC' = c(299:310), 'Memory_B' = c(497:510), 'Treg' = c(868:900))
#' estimated_params <- sc_param_estimation(omicsList, cell_types)
#' 
sc_param_estimation <- function(omics, cellTypes, diffGenes = c(0.2, 0.2), minFC = 0.25, 
                                maxFC = 4, numberCells = NULL, mean = NULL, 
                                sd = NULL, noiseGroup = 0.5, group = 1){
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
  
  # for (e in 1:N_omics){
  #   # Add variability to groups in comparison to group 1
  #   var_group <- stats::rnorm(nrow(as.data.frame(omics[[e]])), 0, noiseGroup)
  #   # Add variability to each omic compared to group 1
  #   # transform the matrix into TRUE when > 0 false when 0
  #   sim_trueFalse <- (omics[[e]] > 0)
  #   # Multiply the variability vector by a 1/0 to keep the zeros.
  #   for (c in 1:length(colnames(omics[[e]]))){
  #     sim_trueFalse[, c] <- as.integer(as.logical(sim_trueFalse[,c]))
  #     omics[[e]][,c] <- omics[[e]][,c] + (sim_trueFalse[, c] * var_group)
  #   }
  #   omics[[e]] <- abs(omics[[e]])
  # }
  
  # Normalize using scran method, we had to suppress warnings because
  # ATAC has too many zeroes for the normalization to be super comfortable
  norm <- function(om) {
    o <- SingleCellExperiment::SingleCellExperiment(assays=list(counts=as.matrix(om)))
    o <- suppressWarnings(scran::computeSumFactors(o, sizes = seq(10, 100, 5),
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
  
  for(i in 1:N_omics){
    message(paste0("Estimating distribution from original data type: ", i))
    param_est <- SPARSim::SPARSim_estimate_parameter_from_data(raw_data = omics[[i]],
                                                               norm_data = norm_list[[i]],
                                                               conditions = cellTypes)
    param_est_list[[paste0("param_est_", names(omics)[i])]] <- param_est
    
  }
  
  if (group > 1){
    FClist <- list()
    
    for(i in 1:N_omics){
      ## Format diffGenes
      if (diffGenes[1] < 1) {
        ## Relative
        up <- round(diffGenes[1]*length(param_est_list[[i]][[1]][[1]]), digits = 0)
        down <- round(diffGenes[2]*length(param_est_list[[i]][[1]][[1]]), digits = 0)
        
      } else {
        ## Absolute
        up <- diffGenes[1]
        down <- diffGenes[2]
      }
      NE <- length(param_est_list[[i]][[1]][[1]]) - up - down
      message(paste0("Up: ", up, " Down: ", down, " NE: ", NE))
      
      # Now we make the FC vector
      notDE_FCvec <- runif(n = NE, min = minFC + 0.001, max = maxFC - 0.001)
      DE_FCvec <- c(runif(n = up, min = maxFC, max = 100), runif(n = down, 
                                                                 min = 0.0001, max = minFC))
      
      FCvec <- c(DE_FCvec, notDE_FCvec)
      FClist[[paste0("FC_est_", names(omics)[i])]] <- FCvec
      
    } 
  } else if (group == 1){
    FClist <- list()
    
    for(i in 1:N_omics){
      FCvec <- rep(1, length(param_est_list[[i]][[1]][[1]]))
      FClist[[paste0("FC_est_", names(omics)[i])]] <- FCvec
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
        
      if(all_missing){
        libs_param <- param_est_list[[i]][[j]][["lib_size"]]
      } else if (all_specified){
        libs_param <- round(stats::rnorm(n = numberCells[j], mean = mean[j], sd = sd[j]))
      } else if (only_cellnum){
        libs_param <- sample(param_est_list[[i]][[j]][["lib_size"]], 
                             size = numberCells[j], replace = TRUE)
      }
        
      cond_param <- SPARSim::SPARSim_create_simulation_parameter(
        intensity = param_est_list[[i]][[j]][["intensity"]] * FClist[[i]],
        variability = param_est_list[[i]][[j]][["variability"]],
        library_size = libs_param,
        condition_name = param_est_list[[i]][[j]][["name"]],
        feature_names = names(param_est_list[[i]][[j]][["intensity"]]))
        
      cell_type_list[[names(cellTypes)[j]]] <- cond_param
        
    }
    param_est_list_mod[[paste0("param_est_", names(omics)[i])]] <- cell_type_list
  }
    
  return(param_est_list_mod)
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
#'    Can be a vector with absolute genes for Up and Down ex: c(250, 500) or a 
#'    percentage for up, down by default: c(0.2, 0.2). The rest will be NE
#' @param minFC OPTIONAL. Threshold of FC below which are downregulated, by 
#'    default 0.25
#' @param maxFC OPTIONAL. Threshold of FC abofe which are upregulated, by default 4
#' @param numberCells OPTIONAL. Vector of numbers. The numbers correspond to the number of 
#'     cells the user wants to simulate per each cell type. The length of the 
#'         vector must be the same as length of \code{cellTypes}.
#' @param mean OPTIONAL. Vector of numbers of mean per each cell type. Must be specified 
#'     just if \code{numberCells} is specified.The length of the vector must be 
#'         the same as length of \code{cellTypes}.
#' @param sd OPTIONAL. Vector of numbers of standard deviation per each cell type. Must be 
#'     specified just if \code{numberCells} is specified.The length of the vector 
#'         must be the same as length of \code{cellTypes}.
#' @param noiseRep OPTIONAL. Number indicating the desired standard deviation 
#'      between biological replicates.
#' @param noiseGroup OPTIONAL. Number indicating the desired standard deviation
#'      between treatment groups
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
#'                     numberGroups = 2, diffGenes = c(0.2, 0.2), 
#'                     minFC = 0.25, maxFC = 4, numberCells = c(10,20),
#'                     mean = c(2000000, 100000), sd = c(10^3, 10^3), 
#'                     noiseRep = 0.1, noiseGroup = 0.5)
#'
scMOSim <- function(omics, cellTypes, numberReps = 1, numberGroups = 1, 
                    diffGenes = NULL, minFC = 0.25, maxFC = 4,
                    numberCells = NULL, mean = NULL, sd = NULL, noiseRep = 0.1 , 
                    noiseGroup = 0.5){
  # Check for mandatory parameters
  if (missing(omics)){
    stop("You must provide the vector of omics to simulate.")
  }
  
  if (missing(cellTypes)){
    stop("You must provide the correspondence of cells and celltypes")
  }
  
  ## Check that number of groups and number of differentially expressed
  # probabilities makes sense
  if (numberGroups > 1){
    if (is.null(diffGenes) || length(diffGenes) != (numberGroups - 1)){
      stop(paste0("Number of elements in diffGenes must have a length equal to", 
                  " numberGroups -1"))
    }
  }

  seu_groups <- list()
  
  for (g in 1:numberGroups){
    seu_replicates <- list()

    message(paste0("Estimating parameters for experimental group ", g))
    

    param_list <- sc_param_estimation(omics, cellTypes, diffGenes, minFC, maxFC, 
                                      numberCells, mean, sd, noiseGroup, g)
    
    
    for (r in 1:numberReps){
      message(paste0("Simulating parameters for replicate ", r))
      
      N_omics <- length(omics)
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
        seu <- Seurat::CreateSeuratObject(counts = sim_list[[i]],
                                          assay = assay_name,
                                          rownames = rownames(sim_list[[i]]), #explicitly specifying rownames
                                          colnames = colnames(sim_list[[i]])) #and colnames for Seurat obj
        seu_obj[[names(sim_list)[i]]] <- seu
        
      }
      seu_replicates[[paste0("Rep_", r)]] <- seu_obj
    }
    
    seu_groups[[paste0("Group_", g)]] <- seu_replicates
    
  } 
  
  return(seu_groups)
  
}




#' make_cluster_patterns
#' 
#' Function to make the tibble with cluster combinations for the gene expression
#' patterns along the cells
#'
#' @param numcells Number of different celltypes we are simulating
#' @param clusters Number of clusters the user wants to simulate
#'
#' @return A tibble with number of columns equal to number of celltypes, rows
#'  according to the number of TRUE/FALSE combinations corresponding to the
#'  gene expression patterns along the cells
#'
#' @examples
#' patterns <- make_cluster_patterns(numcells = 4, clusters = 4)
#' # patterns <- make_cluster_patterns(numcells = length(cell_types), 
#' # clusters = 8)
#' 
make_cluster_patterns <- function(numcells = 4, clusters = 4){
  
  patterns <- tibble()
  col_names <- paste0("Var", 1:numcells)
  
  for (i in 0:(2^numcells - 1)) {
    boolArr <- vector(mode = "logical", length = numcells)
    
    # Increasing or decreasing depending on which direction
    # you want your array to represent the binary number
    for (j in (numcells - 1):0) {
      boolArr[j + 1] <- as.logical(bitwAnd(i, 2^j))
    }
    
    patterns <- rbind(patterns, boolArr)
  }
  
  colnames(patterns) <- col_names
  # Subset to number of clusters the user wants
  patterns <- dplyr::slice_sample(patterns, n = clusters)
  return(patterns)
}

#' simulate coexpression
#' 
#' Adapted from ACORDE (https://github.com/ConesaLab/acorde) to adapt to our
#' data input type
#'
#' @param feature_ids 
#' @param group_pattern 
#' @param ngroups 
#' @param sim_data 
#'
#' @return the simulated coexpression
#' @export
#'
#' @examples
#' simulate_coexpression(sparsim_sce, 
#' feature_no = 3200, 
#' patterns, cluster_size = 400)
#' 
simulate_coexpression <- function(){}



#' shuffle_group_matrix, Reorder cell type-specific expression matrix during 
#' co-expression simulation. Copied from ACORDE (https://github.com/ConesaLab/acorde)
#' to facilitate stability and running within our scripts
#' 
#' @description This function is used internally by \code{acorde} to perform
#' the shuffling of simulated features for an individual cell type, as part of
#' the co-expression simulation process. The function is called recursively by
#' \code{\link[MOSim:simulate_coexpression]{simulate_coexpression()}} to
#' perform the simulation on a full scRNA-seq matrix.
#'
#' @param sim_data A count matrix with features as rows and cells as columns.
#' Feature IDs must be included in an additional column named \code{feature}.
#' @param feature_ids A two-column \code{tibble} containing \code{top} and \code{bottom}
#' columns, each including the feature IDs of features to be used as highly or
#' lowly expressed when shuffling by the indicated expression pattern.
#' @param group_pattern A logical vector, containing \code{TRUE} to indicate that
#' high expression in that cell type is desired and \code{FALSE} if the opposite.
#' The vector must be ordered as the cell types in \code{sim_data}.
#' @param ngroups An integer indicating the number of groups that top and bottom
#' features should be divided into. It is computed by dividing the number
#' of features selected as highly/lowly expressed by the size of the clusters
#' that are to be generated.
#'
#' @return An expression matrix, with the same characteristics as \code{sim_data},
#' and a number of features defined as the total amount of top/bottom features
#' selected divided by the number of clusters for which co-expression patterns
#' where supplied.

shuffle_group_matrix <- function(sim_data, feature_ids, group_pattern, ngroups){
  
  # select top and bottom features in group
  top <- dplyr::select(feature_ids, top) %>% unlist
  bottom <- dplyr::select(feature_ids, bottom) %>% unlist
  
  # random partitioning of features
  # top
  top.shuffle <- sample(length(top))
  top <- top[top.shuffle]
  top.list <- split(top, cut(seq(1, length(top)), breaks = ngroups, labels = FALSE))
  # bottom
  bottom.shuffle <- sample(length(bottom))
  bottom <- bottom[bottom.shuffle]
  bottom.list <- split(bottom, cut(seq(1, length(bottom)), breaks = ngroups, labels = FALSE))
  
  # bind features following pattern
  features_bound <- vector(mode = "list", length = length(group_pattern))
  # Remove harmless warning
  suppressWarnings(features_bound[group_pattern] <- top.list)
  suppressWarnings(features_bound[!group_pattern] <- bottom.list)
  features_bound <- unlist(features_bound)
  
  # build expression matrix for group
  sim_data.mod <- sim_data %>%
    dplyr::filter(feature %in% features_bound) %>%
    tibble::column_to_rownames("feature")
  sim_data.mod <- sim_data.mod[features_bound,] %>% tibble::rownames_to_column("feature")
  
  return(sim_data.mod)
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
      
      message(paste("the number of totalFeatures you have inserted is higher than what's possible to be generated,", nrow(sim[["scATAC-seq"]]@assays[["ATAC"]]@counts), "peaks were generated instead." ))
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