#' @import SPARSim
#' @import dplyr
#' @import Seurat
#' @import Signac
#' @import stringr
NULL

# Avoid note with R CMD check
if(getRversion() >= "2.15.1")  utils::globalVariables(c(".", "n", "sc_sampleData"))


#' sc_omicData
#' 
#' Checks if the user defined data is in the correct format, or loads
#' the pbmc dataset from https://satijalab.org/seurat/articles/pbmc3k_tutorial.html
#' @param omic A string which can be either "scRNA-seq" or "scATAC-seq"
#' @param data A user input matrix with genes (peaks in case of scATAC-seq) 
#'    as rows and cells as columns. Alternatively MOSim allows user to estimate 
#'    the input parameters from an existing count table by typing 'example_matrix'
#' @return a named list with omics type as name and the count matrix as value
#' @export
#' 
#' @examples
#' 
#' scRNAseq <- sc_omicData("scRNA-seq")
#' scATACseq <- sc_omicData("scATAC-seq")
#' scRNAseq_user <- sc_omicData("scRNA-seq", count_matrix) 

sc_omicData <- function(omics, data = NULL){
  
  if (omics != "scRNA-seq" && omics != "scATAC-seq"){
    
    print("omics must be a either 'scRNA-seq' or 'scATAC-seq'")
    return(NA)
    
  }
  
  if (is.null(data)){ 
    load(file='data/sc_sampleData.rda')
    
    if (omics == "scRNA-seq"){ 
      
      ##scRNA##
      omics_list <- list("scRNA-seq" = rna_orig_counts)
      return(omics_list)
      
    } else if (omics =="scATAC-seq"){
      
      ##scATAC##
      omics_list <- list("scATAC-seq" = atac_orig_counts)
      return(omics_list)
      
    } else if(omics == c("scRNA-seq", "scATAC-seq")){
      
      omics_list <- list("scRNA-seq" = rna_orig_counts, "scATAC-seq" = atac_orig_counts)
      return(omics_list)
    }
  } 
  
  if (! is.matrix(data) && class(data) != "Seurat"){
    
    print("data must be either matrix or a Seurat object")
    return(NA)
    
  } else if (is.matrix(data)){
    
    omics_list <- list()
    omics_list[[omics]] <- data 
    return(omics_list)
    
  } else if (class(data) == "Seurat" && omics == "scRNA-seq"){
    
    omics_list <- list()
    counts <- data@assays[["RNA"]]@counts
    counts_matrix <- as.matrix(counts)
    omics_list[[omics]] <- counts_matrix
    return(omics_list)
    
  } else if (class(data) == "Seurat" && omics == "scATAC-seq"){
    
    omics_list <- list()
    counts <- data@assays[["ATAC"]]@counts 
    counts_matrix <- as.matrix(counts)
    omics_list[[omics]] <- counts_matrix
    return(omics_list)
    
  }
}


#' param_estimation
#' 
#' Evaluate the users parameters for single cell simulation and use SPARSim
#' to simulate the main dataset
#' @param omics named list containing the omics to simulate as names, which can 
#'    be "scRNA-seq" or "scATAC-seq, and the input count matrix
#' @param cellTypes list where the i-th element of the list contains the column 
#'    indices for i-th cell type. List must be a named list.
#' @return a named list with simulation parameters for each omics as values

param_estimation <- function(omics, cellTypes, numberCells = NULL, mean = NULL, sd = NULL){
  
  all_missing <- is.null(numberCells) && is.null(mean) && is.null(sd)
  all_specified <- !is.null(numberCells) && !is.null(mean) && !is.null(sd)
  
  if( !(all_missing || all_specified )){
    
    print("the user must either not provide the optional arguments or provide them all")
    return(NA)
    
  }
  
  # Use SPARSim to normalize the original dataset to use as seed
  N_omics <- length(omics)
  norm_list <- lapply(omics, scran_normalization)
  param_est_list <- list()
  
  for(i in 1:N_omics){
    # Estimate parameters for single cell estimation in each omic of interest
    # using SPARSim
    param_est <- SPARSim_estimate_parameter_from_data(raw_data = omics[[i]],
                                                      norm_data = norm_list[[i]],
                                                      conditions = cellTypes)
    param_est_list[[paste0("param_est_", names(omics)[i])]] <- param_est
    
  }
  
  if(all_missing){
    # If the user didnt input parameters, give back the empty list,
    # simulation will be run with default parameters by SPARSim
    return(param_est_list)
    
  } else if (all_specified){
    
    N_param_est_list<- length(param_est_list)
    N_cellTypes <- length(cellTypes)
    param_est_list_mod <- list()
    
    for(i in 1:N_param_est_list){
      cell_type_list <- list()
      
      for(j in 1:N_cellTypes){
        
        cond_param <- SPARSim_create_simulation_parameter(intensity = param_est_list[[i]][[j]][["intensity"]],
                                                          variability = param_est_list[[i]][[j]][["variability"]],
                                                          library_size = round(rnorm(n = numberCells[j], mean = mean[j], sd = sd[j])),
                                                          condition_name = param_est_list[[i]][[j]][["name"]],
                                                          feature_names = names(param_est_list[[i]][[j]][["intensity"]]))
        cell_type_list[[names(cellTypes)[j]]] <- cond_param 
        
      }
      param_est_list_mod[[paste0("param_est_", names(omics)[i])]] <- cell_type_list 
    }
    
    return(param_est_list_mod)
  } 
  
}

#' scMOSim
#' 
#' Performs multiomic simulation of single cell datasets
#' 
#' @param omics named list containing the omic to simulate as names, which can be "scRNA-seq" or "scATAC-seq, and the input count matrix as 
#' @param cellTypes list where the i-th element of the list contains the column indices for i-th experimental conditions. List must be a named list.
#' @param numberCells vector of numbers. The numbers correspond to the number of cells the user wants to simulate per each cell type. The length of the vector must be the same as length of \code{cellTypes}.
#' @param mean vector of numbers of mean per each cell type. Must be specified just if \code{numberCells} is specified.
#' @param sd vector of numbers of standard deviation per each cell type. Must be specified just if \code{numberCells} is specified.
#' @return a list of Seurat object, one per each omic. 
#' @export
#' 
#' @examples 
#' 
#' cellTypes <- list(cellA = c(1:20), cellB = c(161:191))
#' sim <- scMOSim(omicsList, cellTypes)
#' or
#' sim_with_arg <- scMOSim(omicsList, cellTypes, numberCells = c(10,20), 
#'       mean = c(2000000, 100000), sd = c(10^3, 10^3))

scMOSim <- function(omics, cellTypes, numberCells = NULL, mean = NULL, sd = NULL){
  # Estimate the parameters for simulation
  param_list <- param_estimation(omics, cellTypes, numberCells, mean, sd)
  
  N_param <- length(param_list)
  sim_list <- list()
  
  for(i in 1:N_param){
    # Use SPARSim to simulate each data type requested
    sim <- SPARSim_simulation(dataset_parameter = param_list[[i]])
    # Save into a list of simulated datasets
    sim <- sim[["count_matrix"]]
    sim_list[[paste0("sim_", names(omics)[i])]] <- sim
    
  }
  
  seu_obj <- list()
  N_sim <- length(sim_list)
  
  for(i in 1:N_sim){
    # Save the list of simulated datasets as a filled seurat object
    assay_name <- str_split(names(sim_list)[i], "-")[[1]][1]
    assay_name <- sub("sim_sc","",assay_name)
    seu <- CreateSeuratObject(counts = sim_list[[i]], assay = assay_name)
    seu_obj[[names(omics)[i]]] <- seu
    
  }
  
  return(seu_obj)
  
}