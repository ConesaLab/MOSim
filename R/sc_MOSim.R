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
#' pbmc dataset from https://satijalab.org/seurat/articles/pbmc3k_tutorial.html
#'
#' @param omics_types A list of strings which can be either "scRNA-seq" or "scATAC-seq"
#' @param data A user input matrix with genes (peaks in case of scATAC-seq) as 
#'    rows and cells as columns. Alternatively scMOSim allows you to use a default 
#'    dataset (PBMC) by not specifying the argument.
#' @return a named list with omics type as name and the count matrix as value
#' @export
#'
#' @examples
#'
#' scOmics <- sc_omicData(list("scRNA-seq", "scATAC-seq"))
#' scRNAseq_user <- sc_omicData(omics_type = list("scRNA-seq"), data = list(count_matrix)) 
#'
sc_omicData <- function(omics_types, data = NULL){
  # Check for mandatory parameters
  if (missing(omics_types)){
    stop("You must provide the vector of omics to simulate.")
  }
  
  use_default_data <- function(){
    if (omics == "scRNA-seq"){
      rna_orig_counts <- sc_sampleData$rna_orig_counts
      return(list("scRNA-seq" = rna_orig_counts))
      
    } else if (omics =="scATAC-seq"){
      atac_orig_counts <- sc_sampleData$atac_orig_counts
      return(list("scATAC-seq" = atac_orig_counts))
      
    }
    return(NA)
  }
  
  use_provided_data <- function(){
    if (!is.list(data) || length(data) != 1 && length(data) != 2){
      print(paste0("the length of data is", length(data)))
      print("data must be a list of 1 or 2 elements")
      return(NA)
    } 
    
    omics_list <- list()
    N_data <- length(data)
    
    for (i in 1:N_data){
      
      if (! is.matrix(data[[i]]) && class(data[[i]]) != "Seurat"){
        print("Each element of data must be either matrix or a Seurat object")
        return(NA)
        
      } else if (is.matrix(data[[i]])){
        
        omics_list[[i]] <- data[[i]]
        
      } else if (class(data[[i]]) == "Seurat" && omics == "scRNA-seq"){
        
        counts <- data[[i]]@assays[["RNA"]]@counts
        counts_matrix <- as.matrix(counts)
        omics_list[[i]] <- counts_matrix
        
      } else if (class(data[[i]]) == "Seurat" && omics == "scATAC-seq"){
        
        counts <- data[[i]]@assays[["ATAC"]]@counts
        counts_matrix <- as.matrix(counts)
        omics_list[[i]] <- counts_matrix

      } else {
        print(paste0("Invalid class for data element ", i))
        return(NA)
      }
    }
    names(omics_list) <- omics_types[1:length(omics_list)]
    return(omics_list)
  }
  
  count_matrix_list<-list()
  
  for(omics in omics_types) {
    if (omics != "scRNA-seq" && omics != "scATAC-seq"){
      print("omics must be a either 'scRNA-seq' or 'scATAC-seq'")
      return(NA)
    }
    if (is.null(data)){
      count_matrix_list <- c(count_matrix_list, use_default_data())
    }
    else {
      count_matrix_list <- use_provided_data()
    }
  }
  return(count_matrix_list)
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
#' @param numberCells vector of numbers. The numbers correspond to the number 
#'    of cells the user wants to simulate per each cell type. The length of the 
#'       vector must be the same as length of \code{cellTypes}.
#' @param mean vector of numbers of mean per each cell type. Must be specified 
#'    just if \code{numberCells} is specified.
#' @param sd vector of numbers of standard deviation per each cell type. Must be 
#'    specified just if \code{numberCells} is specified.
#' @return a named list with simulation parameters for each omics as values.
#'
sc_param_estimation <- function(omics, cellTypes, numberCells = NULL, mean = NULL, sd = NULL){
  # Check for mandatory parameters
  if (missing(omics)){
    stop("You must provide the vector of omics to simulate.")
  }
  
  if (missing(cellTypes)) {
    stop("You must provide the correspondence of cells and celltypes")
  }
  
  all_missing <- is.null(numberCells) && is.null(mean) && is.null(sd)
  all_specified <- !is.null(numberCells) && !is.null(mean) && !is.null(sd)
 
  if( !(all_missing || all_specified )){
  
    print("the user must either not provide the optional arguments or provide them all")
    return(NA)
    
  }
  
  N_omics <- length(omics)
  norm_list <- lapply(omics, scran_normalization)
  param_est_list <- list()
  
  for(i in 1:N_omics){

    param_est <- SPARSim::SPARSim_estimate_parameter_from_data(raw_data = omics[[i]],
                                                    norm_data = norm_list[[i]],
                                                    conditions = cellTypes)
    param_est_list[[paste0("param_est_", names(omics)[i])]] <- param_est
    
  }
  
  N_param_est_list<- length(param_est_list)
  N_cellTypes <- length(cellTypes)
  param_est_list_mod <- list()
  
  if(all_missing){
    
    for (i in 1:N_param_est_list){
      cell_type_list <- list()
      
      for(j in 1:N_cellTypes){
        
        cond_param <- SPARSim::SPARSim_create_simulation_parameter(
                          intensity = param_est_list[[i]][[j]][["intensity"]],
                          variability = param_est_list[[i]][[j]][["variability"]],
                          library_size = param_est_list[[i]][[j]][["lib_size"]],
                          condition_name = param_est_list[[i]][[j]][["name"]],
                          feature_names = names(param_est_list[[i]][[j]][["intensity"]]))
        cell_type_list[[names(cellTypes)[j]]] <- cond_param
        
      }
      param_est_list_mod[[paste0("param_est_", names(omics)[i])]] <- cell_type_list 
    }
    
  } else if (all_specified){
    
      for(i in 1:N_param_est_list){
        cell_type_list <- list()
        
        for(j in 1:N_cellTypes){
        
          cond_param <- SPARSim::SPARSim_create_simulation_parameter(
                            intensity = param_est_list[[i]][[j]][["intensity"]],
                            variability = param_est_list[[i]][[j]][["variability"]],
                            library_size = round(rnorm(n = numberCells[j], mean = mean[j], sd = sd[j])),
                            condition_name = param_est_list[[i]][[j]][["name"]],
                            feature_names = names(param_est_list[[i]][[j]][["intensity"]]))
        
          cell_type_list[[names(cellTypes)[j]]] <- cond_param
        
        }
        param_est_list_mod[[paste0("param_est_", names(omics)[i])]] <- cell_type_list
      }
      
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
#' @param numberCells vector of numbers. The numbers correspond to the number of 
#'     cells the user wants to simulate per each cell type. The length of the 
#'         vector must be the same as length of \code{cellTypes}.
#' @param mean vector of numbers of mean per each cell type. Must be specified 
#'     just if \code{numberCells} is specified.The length of the vector must be 
#'         the same as length of \code{cellTypes}.
#' @param sd vector of numbers of standard deviation per each cell type. Must be 
#'     specified just if \code{numberCells} is specified.The length of the vector 
#'         must be the same as length of \code{cellTypes}.
#' @return a list of Seurat object, one per each omic.
#' @export
#'
#' @examples
#'
#' cellTypes <- list(cellA = c(1:20), cellB = c(161:191))
#' sim <- scMOSim(omicsList, cellTypes)
#' # or
#' sim_with_arg <- scMOSim(omicsList, cellTypes, numberCells = c(10,20),
#'       mean = c(2000000, 100000), sd = c(10^3, 10^3))
#'
scMOSim <- function(omics, cellTypes, numberCells = NULL, mean = NULL, sd = NULL){
    # Check for mandatory parameters
    if (missing(omics)){
      stop("You must provide the vector of omics to simulate.")
    }
  
    if (missing(cellTypes)){
      stop("You must provide the correspondence of cells and celltypes")
    }
    
    param_list <- sc_param_estimation(omics, cellTypes, numberCells, mean, sd)
    
    N_param <- length(param_list)
    sim_list <- list()
    
    for(i in 1:N_param){
      
      sim <- SPARSim::SPARSim_simulation(dataset_parameter = param_list[[i]])
      sim <- sim[["count_matrix"]]
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
    
    return(seu_obj)
    
}


#' sc_omicSim
#'
#' Defines regulatory functions of features for regulatory omics
#'
#' @param sim named list containing the omic simulated as names ("scRNA-seq" and 
#'     "scATAC-seq") , and seurat objects as values.
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
#'     cell_type), one per each couple of \code{cellTypes}.
#' @export
#'
#' @examples
#'
#'
#' omic_list <- sc_omicData(c("scRNA-seq","scATAC-seq"))
#' cell_types <- list(cellA = c(1:20), cellB = c(161:191))
#' sim <-scMOSim(omic_list, cell_types, numberCells = c(10,20), mean = c(2*10^6, 2*10^3), sd = c(10^3, 10^2))
#' cell_types <- list(cellA= c(1:10), cellB = c(11:30))
#' regulatorEffect = list('activator' = 0.8,'repressor' = 0.1,'NE' = 0.1)
#' sc_omicSim(sim, cell_types, totalFeatures = 500, regulatoreEffect = regulatorEffect)
#'
sc_omicSim <- function(sim, cellTypes, totalFeatures = NULL, 
                       regulatoreEffect = NULL, associationList = NULL ){
  # Check for mandatory parameters
  if (missing(sim)){
    stop("You must provide the named list of simulated omics.")
  }
    
  if (missing(cellTypes)){
    stop("You must provide the correspondence of cells and celltypes")
  }
  
  if(is.null(totalFeatures)){
    
    atac_counts <- sim[["scATAC-seq"]]@assays[["ATAC"]]@counts
    
  } else if (is.numeric(totalFeatures)){
    
    if(totalFeatures <= nrow(sim[["scATAC-seq"]]@assays[["ATAC"]]@counts)){
      
      # Set the number of features you want to select
      num_features <- totalFeatures
      # Get the data from the ATAC assay of the scATAC-seq object
      atac_data <- sim[["scATAC-seq"]]@assays[["ATAC"]]@counts
      # Subset the ATAC data to retain only the selected number of features
      atac_data <- atac_data[sample(nrow(atac_data), num_features, replace = FALSE), ]
      # Update the ATAC assay with the subsetted data
      sim[["scATAC-seq"]]@assays[["ATAC"]]@counts<- atac_data
      atac_counts <- sim[["scATAC-seq"]]@assays[["ATAC"]]@counts
      
    } else if (totalFeatures > nrow(sim[["scRNA-seq"]]@assays[["RNA"]]@counts)){
      
      print(paste("the number of totalFeatures you have inserted is higher than 
                  what's possible to be generated, ", 
                  nrow(sim[["scATAC-seq"]]@assays[["ATAC"]]@counts), 
                  " peaks were generated instead." ))
      atac_counts <- sim[["scATAC-seq"]]@assays[["ATAC"]]@counts
      
    }
  }
  
  rna_counts <- sim[["scRNA-seq"]]@assays[["RNA"]]@counts
  
  #calculate gene expression for each cellTypes
  gene_expression_list <- lapply(names(cellTypes), function(cell_type) {
    
    gene_expression <- rowSums(rna_counts[, cellTypes[[cell_type]]])
    names(gene_expression) <- rownames(rna_counts)
    return(gene_expression)
    
  })
  
  names(gene_expression_list) <- paste0("gene_expr_", names(cellTypes))
  
  if(length(cellTypes) > 2){
    print(paste0("the length of cellTypes is ", length(cellTypes)))
    
    da_peaks_atac_list <- lapply(seq_along(cellTypes), function(i) {
      
      j <- i %% length(cellTypes) + 1
      ident1 <- names(cellTypes)[i]
      ident2 <- names(cellTypes)[j]
      FindMarkers(object = sim[["scATAC-seq"]], ident.1 = ident1, ident.2 = ident2, 
                  min.pct = 0.05)
      
    })
    
    names(da_peaks_atac_list) <- paste0(names(cellTypes), "markers_", 
                                c(names(cellTypes)[-1], names(cellTypes)[1]))
    
  } else if(length(cellTypes) ==2 ){
    
    da_peaks_atac_list <- list()
    da_peaks_atac <- FindMarkers(object = sim[["scATAC-seq"]], 
                                 ident.1 = names(cellTypes[1]), 
                                 ident.2 = names(cellTypes[2]) , min.pct = 0.05)
    da_peaks_atac_list[[paste0("markers_", names(cellTypes[1]), 
                               "_", names(cellTypes[2]) )]] <- da_peaks_atac
    
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
    
  association_list <- sc_sampleData$seurat_association_list
  
  } else if (is.list(associationList)){
    
    if (length(associationList) == 2){
      
      association_list <- associationList
      
    } else {
      
      print("the associationList must be a dataframe having two columns reporting 
            peak ids and gene names")
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
        
        # Check if the peak is upregulated in cell type A and the gene expression 
        #of gene i for cell type A is >0
        if(peak %in% subset_list[[j]][[1]] && gene_expression_list[[j]][gene] > 0) {
          
          activity <- "activator"
          peak_df <- rbind(peak_df, data.frame(peak_id = peak, activity = activity, 
                            cell_type = strsplit(names(subset_list[[j]]),"_")[[1]][3], 
                            stringsAsFactors = FALSE))
          
          #check if the peak is upregulated in cell type A and the gene expression 
          #of gene i for cell type A is ==0 or
          # if the peak is downregulated in cell type A (up in B) and the gene 
          #expression of gene i for cell type A is >0
        } else if((peak %in% subset_list[[j]][[1]] && gene_expression_list[[j]][gene] == 0) ||
                  (peak %in% subset_list[[j]][[2]] && gene_expression_list[[j]][gene] > 0)) {
          
          activity <- "repressor"
          peak_df <- rbind(peak_df, data.frame(peak_id = peak, activity = activity, 
                            cell_type = strsplit(names(subset_list[[j]]),"_")[[1]][3], 
                            stringsAsFactors = FALSE))
          
          
        } 
        next_j <- (j %% length(subset_list)) +1
        #check if the peak is upregulated in cell type B and the gene expression 
        #of gene i for cell type B is >0
        if(peak %in% subset_list[[j]][[2]] && gene_expression_list[[next_j]][gene] > 0) {
          
          activity <- "activator"
          peak_df <- rbind(peak_df, data.frame(peak_id = peak, activity = activity, 
                            cell_type = strsplit(names(subset_list[[j]]),"_")[[2]][3], 
                            stringsAsFactors = FALSE))
          
          #check if the peak is upregulated in cell type B and the gene expression 
          #of gene i for cell type B is ==0 or
          #check if the peak is downregulated in cell type B (up in A) and the 
          #gene expression of gene i for cell type B is >0
        } else if((peak %in% subset_list[[j]][[2]] && gene_expression_list[[next_j]][gene] == 0) ||
                  (peak %in% subset_list[[j]][[1]] && gene_expression_list[[next_j]][gene] > 0)) {
          
          activity <- "repressor"
          peak_df <- rbind(peak_df, data.frame(peak_id = peak, activity = activity, 
                            cell_type = strsplit(names(subset_list[[j]]),"_")[[2]][3], 
                            stringsAsFactors = FALSE))
          
          #otherwise the regulator has a No effect Activity "NE"
        } 
        
        if (activity == "NE") {
          
          activity <- "NE"
          peak_df <- rbind(peak_df, data.frame(peak_id = peak, activity = activity, 
                            cell_type = "NA", stringsAsFactors = FALSE))
          
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
    
    
  } # if regulatorEffect is not null then it returns the sub-setted dataframe 
  #according to percentages passed
  
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
      
      # Apply the function to each dataframe in result_list and store the results 
      #in a new named list
      sub_result_list <- lapply(result_list, subset_df)
      names(sub_result_list) <- paste0(names(result_list), "_subset")
      return(sub_result_list)
      
    } else {
      
      print("regulatorEffect must be a named list with percentages as values of 
            length 3. The names must be 'activator', 'repressor', 'NE' ")
      return(NA)
      
    }
    
  }
  
}
