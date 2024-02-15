# Since SPARSim is neiher on CRAN nor bioconductor, we need to get the functions
# We need from them. These are
#    - SPARSIM_estimate_parameter_from_data
#    - SPARSIM_create_simulator_parameter
#    - SPARSIM_simulation
# SPARSim can be downloaded and installed from the sysbiobig GITLAB
# <https://gitlab.com/sysbiobig/sparsim>
# When using SPARSIM, please cite DOI: 10.1093/bioinformatics/btz752
# Baruzzo G, Patuzzi I, Di Camillo B. (2020) SPARSim single cell: a count data  
# simulator for scRNA-seq data


#' random_unif_interval
#' Function to call the C code
#'
#' @param size from sparsim
#' @param max_val from sparsim
#'
#' @export

random_unif_interval <- function(size, max_val) {
  .Call('_sparsim_random_unif_interval', PACKAGE = 'MOSim', size, max_val)
}


#' Create SPARSim simulation parameter
#'
#' Function to create a SPARSim simulation parameter.
#'
#' To simulate N feature (e.g. genes), user must specify N values of gene expression level and gene expression variability in the function input parameters \code{intensity} and \code{variability}, respectively.
#' To simulate M samples (i.e. cells), user must specify M values of sample library size in the function input parameter \code{library_size}.
#'
#' User can optionally specify the names to assign at the single feature and sample to simulate (function input parameters \code{feature_names} and \code{sample_names}, respectively,
#' as well as the name of the experimental condition (function input parameter \code{condition_name}). If the user does not specify such information, the function will set some default values.
#'
#' To simulate T different experimental conditions in a single count table, then T different simulation parameters must be created.
#'
#' @param intensity Array of gene expression intensity values
#' @param variability Array of gene expression variability values
#' @param library_size Array of library size values
#' @param feature_names Array of feature names. It must be of the same length of \code{intensity} array. If NA (default), feature will be automatically named "gene_1", "gene_2", ... "gene_<N>", where N = length(intensity)
#' @param sample_names Array of sample names. It must be of the same length of \code{library_size} array. If NA (defatul), sample will be automatically named "<condition_name>_cell1", "<condition_name>_cell2", ..., "<condition_name>_cell<M>", where M = length(library_size)
#' @param condition_name Name associated to the current experimental condition. If NA (default), it will be set to "cond<l1><l2>", where l1 and l2 are two random letters.
#' @param intensity_2 Array of gene expression intensity values for the second expression mode, if simulating genes with bimodal gene expression. Entries containing \code{NAs} will be ignored. If NULL (default), no bimodal gene expression is simulated.
#' @param variability_2 Array of gene expression variability values for the second expression mode, if simulating genes with bimodal gene expression. If NULL (default), no bimodal gene expression is simulated.
#' @param p_bimod Array of bimodal gene expression probabilities; the i-th value indicates the probability \code{p} of the i-th gene to be expressed in the first mode
#' (i.e. the one specified in the i-th entries of parameters \code{intensity} and \code{variability}); with probability \code{1-p} the i-th gene will be expressed in the second mode (i.e. the one specified in the i-th entries of parameters \code{intensity_2} and \code{variability_2})
#' @return SPARSim simulation parameter describing one experimental condition
#' @export
sparsim_create_simulation_parameter <- function(intensity, variability, library_size, feature_names = NA, sample_names = NA, condition_name = NA,
                                                intensity_2 = NULL, variability_2 = NULL, p_bimod = NULL){
  
  # assign feature name to "intensity" and "variability"
  if(length(feature_names)==1){ # potential NA
    if(is.na(feature_names)){
      feature_names <- paste0("gene_", c(1:length(intensity)) )
    }
  }
  names(intensity) <- feature_names
  names(variability) <- feature_names
  
  # assign condition name
  if(is.na(condition_name)){
    condition_name <- paste0("cond_", paste0(sample(LETTERS, size = 2), collapse = "") )
  }
  
  # assign sample names to "library_size"
  if(length(sample_names)==1){
    if(is.na(sample_names)){
      sample_names <- paste0(condition_name,"_cell",c(1:length(library_size)))
    }
  }
  names(library_size) <- sample_names
  
  cond_param <- list()
  cond_param$intensity <- intensity
  cond_param$variability <- variability
  cond_param$lib_size <- library_size
  cond_param$name <- condition_name
  
  # if bimodal gene expression is required, then set the related fields
  if(!is.null(intensity_2) & !is.null(variability_2) & !is.null(p_bimod)){
    names(intensity_2) <- feature_names
    names(variability_2) <- feature_names
    cond_param$intensity_2 <- intensity_2
    cond_param$variability_2 <- variability_2
    cond_param$p_bimod <- p_bimod
  }
  
  return(cond_param)
}


#' Estimate SPARSIm "intensity" parameter
#'
#' Function to estimate the intensity values from the genes in \code{data}. The intensity is computed as mean of normalized counts for each gene.
#'
#' This function is used in \code{sparsim_estimate_parameter_from_data} to compute SPARSim "intensity" parameter, given a real count table as input.
#' If the count table contains more than one experimental condition, then the function is applied to each experimental conditions.
#'
#' @param data normalized count data matrix (gene on rows, samples on columns). \code{rownames(data)} must contain gene names.
#' @return An array of intensity values having \code{N_genes} elements (\code{N_genes = nrow(data)}). Array entries are named with gene names.
#' @export
sparsim_estimate_intensity <- function(data){
  
  N_genes <- nrow(data)
  intensity_val <- array(0, dim = N_genes)
  names(intensity_val) <- rownames(data)
  
  # intesity as mean of normalized counts
  intensity_val <- apply(data, 1, mean)
  rm(data); gc()
  return (intensity_val)
  
}


#' Estimate SPARSim "variability" parameter
#'
#' Function to estimate the variability values from the genes in \code{data}.
#'
#' This function is used in \code{sparsim_estimate_parameter_from_data} to compute SPARSim "variability" parameter, given a real count table as input.
#' If the count table contains more than one experimental condition, then the function is applied to each experimental conditions.
#'
#' @param data raw count data matrix (gene on rows, samples on columns)
#' @return An array of variability values having \code{N_genes} elements (\code{N_genes = nrow(data)})
#' @export
sparsim_estimate_variability <- function (data){
  
  # compute percentage of entry >=0 for each gene
  perc_not_zeros_per_gene <- rowSums(data>0)/ncol(data)
  
  perc_not_zeros = 0
  ind_pass_filter <- (perc_not_zeros_per_gene > perc_not_zeros)
  N_gene_pass_filter <- sum(ind_pass_filter)
  
  variability <- array(NA, dim = nrow(data))
  
  # variability as edgeR dispersion of RAW counts in data
  if(N_gene_pass_filter > 1){
    
    #library(edgeR)
    f <- edgeR::DGEList(as.matrix(data[ind_pass_filter,]))
    gc()
    
    f <- edgeR::estimateGLMCommonDisp(f)
    gc()
    
    f <- edgeR::estimateGLMTrendedDisp(f)
    gc()
    
    if(any(is.na(f$trended.dispersion))){
      f$trended.dispersion <- NULL
      print("Trended dispresion is NA")
    }
    
    f <- edgeR::estimateGLMTagwiseDisp(f)
    gc()
    
    variability[ind_pass_filter] <- f$tagwise.dispersion
  }
  rm(data); gc()
  return(variability)
  
}



#' Estimate SPARSim "library size" parameter
#'
#' Function to estimate the library sizes from the samples in \code{data}.
#'
#' This function is used in \code{sparsim_estimate_parameter_from_data} to compute SPARSim "library size" parameter, given a real count table as input.
#' If the count table contains more than one experimental condition, then the function is applied to each experimental conditions.
#'
#' @param data raw count data matrix (gene on rows, samples on columns)
#' @return An array of library size values having \code{N_samples} elements (\code{N_samples = ncol(data)})
#' @export
sparsim_estimate_library_size <- function (data){
  return(colSums(data))
}



#' Estimate SPARSim simulation parameter from a given count table
#'
#' Function to estimate SPARSim simulation parameters (intensity, variability and library sizes) from a real count table.
#' If the real count table contains more than one experimental condition, it is possible to estimate the parameters for each experimental condition.
#'
#' @param raw_data count matrix (gene on rows, samples on columns) containing raw count data
#' @param norm_data count matrix (gene on rows, samples on columns) containing normalized count data
#' @param conditions list where the i-th element of the list contains the column indices for i-th experimental conditions. List must be a named list.
#' @return A SPARSim simulation parameters
#' @export
sparsim_estimate_parameter_from_data <- function (raw_data, norm_data, conditions){
  
  # Number of experimental conditions
  N_cond <- length(conditions)
  
  # Number of genes
  N_feature <- nrow(raw_data)
  
  
  dataset_parameter <- list()
  
  cond_i <- 1
  for (cond in conditions){# for each experimental condition
    
    #print(paste0("Experimental condition ", cond_i))
    
    gene_names <- rownames(raw_data)
    sample_names <- colnames(raw_data)
    
    #print("...estimating gene intensity")
    # estimate intesity

    N_genes <- length(gene_names)
    estim_intensity <- array(0, dim = N_genes)
    names(estim_intensity) <- gene_names
    for(i in 1:N_genes){
      estim_intensity[i] <- mean(norm_data[i,cond])
    }

    gc()
    
    #print("...estimating gene variability")
    # estimate variability
    estim_variability <- sparsim_estimate_variability (data = raw_data[,cond])
    
    #print("...estimating library size")
    # estimate library size
    estim_lib_size <- sparsim_estimate_library_size(data = raw_data[,cond])
    
    # get feature names, if present
    estim_feature_names <- NA
    if(!is.null(gene_names)){
      estim_feature_names <- gene_names
    }
    
    # get sample names, if present
    estim_sample_names <- NA
    if(!is.null(sample_names)){
      estim_sample_names <- sample_names[cond]
    }
    
    # get condition name, if present
    estim_cond_name <- paste0("cond_",cond_i)
    if(!is.null(names(conditions))){
      estim_cond_name <- names(conditions)[cond_i]
    }
    
    #print("...creating SPARSim simulation parameter")
    cond_parameter <- sparsim_create_simulation_parameter (intensity = estim_intensity,
                                                           variability = estim_variability,
                                                           library_size = estim_lib_size,
                                                           feature_names = estim_feature_names,
                                                           sample_names = estim_sample_names,
                                                           condition_name = estim_cond_name)
    
    dataset_parameter[[cond_i]] <- cond_parameter
    cond_i <- cond_i + 1
    
  }
  
  return(dataset_parameter)
}

#' Simulate technical variability
#'
#' Function to simulate the technical variability (i.e. a multivariate hypergeometric on a gamma expression value array)
#'
#' @param avgAbund array containing the intensity values for each feature. It describes the intensity of a single sample
#' @param seqdepth sequencing depth (i.e. sample size of the MH)
#' @param digits number of digits for random number generation
#' @param max_val max value for random number generation
#' @return An array of \code{length(avgAbund)} elements representing the count values for the current sample
simulate_hyper <- function(avgAbund, seqdepth = NULL, digits, max_val) {
  
  if(seqdepth > 10^digits){
    print("seqdepth must be <= than 10^digits")
    return(NA)
  }
  
  index <- MOSim::random_unif_interval(round(seqdepth*1.1),max_val)
  
  while_i <- 1
  index <- unique(index)
  while (length(index)<seqdepth){
    new_ind <- MOSim::random_unif_interval( (seqdepth-length(index))*(while_i*10),max_val)
    index<-c(index,new_ind)
    while_i <- while_i +1
    index<-unique(index)
  }
  index <- index[1:seqdepth]
  
  cumulative_orig<-cumsum(avgAbund)
  check_orig<-findInterval(index,cumulative_orig,left.open = T)+1
  tmp_tabulate <- tabulate(check_orig, nbins = length(avgAbund))
  names(tmp_tabulate)<-names(avgAbund)
  
  return(tmp_tabulate)
}

#' Function to simulate a raw count table
#'
#' @param dataset_parameter list containing, the intensity, variability and lib sizes of each experimental condition. It is the return value of "estimate_parameter_from_data" or could be created by the users
#' @param output_sim_param_matrices boolean flag. If TRUE, the function will output two additional matrices, called abundance_matrix and variability_matrix, containing the gene intensities and gene variabilities used as simulation input. (Default: FALSE)
#' @param output_batch_matrix boolean flag. If TRUE, the function will output an additional matrix, called batch_factors_matrix, containing the multiplicative factors used in batch effect simulation. (Default: FALSE)
#' @param count_data_simulation_seed inherited from sparsim
#' @return A list of 5 elements:
#'
#'   - \code{count_matrix}: the simulated count matrix (genes on rows, samples on columns)
#'
#'   - \code{gene_matrix}: the simulated gene expression levels (genes on rows, samples on columns)
#'
#'   - \code{abundance_matrix}: the input gene intensity values provided as input (genes on rows, samples on columns), if \code{output_sim_param_matrices} = TRUE. NULL otherwise.
#'
#'   - \code{variability_matrix}: the input gene variability values provided as input (genes on rows, samples on columns), if \code{output_sim_param_matrices} = TRUE. NULL otherwise.
#'
#'   - \code{batch_factors_matrix}: the multiplicative factor used in batch generation (genes on rows, samples on columns), if \code{output_batch_matrix} = TRUE. NULL otherwise.
#'
#' @export
sparsim_simulation <- function(dataset_parameter,
                               output_sim_param_matrices = FALSE,
                               output_batch_matrix = FALSE,
                               count_data_simulation_seed = NULL){
  
  
  
  # number of experimental condition
  N_cond <- length(dataset_parameter)
  
  # number of genes
  N_genes <- length(dataset_parameter[[1]]$intensity)
  
  # total number of cells
  N_cell <- sum(unlist(lapply( dataset_parameter, function(x){return(length(x$lib_size))})))
  
  cat("Number of experimental conditions: ", N_cond, "\n")
  cat("Number of genes: ", N_genes, "\n")
  cat("Number of cells: ", N_cell, "\n")
  
  
  cat("Setting gene expression intensity... ", "\n")
  # initialize gene expression level matrix
  gene_expression_matrix <- matrix(0, nrow = N_genes, ncol = N_cell)
  rownames(gene_expression_matrix) <- names(dataset_parameter[[1]]$intensity)
  column_index <- 1
  new_col_names <- ""
  for(cond in 1:N_cond){
    
    cell_cond_name <- names(dataset_parameter[[cond]]$lib_size)
    new_col_names <- c( new_col_names ,  cell_cond_name)
    
  }
  new_col_names <- new_col_names[-1]
  colnames(gene_expression_matrix) <- new_col_names
  
  
  # fill gene expression level matrix with initial (still without biological variability) gene expression values
  for(cond in 1:N_cond){
    tmp_cell_ids <- names(dataset_parameter[[cond]]$lib_size)
    tmp_N_cell <- length(tmp_cell_ids)
    
    N_cell_mem <- 5000
    if(tmp_N_cell < N_cell_mem){
      gene_expression_matrix[, tmp_cell_ids] <- rep(dataset_parameter[[cond]]$intensity, tmp_N_cell)
    }else{
      tmp_data <- rep(dataset_parameter[[cond]]$intensity, N_cell_mem)
      N_mult <- floor(tmp_N_cell/N_cell_mem)
      
      for(ind in c(1:N_mult)){
        ind_range <- c( ((ind-1)*N_cell_mem+1) : (ind*N_cell_mem) )
        gene_expression_matrix[, tmp_cell_ids[ind_range]  ] <- tmp_data
      }
      if(N_mult*N_cell_mem < tmp_N_cell){
        ind_range <- c( (N_mult*N_cell_mem+1) : tmp_N_cell)
        gene_expression_matrix[, tmp_cell_ids[ind_range]  ] <- tmp_data[c(1: (length(ind_range)*N_genes) )]
      }
      
      rm(tmp_data); gc(verbose = FALSE)
    }
    
  }
  
  
  cat("Setting gene expression variability ... ", "\n")
  # initialize gene expression variability matrix
  gene_expression_var_matrix <- matrix(0, nrow = N_genes, ncol = N_cell)
  rownames(gene_expression_var_matrix) <- rownames(gene_expression_matrix)
  colnames(gene_expression_var_matrix) <- colnames(gene_expression_matrix)
  
  # fill gene expression variability matrix
  for(cond in 1:N_cond){
    tmp_cell_ids <- names(dataset_parameter[[cond]]$lib_size)
    tmp_N_cell <- length(tmp_cell_ids)
    
    # get variability values and fix possible NA or negative values
    variability_values <- dataset_parameter[[cond]]$variability
    variability_values[is.na(variability_values)]<-0
    variability_values[variability_values<0]<-min(variability_values[variability_values>0])
    
    N_cell_mem <- 5000
    if(tmp_N_cell < N_cell_mem){
      gene_expression_var_matrix[, tmp_cell_ids] <- rep(variability_values, tmp_N_cell)
    }else{
      tmp_data <- rep(variability_values, N_cell_mem)
      N_mult <- floor(tmp_N_cell/N_cell_mem)
      
      for(ind in c(1:N_mult)){
        ind_range <- c( ((ind-1)*N_cell_mem+1) : (ind*N_cell_mem) )
        gene_expression_var_matrix[, tmp_cell_ids[ind_range]  ] <- tmp_data
      }
      if(N_mult*N_cell_mem < tmp_N_cell){
        ind_range <- c( (N_mult*N_cell_mem+1) : tmp_N_cell)
        gene_expression_var_matrix[, tmp_cell_ids[ind_range]  ] <- tmp_data[c(1: (length(ind_range)*N_genes) )]
      }
      
      rm(tmp_data); gc(verbose = FALSE)
    }
    
  }
  

  for(cond in 1:N_cond){
    
    if(!is.null(dataset_parameter[[cond]]$intensity_2)){ # check if bimodal gene expression is present in this experimental condition
      tmp_cell_ids <- names(dataset_parameter[[cond]]$lib_size)
      tmp_N_cell <- length(tmp_cell_ids)
      
      # get intensity values for mode 2
      intensity_values <- dataset_parameter[[cond]]$intensity_2
      
      # get variability values for mode 2 and fix possible NA or negative values
      variability_values <- dataset_parameter[[cond]]$variability_2
      variability_values[is.na(variability_values)]<-0
      variability_values[variability_values<0]<-min(variability_values[variability_values>0])
      
      # get probability for mode 2 (as 1 - p_mode)
      p_mode_2 <- 1 - dataset_parameter[[cond]]$p_bimod
      
      # identify bimodal genes (i.e. having intensities != NA and probability of mode 2 beign greater than 0)
      bimodal_genes_id <- ( ( !is.na(intensity_values) ) & (p_mode_2>0) )
      bimodal_genes_name <- rownames(gene_expression_matrix)[bimodal_genes_id]
      
      intensity_values_tmp <- intensity_values[bimodal_genes_id]
      variability_values_tmp <- variability_values[bimodal_genes_id]
      p_mode_2_tmp <- p_mode_2[bimodal_genes_id]
      
      N_bimod_gene <- length(intensity_values_tmp)
      
      # simulate bimodal genes in the cells
      for(cell in tmp_cell_ids){
        
        # use the probability to be in mode 2 to select which genes are in the second expression mode in the current cell
        mode2_ind <- as.logical(rbinom(N_bimod_gene, 1, p_mode_2_tmp))
        
        # only for the selected genes, change the gene expression intensity matrix and the gene expression variability matrix
        # with the values describing the second expression mode
        gene_expression_matrix[bimodal_genes_name[mode2_ind],cell] <- intensity_values_tmp[mode2_ind]
        gene_expression_var_matrix[bimodal_genes_name[mode2_ind],cell] <- variability_values_tmp[mode2_ind]
      }
      
    }
    
  }
  
  
  
  # get samples library sizes
  sample_lib_size <- unlist(lapply( dataset_parameter, function(x){return(x$lib_size)} ), use.names = FALSE)
  names(sample_lib_size) <- colnames(gene_expression_matrix)
  
  
  batch_factor_matrix <- NULL
  
  # if user does not require the batch factor matrix as output, free the memory and return a NULL matrix
  if(output_batch_matrix == FALSE){
    rm(batch_factor_matrix); gc(verbose = FALSE);
    batch_factor_matrix <- NULL
  }
  
  
  
  ### Simulate biological variability
  cat("Simulating biological variability ... ", "\n")
  
  gene_expression_matrix_bio_var <- gene_expression_matrix
  
  # simulate biological variability using a gamma
  gene_expression_matrix_bio_var <- matrix(
    stats::rgamma(n = nrow(gene_expression_matrix_bio_var)*ncol(gene_expression_matrix_bio_var),
           shape = 1/gene_expression_var_matrix,
           scale = gene_expression_var_matrix*gene_expression_matrix),
    ncol = ncol(gene_expression_matrix)
  )
  
  # feature having null (i.e. zero) variability have no variability, so "undo" the gamma
  zero_var_index <- (gene_expression_var_matrix == 0)
  gene_expression_matrix_bio_var [zero_var_index] <- gene_expression_matrix [zero_var_index]
  
  rownames(gene_expression_matrix_bio_var) <- rownames(gene_expression_matrix)
  colnames(gene_expression_matrix_bio_var) <- colnames(gene_expression_matrix)
  
  # if user does not require the gene intensity and gene variability matrices as output, free the memory and return NULL matrices
  if(output_sim_param_matrices == FALSE){
    rm(gene_expression_matrix, gene_expression_var_matrix); gc(verbose = FALSE)
    gene_expression_matrix <- NULL
    gene_expression_var_matrix <- NULL
  }
  
  # identify the maximum library size
  max_lib_size <- max(sample_lib_size)
  
  # set the input fragment library size (the population size parameter for MH)
  input_fragment_lib_size <- max_lib_size * 10^2
  digits<-floor(log10(input_fragment_lib_size))+1
  new_libsize<-10^digits
  
  gene_expression_matrix_bio_var_scaled <- round(t(t(gene_expression_matrix_bio_var)*new_libsize/colSums(gene_expression_matrix_bio_var)))
  num_fragment <- colSums(gene_expression_matrix_bio_var_scaled)
  
  
  ### Simulate technical variability
  cat("Simulating technical variability ... ", "\n")
  
  sim_count_matrix <- matrix(0, ncol = ncol(gene_expression_matrix_bio_var_scaled), nrow = nrow(gene_expression_matrix_bio_var_scaled) )
  rownames(sim_count_matrix) <- rownames(gene_expression_matrix_bio_var_scaled)
  colnames(sim_count_matrix) <- colnames(gene_expression_matrix_bio_var_scaled)
  
  # if the user set NO seed for count data simulation
  if(is.null(count_data_simulation_seed)){
    for(sample in colnames(sim_count_matrix)){
      sim_count_matrix[, sample] <- simulate_hyper(avgAbund = gene_expression_matrix_bio_var_scaled[, sample],
                                                   seqdepth = sample_lib_size[sample],
                                                   digits = digits,
                                                   max_val = num_fragment[sample])
    }
  }

  
  ### collect and return the results
  return(list(count_matrix = sim_count_matrix,
              gene_matrix = gene_expression_matrix_bio_var,
              abundance_matrix = gene_expression_matrix,
              variability_matrix = gene_expression_var_matrix,
              batch_factors_matrix = batch_factor_matrix )
  )
}




