#' @param omic A string which can be either "scRNA-seq" or "scATAC-seq"
#' @param data A user input matrix with genes (peaks in case of scATAC-seq) as rows and cells as columns. Alternatively MOSim allows user to estimate the input parameters from an existing count table by typing 'example_matrix'
#' @return a named list with omic type as name and the count matrix as value

length(rna_orig_counts)

sc_omicData <- function(omic, data = NULL){

  if (omic != "scRNA-seq" && omic != "scATAC-seq"){

    print("omic must be a either 'scRNA-seq' or 'scATAC-seq'")
    return(NA)
  }


  if (is.null(data)){ 
      Seurat_obj <- readRDS("pbmc_idents_SeuratObject.rds")

      #Selecting pDC and Bmemory cell from the larger dataset
      Seurat_Bmemory<- subset(x = Seurat_obj, idents = c("B memory"))
      Seurat_pDC <- subset(x = Seurat_obj, idents = c("pDC"))
      example_matrix <- merge(x= Seurat_Bmemory, y= Seurat_pDC)

      if (omic == "scRNA-seq"){ 

        ##scRNA##
        rna_obj <- example_matrix@assays[["RNA"]]
        dgCmatrix_counts <- rna_obj@counts
        rna_orig_counts <- as.matrix(dgCmatrix_counts)

        omic_list <- list("scRNA-seq" = rna_orig_counts)
        return(omic_list)

      } else if (omic =="scATAC-seq"){

        ##scATAC##
        atac_obj <- example_matrix@assays[["ATAC"]]
        dgCmatrix_counts <- atac_obj@counts
        atac_orig_counts <- as.matrix(dgCmatrix_counts)

         omic_list <- list("scATAC-seq" = atac_orig_counts)
        return(omic_list)

      }
  } 

  if (! is.matrix(data)){

    print("data must be a matrix")
    return(NA)

  } else if (is.matrix(data)){
    print(omic)
    omic_list <- list()
    omic_list[[omic]] <- data  #rotto qui
    return(omic_list)
  }
}

prova1 <- sc_omicData("scR-seq", rna_orig_counts)
prova2 <- sc_omicData("scATAC-seq")
prova3 <- sc_omicData("scATAC-seq", c)
ca <- 2
pa <- 4
!is.matrix(c)
rna_orig_counts <- as.matrix(dgCmatrix_counts)

omic <- "scRNA-seq"
provona <- tibble::lst(omic, rna_orig_counts)


#@param omics character vector containing the names of the omics to simulate,
#       which can be "scRNA-seq" or "scATAC-seq



sc_MOSim <- function(omics, numberCellTypes, numberCells = FALSE, mean = FALSE, sd = FALSE, sim_parameter = FALSE ){
}
