#' Data to showcase scRNA and scATAC-seq association
#' @format A dataframe with two columns and rows according to gene/feature relationships
#'  \describe{
#'      \item{Peak_ID}{ATAC chromosomic positions associated to genes}
#'      \item{Gene_ID}{RNA genes associated to peaks}
#'      }
#'      
#'  @source {Created in-house to serve as an example}
#' @name associationList
#' @docType data
#' @usage data("associationList")
data(associationList)

#' Data to test scMOSim
#' @format A seurat Object, subset from seuratData with RNA
#'  \describe{
#'      \item{assays}{RNA expression values}
#'      \item{meta.data}{annotations of celltypes}
#'      }
#'      
#'  @source {https://github.com/satijalab/seurat-data, we took 11 cells 
#'  from each of 4 celltypes}
#'  This is how:
#'  dat <- pbmcMultiome.SeuratData::pbmc.rna
#'  dat <- subset(x = dat, subset = seurat_annotations %in% c("CD4 TEM", 
#'                                                            "cDC", "Memory B", "Treg"))
#'  unique_cell_types <- unique(datATmeta.data$seurat_annotations)
#'  extracted_cells <- list()
#'  cellnames <- c()
#'  for (cell_type in unique_cell_types) {
#'    type_cells <- subset(dat, subset = seurat_annotations %in% cell_type)
#'    counts <- as.matrix(type_cellsATassays[["RNA"]]ATcounts)
#'    extracted_cells[[cell_type]] <- counts[, 1:10]
#'    cellnames <- append(cellnames, replicate(11, cell_type))
#'  }
#'  scrna <- Reduce(cbind, extracted_cells)
#' @name scrna
#' @docType data
#' @usage data("scrna")
data(scrna)

#' Data to test scMOSim
#' @format A seurat Object, subset from seuratData with ATAC
#'  \describe{
#'      \item{assays}{ATAC expression values}
#'      \item{meta.data}{annotations of celltypes}
#'      }
#'      
#'  @source {https://github.com/satijalab/seurat-data, we took 11 cells 
#'  from each of 4 celltypes}
#' @name scatac
#' @docType data
#' @usage data("scatac")
data(scatac)