## Test scMOSim

suppressPackageStartupMessages(library(testthat))

# sc_omicData function
testthat::test_that("Passing a wrong string in 'omics' returns Error", {
  testthat::expect_error(MOSim::sc_omicData(c("scR-seq"), rna_orig_counts))
})

testthat::test_that("Passing an object in data which is neither a matrix or Seurat obj returns error", {
  vector <- c(1,2,3)
  testthat::expect_error(MOSim::sc_omicData(c("scATAC-seq"), vector))
})

testthat::test_that("Passing 'scRNA-seq' as omic and a Seurat obj as data returns a list", {
  scRNA <- MOSim::sc_omicData("scRNA-seq")
  count <- scRNA[["scRNA-seq"]]
  Seurat_obj <- Seurat::CreateSeuratObject(counts = count, assay = 'RNA')
  testthat::expect_type(sc_omicData(c("scRNA-seq"), c(Seurat_obj)), "list")
})

testthat::test_that("Passing an array ('scRNA-seq','scATAC-seq') as omic returns a list of length 2", {
  res<-MOSim::sc_omicData(list("scRNA-seq","scATAC-seq"))
  testthat::expect_equal(length(res), 2)
})

#scMOSim function
testthat::test_that("param_estimation returns a list", {
  omic_list <- MOSim::sc_omicData(c("scRNA-seq", "scATAC-seq"))
  conditions <- list('CD4_TEM' = c(1:60), 'cDC' = c(299:310))
  testthat::expect_type(MOSim::scMOSim(omic_list, conditions, numberCells = c(10,20), mean = c(2*10^6, 2*10^3), sd = c(10^3, 10^2)),"list")
})

testthat::test_that("Not passing all optional arguments at once returns NA", {
  omic_list <- MOSim::sc_omicData(c("scRNA-seq","scATAC-seq"))
  conditions <- list('CD4_TEM' = c(1:60), 'cDC' = c(299:310))
  testthat::expect_error(MOSim::scMOSim(omic_list, conditions, numberCells = c(10,20), sd = c(10^3, 10^2)))
})


testthat::test_that("scMOSim returns a list with S4 obj as values", {
  omic_list <- MOSim::sc_omicData(c("scRNA-seq","scATAC-seq"))
  cell_types <- list('CD4_TEM' = c(1:60), 'cDC' = c(299:310))
  sim <-MOSim::scMOSim(omic_list, cell_types, numberCells = c(10,20), mean = c(2*10^6, 2*10^3), sd = c(10^3, 10^2))
  testthat::expect_type(sim[[1]][[1]][[1]], "S4")
})

testthat::test_that("checking that scMOSim is able to simulate groups and replicates", {
  omicsList <- MOSim::sc_omicData(list("scRNA-seq", "scATAC-seq"))
  cell_types <- list('CD4_TEM' = c(1:60), 'cDC' = c(299:310), 'Memory_B' = c(497:510), 'Treg' = c(868:900))
  testing_groupsreps <- MOSim::scMOSim(omicsList, cell_types, numberReps = 2, numberGroups = 3, 
                                diffGenes = list(c(0.5, 0.4), c(0.3, 0.3)), minFC = 0.25, maxFC = 4,
                                numberCells = NULL, mean = NULL, sd = NULL, regulatorEffect = list(c(0.2, 0.5), c(0.5, 0.2)))
  testthat::expect_type(testing_groupsreps, "list")
})

testthat::test_that("If numberGroups > 1, length(diffGenes) must be == (numberGroups -1)", {
  omicsList <- MOSim::sc_omicData(list("scRNA-seq", "scATAC-seq"))
  cell_types <- list('CD4_TEM' = c(1:60), 'cDC' = c(299:310), 'Memory_B' = c(497:510), 'Treg' = c(868:900))
  testthat::expect_error(MOSim::scMOSim(omicsList, cell_types, numberReps = 2, numberGroups = 2, 
                       diffGenes = c(0.2, 0.2), minFC = 0.25, maxFC = 4,
                       numberCells = NULL, mean = NULL, sd = NULL))
})

#make_cluster_patterns function
testthat::test_that("make_cluster_patterns returns a tibble of dimensions 
          numcells^2xnumcells",{
            patterns <- MOSim::make_cluster_patterns(4, 4)
            expected <- c(4, 4)
            testthat::expect_equal(dim(patterns), expected)
          })

## simulate_coexpression function
testthat::test_that("The number of patterns we want are generated", {
  omic_list <- MOSim::sc_omicData(c("scRNA-seq"))
  cell_types <- list('CD4_TEM' = c(1:60), 'cDC' = c(299:310), 
                     'Memory_B' = c(497:520), 'Treg' = c(868:900))
  sim <-MOSim::scMOSim(omic_list, cell_types, numberCells = c(100, 100, 100, 100))
  
  sim_matrix <- sim$Group_1$Rep_1$`sim_scRNA-seq`@assays$RNA@counts
  numberCells <- c(100, 100, 100, 100)

  patterns <- MOSim::make_cluster_patterns(length(cell_types), clusters = 8)
  coexpr_results <- MOSim::simulate_coexpression(sim_matrix, numberCells,
                                          feature_no = 8000, 
                                          patterns)
  testthat::expect_type(coexpr_results, "list")
})

## Test the patterns
testthat::test_that("when there are two opposite patterns, it gives them back", {
  patterns <- tibble(one = c(TRUE, FALSE, TRUE, FALSE), 
                     two = c(TRUE, FALSE, TRUE, TRUE), 
                     three = c(FALSE, TRUE, FALSE, TRUE), 
                     four = c(FALSE, TRUE, TRUE, TRUE))
  opposite_indices <- check_patterns(patterns)
  testthat::expect_equal(length(opposite_indices), 2)
})


#sc_omicSim function
test_that("sc_omicSim returns a list", {
  omic_list <- sc_omicData(c("scRNA-seq","scATAC-seq"))
  cell_types <- list('CD4_TEM' = c(1:60), 'cDC' = c(299:310))
  sim <-scMOSim(omic_list, cell_types, numberCells = c(10,20), mean = c(2*10^6, 2*10^3), sd = c(10^3, 10^2))
  cell_types <- list('CD4_TEM' = c(1:60), 'cDC' = c(299:310))
  integration <- sc_omicSim(sim, cell_types, totalFeatures = 500)
  expect_type(integration, "list")
})

test_that("sc_omicSim returns the expected subarray for the activity column", {
  expected <- c("activator","repressor","NE","NE","NE","activator","repressor","activator","repressor","activator")
  omic_list <- sc_omicData(c("scRNA-seq","scATAC-seq"))
  cell_types <- list('CD4_TEM' = c(1:60), 'cDC' = c(299:310))
  sim <-scMOSim(omic_list, cell_types, numberCells = c(10,20), mean = c(2*10^6, 2*10^3), sd = c(10^3, 10^2))
  cell_types <- list('CD4_TEM' = c(1:60), 'cDC' = c(299:310))
  integration <- sc_omicSim(sim, cell_types, totalFeatures = 500)
  expect_equal(integration[["markers_cellA_cellB"]]$activity[1:10], expected)
})


test_that("sc_omicSim also works for a simulation with 4 cell types",{
  rna_counts <- readRDS("./data/RNA_4CellTypes.rds")
  atac_counts <- readRDS("./data/ATAC_4CellTypes.rds")
  omicData_4cell <- sc_omicData(omics_types =  list("scRNA-seq","scATAC-seq"))
  cell_types <- list('CD4_TEM' = c(1:60), 'cDC' = c(299:310), 'Memory_B' = c(497:510), 'Treg' = c(868:900))
  sim_4cell <- scMOSim(omicData_4cell, cell_types, numberCells = c(10,20,10,20), mean = c(2*10^6, 2*10^3, 2*10^5, 2*10^6), sd = c(10^3, 10^2, 10^2, 10^3))
  cell_types_integration <- list(cellA = c(1:10), cellB = c(11:30), cellC =c(31:40), cellD=c(41:60))
  expect_type(sc_omicSim(sim_4cell, cell_types_integration, totalFeatures = 500), "list")
})


test_that("checking if sc_omicSim does not return two identical list for markers_cellA_cellB and markers_cellB_cellC",{
  rna_counts <- readRDS("./data/RNA_4CellTypes.rds")
  atac_counts <- readRDS("./data/ATAC_4CellTypes.rds")
  omicData_4cell <- sc_omicData(omics_types =  list("scRNA-seq","scATAC-seq"), data = list(rna_counts,atac_counts))
  cell_types <- list('CD4_TEM' = c(1:60), 'cDC' = c(299:310), 'Memory_B' = c(497:510), 'Treg' = c(868:900))
  sim_4cell <- scMOSim(omicData_4cell, cell_types, numberCells = c(10,20,10,20), mean = c(2*10^6, 2*10^3, 2*10^5, 2*10^6), sd = c(10^3, 10^2, 10^2, 10^3))
  cell_types_integration <- list(cellA = c(1:10), cellB = c(11:30), cellC =c(31:40), cellD=c(41:60))
  omicSim <- sc_omicSim(sim_4cell, cell_types_integration, totalFeatures = 500)
  expect_true(!identical(cell_types_integration[[1]][[1]], cell_types_integration[[2]][[1]]))
})


test_that("checking that passing the argument regulatorEffect to sc_omicSim subsets the activity column in the correct way",{
  rna_counts <- readRDS("./data/RNA_4CellTypes.rds")
  atac_counts <- readRDS("./data/ATAC_4CellTypes.rds")
  omicData_4cell <- sc_omicData(omics_types =  list("scRNA-seq","scATAC-seq"), data = list(rna_counts,atac_counts))
  cell_types <- list('CD4_TEM' = c(1:60), 'cDC' = c(299:310), 'Memory_B' = c(497:510), 'Treg' = c(868:900))
  sim_4cell <- scMOSim(omicData_4cell, cell_types, numberCells = c(10,20,10,20), mean = c(2*10^6, 2*10^3, 2*10^5, 2*10^6), sd = c(10^3, 10^2, 10^2, 10^3))
  
  cell_types_integration <- list('CD4_TEM' = c(1:60), 'cDC' = c(299:310), 'Memory_B' = c(497:510), 'Treg' = c(868:900))
  regulatorEffect <- list('activator' = 0.8,'repressor' = 0.1,'NE' = 0.1)
  omicSim <- sc_omicSim(sim_4cell, cell_types_integration, totalFeatures = 500)
  omicSim_subset <- sc_omicSim(sim_4cell, cell_types_integration, totalFeatures = 500, regulatorEffect = regulatorEffect)
  NE <- sum(omicSim[["markers_cellA_cellB"]][["activity"]]=="NE")
  NE_rounded <- round(NE*0.1)
  NE_subset <- sum(omicSim_subset[["markers_cellA_cellB_subset"]][["activity"]] =="NE")
  expect_true(NE_rounded == NE_subset)
})
