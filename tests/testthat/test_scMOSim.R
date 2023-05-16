## Test scMOSim

library(testthat)

# sc_omicData function
test_that("Passing a wrong string in 'omics' returns NA", {
  expect_message(sc_omicData(c("scR-seq"), rna_orig_counts), NA)
})

test_that("Passing an object in data which is neither a matrix or Seurat obj returns NA", {
  vector <- c(1,2,3)
  expect_message(sc_omicData(c("scATAC-seq"), vector), NA)
})

test_that("Passing 'scATAC-seq' as omic returns the expected subarray", {
  res<-sc_omicData("scATAC-seq")
  expected<-c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0,
              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  expect_equal(res$`scATAC-seq`[1:50], expected)
})

test_that("Passing 'scRNA-seq' as omic and a Seurat obj as data returns a list", {
  scRNA <- sc_omicData("scRNA-seq")
  count <- scRNA[["scRNA-seq"]]
  Seurat_obj <- Seurat::CreateSeuratObject(counts = count, assay = 'RNA')
  expect_type(sc_omicData(c("scRNA-seq"), c(Seurat_obj)), "list")
})

test_that("Passing an array ('scRNA-seq','scATAC-seq') as omic returns a list of length 2", {
  res<-sc_omicData(c("scRNA-seq","scATAC-seq"))
  expect_equal(length(res), 2)
})

#sc_param_estimation function
test_that("param_estimation returns a list", {
  omic_list <- sc_omicData(c("scRNA-seq","scATAC-seq"))
  conditions <- list(cellA = c(1:20), cellB = c(161:191))
  expect_type(sc_param_estimation(omic_list, conditions, numberCells = c(10,20), mean = c(2*10^6, 2*10^3), sd = c(10^3, 10^2)),"list")
})

test_that("Not passing all optional arguments at once returns NA", {
  omic_list <- sc_omicData(c("scRNA-seq","scATAC-seq"))
  conditions <- list(cellA = c(1:20), cellB = c(161:191))
  expect_message(sc_param_estimation(omic_list, conditions, numberCells = c(10,20),sd = c(10^3, 10^2)), NA)
})


#scMOSim function
test_that("scMOSim returns a list with S4 obj as values", {
  omic_list <- sc_omicData(c("scRNA-seq","scATAC-seq"))
  cell_types <- list(cellA = c(1:20), cellB = c(161:191))
  sim <-scMOSim(omic_list, cell_types, numberCells = c(10,20), mean = c(2*10^6, 2*10^3), sd = c(10^3, 10^2))
  expect_type(sim[[1]], "S4")
  
})

#sc_omicSim function
test_that("sc_omicSim returns a list", {
  omic_list <- sc_omicData(c("scRNA-seq","scATAC-seq"))
  cell_types <- list(cellA = c(1:20), cellB = c(161:191))
  sim <-scMOSim(omic_list, cell_types, numberCells = c(10,20), mean = c(2*10^6, 2*10^3), sd = c(10^3, 10^2))
  cell_types <- list(cellA= c(1:10), cellB = c(11:30))
  integration <- sc_omicSim(sim, cell_types, totalFeatures = 500)
  expect_type(integration, "list")
})

test_that("sc_omicSim returns the expected subarray for the activity column", {
  expected <- c("activator","repressor","NE","NE","NE","activator","repressor","activator","repressor","activator")
  omic_list <- sc_omicData(c("scRNA-seq","scATAC-seq"))
  cell_types <- list(cellA = c(1:20), cellB = c(161:191))
  sim <-scMOSim(omic_list, cell_types, numberCells = c(10,20), mean = c(2*10^6, 2*10^3), sd = c(10^3, 10^2))
  cell_types <- list(cellA= c(1:10), cellB = c(11:30))
  integration <- sc_omicSim(sim, cell_types, totalFeatures = 500)
  expect_equal(integration[["markers_cellA_cellB"]]$activity[1:10], expected)
})


test_that("sc_omicSim also works for a simulation with 4 cell types",{
  rna_counts <- readRDS("./data/RNA_4CellTypes.rds")
  atac_counts <- readRDS("./data/ATAC_4CellTypes.rds")
  omicData_4cell <- sc_omicData(omics_types =  list("scRNA-seq","scATAC-seq"), data = list(rna_counts,atac_counts))
  cell_types <- list(cellA = c(1:10), cellB = c(161:171), cellC =c(271:281), cellD=c(290:300))
  sim_4cell <- scMOSim(omicData_4cell, cell_types, numberCells = c(10,20,10,20), mean = c(2*10^6, 2*10^3, 2*10^5, 2*10^6), sd = c(10^3, 10^2, 10^2, 10^3))
  cell_types_integration <- list(cellA = c(1:10), cellB = c(11:30), cellC =c(31:40), cellD=c(41:60))
  expect_type(sc_omicSim(sim_4cell, cell_types_integration, totalFeatures = 500), "list")
})


test_that("checking if sc_omicSim does not return two identical list for markers_cellA_cellB and markers_cellB_cellC",{
  rna_counts <- readRDS("./data/RNA_4CellTypes.rds")
  atac_counts <- readRDS("./data/ATAC_4CellTypes.rds")
  omicData_4cell <- sc_omicData(omics_types =  list("scRNA-seq","scATAC-seq"), data = list(rna_counts,atac_counts))
  cell_types <- list(cellA = c(1:10), cellB = c(161:171), cellC =c(271:281), cellD=c(290:300))
  sim_4cell <- scMOSim(omicData_4cell, cell_types, numberCells = c(10,20,10,20), mean = c(2*10^6, 2*10^3, 2*10^5, 2*10^6), sd = c(10^3, 10^2, 10^2, 10^3))
  cell_types_integration <- list(cellA = c(1:10), cellB = c(11:30), cellC =c(31:40), cellD=c(41:60))
  omicSim <- sc_omicSim(sim_4cell, cell_types_integration, totalFeatures = 500)
  expect_true(!identical(cell_types_integration[[1]][[1]], cell_types_integration[[2]][[1]]))
})


test_that("checking that passing the argument regulatorEffect to sc_omicSim subsets the activity column in the correct way",{
  rna_counts <- readRDS("./data/RNA_4CellTypes.rds")
  atac_counts <- readRDS("./data/ATAC_4CellTypes.rds")
  omicData_4cell <- sc_omicData(omics_types =  list("scRNA-seq","scATAC-seq"), data = list(rna_counts,atac_counts))
  cell_types <- list(cellA = c(1:10), cellB = c(161:171), cellC =c(271:281), cellD=c(290:300))
  sim_4cell <- scMOSim(omicData_4cell, cell_types, numberCells = c(10,20,10,20), mean = c(2*10^6, 2*10^3, 2*10^5, 2*10^6), sd = c(10^3, 10^2, 10^2, 10^3))
  
  cell_types_integration <- list(cellA = c(1:10), cellB = c(11:30), cellC =c(31:40), cellD=c(41:60))
  regulatorEffect <- list('activator' = 0.8,'repressor' = 0.1,'NE' = 0.1)
  omicSim <- sc_omicSim(sim_4cell, cell_types_integration, totalFeatures = 500)
  omicSim_subset <- sc_omicSim(sim_4cell, cell_types_integration, totalFeatures = 500, regulatorEffect = regulatorEffect)
  NE <- sum(omicSim[["markers_cellA_cellB"]][["activity"]]=="NE")
  NE_rounded <- round(NE*0.1)
  NE_subset <- sum(omicSim_subset[["markers_cellA_cellB_subset"]][["activity"]] =="NE")
  expect_true(NE_rounded == NE_subset)
})
