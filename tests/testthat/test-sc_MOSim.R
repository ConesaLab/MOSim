# sc_omicData function
test_that("Passing a wrong string in 'omics' returns NA", {
  expect_message(sc_omicData("scR-seq", rna_orig_counts), NA)
})

test_that("Passing an object in data which is neither a matrix or Seurat obj returns NA", {
  vector <- c(1,2,3)
  expect_message(sc_omicData("scATAC-seq", vector), NA)
})

test_that("Passing 'scATAC-seq' as omic returns list", {
  expect_type(sc_omicData("scATAC-seq"), "list")
})

test_that("Passing 'scRNA-seq' as omic and a Seurat obj as data returns a list", {
  scRNA <- sc_omicData("scRNA-seq")
  count <- scRNA[["scRNA-seq"]]
  Seurat_obj <- CreateSeuratObject(counts = count, assay = 'RNA')
  expect_type(sc_omicData("scRNA-seq", Seurat_obj), "list")
})


#param_estimation function
test_that("param_estimation returns a list", {
  scRNA <- sc_omicData("scRNA-seq")
  scATAC <- sc_omicData("scATAC-seq")
  omic_list <- c(scRNA, scATAC)
  conditions <- list(cellA = c(1:20), cellB = c(161:191))
  expect_type(param_estimation(omic_list, conditions, numberCells = c(10,20), 
                               mean = c(2*10^6, 2*10^3), sd = c(10^3, 10^2)),"list")
})

test_that("Not passing all optional arguments at once returns NA", {
  scRNA <- sc_omicData("scRNA-seq")
  scATAC <- sc_omicData("scATAC-seq")
  omic_list <- c(scRNA, scATAC)
  conditions <- list(cellA = c(1:20), cellB = c(161:191))
  expect_message(param_estimation(omic_list, conditions, 
                                  numberCells = c(10,20),sd = c(10^3, 10^2)), NA)
})


#sc_MOSim function
test_that("scMOSim returns a list with S4 obj as values", {
  scRNA <- sc_omicData("scRNA-seq")
  scATAC <- sc_omicData("scATAC-seq")
  omic_list <- c(scRNA, scATAC)
  conditions <- list(cellA = c(1:20), cellB = c(161:191))
  sim <-scMOSim(omic_list, cell_types, numberCells = c(10,20), 
                mean = c(2*10^6, 2*10^3), sd = c(10^3, 10^2))
  expect_type(sim[[1]], "S4")
  
})

