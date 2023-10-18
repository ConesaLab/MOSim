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
  options(Seurat.object.assay.version = "v3")
  Seurat_obj <- Seurat::CreateAssayObject(counts = count, assay = 'RNA')
  testthat::expect_type(sc_omicData(c("scRNA-seq"), c(Seurat_obj)), "list")
})

testthat::test_that("Passing an array ('scRNA-seq','scATAC-seq') as omic returns a list of length 2", {
  res<-MOSim::sc_omicData(list("scRNA-seq","scATAC-seq"))
  testthat::expect_equal(length(res), 2)
})

#scMOSim function
testthat::test_that("param_estimation returns a list", {
  omic_list <- MOSim::sc_omicData(c("scRNA-seq"))
  conditions <- list('CD4_TEM' = c(1:60), 'cDC' = c(299:310))
  testthat::expect_type(MOSim::scMOSim(omic_list, conditions, numberCells = c(10,20), 
                          mean = c(2*10^6, 2*10^3), sd = c(10^3, 10^2)),"list")
})


testthat::test_that("Not passing all optional arguments at once returns an error", {
  omic_list <- MOSim::sc_omicData(c("scRNA-seq","scATAC-seq"))
  conditions <- list('CD4_TEM' = c(1:60), 'cDC' = c(299:310))
  testthat::expect_error(MOSim::scMOSim(omic_list, conditions, numberCells = c(10,20), sd = c(10^3, 10^2)))
})

testthat::test_that("scMOSim returns a list with S4 obj as values", {
  omic_list <- MOSim::sc_omicData(c("scRNA-seq"))
  cell_types <- list('CD4_TEM' = c(1:60), 'cDC' = c(299:310))
  sim <-MOSim::scMOSim(omic_list, cell_types, numberCells = c(10,20), mean = c(2*10^6, 2*10^3), sd = c(10^3, 10^2))
  testthat::expect_type(sim[[1]][[1]][[1]], "S4")
})


testthat::test_that("If numberGroups > 1, length(diffGenes) must be == (numberGroups -1)", {
  omicsList <- MOSim::sc_omicData(list("scRNA-seq"))
  cell_types <- list('CD4_TEM' = c(1:60), 'cDC' = c(299:310), 'Memory_B' = c(497:510), 'Treg' = c(868:900))
  # For this not to be an error, it should be diffGenes = list(c(0.2, 0.2))
  testthat::expect_error(MOSim::scMOSim(omicsList, cell_types, numberReps = 2, numberGroups = 2, 
                       diffGenes = c(0.2, 0.2), minFC = 0.25, maxFC = 4,
                       numberCells = NULL, mean = NULL, sd = NULL))
})

#make_cluster_patterns function
testthat::test_that("make_cluster_patterns returns a tibble of dimensions 
          numcells^2xnumcells",{
            patterns <- MOSim::make_cluster_patterns(4, 4)
            expected <- c(4, 4)
            testthat::expect_equal(dim(patterns$patterns), expected)
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

## Check we are generating matrices
testthat::test_that("The number of patterns we want are generated", {
  omic_list <- MOSim::sc_omicData(c("scRNA-seq", "scATAC-seq"))
  cell_types <- list('CD4_TEM' = c(1:60), 'cDC' = c(299:310), 
                     'Memory_B' = c(497:520), 'Treg' = c(868:900))
  sim <-MOSim::scMOSim(omic_list, cell_types, regulatorEffect = list(c(0.1, 0.2)))
  
  sim_matrix <- sim$Group_1$Rep_1$`sim_scRNA-seq`$counts
  testthat::expect_type(sim_matrix, "S4")
})

## test passing association list and working with groups and replicates
testthat::test_that("checking that scMOSim is able to simulate groups and replicates", {
  omicsList <- MOSim::sc_omicData(list("scRNA-seq", "scATAC-seq"))
  cell_types <- list('CD4_TEM' = c(1:60), 'cDC' = c(299:310), 'Memory_B' = c(497:510), 'Treg' = c(868:900))
  data("associationList")
  testing_groupsreps <- MOSim::scMOSim(omicsList, cell_types, numberReps = 2, numberGroups = 2, 
                                       diffGenes = list(c(0.1, 0.2)), minFC = 0.25, maxFC = 4,
                                       numberCells = NULL, mean = NULL, sd = NULL, 
                                       regulatorEffect = list(c(0.1, 0.2), c(0.2, 0.3)),
                                       associationList = associationList)
  testthat::expect_type(testing_groupsreps, "list")
})

testthat::test_that("checking that scMOSim is able to simulate more than 2 replicates", {
  omicsList <- MOSim::sc_omicData(list("scRNA-seq", "scATAC-seq"))
  cell_types <- list('CD4_TEM' = c(1:60), 'cDC' = c(299:310), 'Memory_B' = c(497:510), 'Treg' = c(868:900))
  testing_groupsreps <- MOSim::scMOSim(omicsList, cell_types, numberReps = 2, numberGroups = 2, 
                                       diffGenes = list(c(0.1, 0.2)), minFC = 0.25, maxFC = 4,
                                       numberCells = NULL, mean = NULL, sd = NULL, 
                                       regulatorEffect = list(c(0.1, 0.2), c(0.2, 0.3)),
                                       associationList = associationList)
  testthat::expect_type(testing_groupsreps, "list")
})

testthat::test_that("checking that scMOSim is able to simulate more than 2 groups", {
  omicsList <- MOSim::sc_omicData(list("scRNA-seq", "scATAC-seq"))
  cell_types <- list('CD4_TEM' = c(1:60), 'cDC' = c(299:310), 'Memory_B' = c(497:510), 'Treg' = c(868:900))
  testing_groups <- MOSim::scMOSim(omicsList, cell_types, numberReps = 1, numberGroups = 4, 
                                       diffGenes = list(c(0.1, 0.2), c(0.1, 0.2), c(0.1, 0.2)), minFC = 0.25, maxFC = 4,
                                       numberCells = NULL, mean = NULL, sd = NULL, 
                                       regulatorEffect = list(c(0.1, 0.2), c(0.2, 0.3), c(0.1, 0.2), c(0.2, 0.3)),
                                       associationList = associationList)
  testthat::expect_type(testing_groupsreps, "list")
})
