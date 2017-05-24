test_that("simulation with rna-seq completes", {
    expect_s4_class(mosim(omics = "RNA-seq"), "Simulation")
})
