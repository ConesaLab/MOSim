test_that("simulation with rna-seq completes", {
    expect_s4_class(mosim(omics = "RNA-seq", omicsOptions = list("RNA-seq" = list(
        totalFeatures = 200
    ))), "Simulation")
})
