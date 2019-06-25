context("Convenience functions")

# Test for convenience user functions
data(sampleData)

# Test data
# Select the first 10 rows
test_data_rna <- head(sampleData$SimRNAseq$data, 10)
test_data_dnase <- head(sampleData$SimDNaseseq$data, 10)
test_assoc_dnase <- dplyr::filter(sampleData$SimDNaseseq$idToGene,
                                  ID %in% rownames(test_data_dnase))

test_that("omicData with custom data", {
    expect_is(omicData("RNA-seq", data = test_data_rna), "SimRNAseq")
    expect_is(omicData("DNase-seq", data = test_data_dnase, associationList = test_assoc_dnase), "SimDNaseseq")
    expect_error(omicData("DNase-seq", data = test_data_dnase, associationList = NULL))
})

test_that("omicSim with multiple omics", {
    expect_named(omicSim("RNA-seq", totalFeatures = 20, depth = 15), "RNA-seq")
    expect_named(omicSim("DNase-seq", totalFeatures = 20, depth = 15), "DNase-seq")
})


test_that("mosim with custom data", {

    paramsRNAseq <- omicSim("RNA-seq", totalFeatures = 20, depth = 15)
    paramsDNaseseq <- omicSim("DNase-seq", totalFeatures = 20, depth = 15)

    mosim(
        omics = c("RNA-seq", "DNase-seq"),

        omicsOptions = c(paramsRNAseq, paramsDNaseseq),

        numberGroups = 2,
        times = c(1, 3, 7),
        numberReps = 3
    )

})
