test_that("profileprobs settings", {

    no_induction <- list(
        continuous.induction = 0,
        continuous.repression = .2,
        transitory.induction = 0,
        transitory.repression = .2,
        flat = .2
    )

    sim_no_induction <- mosim(omics = "RNA-seq", omicsOptions = list(
        "RNA-seq"=list(totalFeatures=300)), profileProbs = no_induction)

    # Check the presence of induction profiles in settings
    expression_profiles <- sim_no_induction@simSettings$geneProfiles$SimRNAseq

    if (any(grepl("induction", expression_profiles$Group1),
            grepl("induction", expression_profiles$Group2)))
        fail("detected invalid expression profiles")

    succeed()
})
