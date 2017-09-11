#' @include auxFunctions.R
#' @import parallel pbapply
NULL

#########################################################################################
######           Functions to integrate omics data using GLMs                      ######
#########################################################################################


## By Sonia & Monica
## 07-July-2016


# Generalized Lineal Model  -----------------------------------------------
# DETAILS
# In the case that the data matrix contains more variables than samples, a stepwise forward is applied, but the experimental
# design variables are kept.
#
# VALUE
# A list with 2 objects:
# - SummaryTable: A data.frame with one row for each gene supplied for the user. The data.frame has 4 columns:
#   - gene: gene IDs
#   - RP: A tag with 3 possibles values. If the gene is classified corretly in a regulatory program, the RP name,
#     if the gene is wrong classified, "WC", and if the gene is not classified in any of the regulatory program, "NC".
# - FinalResults: A list with one element for each gene studied. Each element of the list is also a list containing:
#  - GLMfinal: an object of class glm  with the best model obtained
#  - GLMorigen: an object of class glm with the starting model if we have enough df for performing at generalized linear models.
#    In other case, GLMfinal=GLMorigen
#  - SummaryStepwise: A named (by methologies stepwise applied) list. The list names could be one or some of the following
#    "stepfor","two.ways.stepfor" or "two.ways.stepback“. It depends on "stepwise" selection. Each one of these objects contain
#    a list with some summary parameters, "sol","coefficients","t.score","variables" and "edesign"
#    - sol, matrix for summary results of the stepwise regression. The following values are given:
#    p-value of the regression ANOVA, R-squared of the model, AIC of the model and p-value of the regression coefficients
#    of the selected variables.
#    - coefficients, regression coefficients for significant regulators
#    - t.score, value of the t statistics of significant regulators
#    - variables: variables in the final model
#    - edesign: matrix of experimental design. ¿LO QUITAMOS DE AQUI?
#  - RegulOrig: a list with 3 objects, with some summary information about original regulators in the model.
#    - OrigRegulators: The start regulators of the model (previous to make any stepwise). Ponemos los mismos que en RegulatorsValue??
#    - RegulatorsValue: The start regulators value of the model after applying CleanPredictors (previous to make any stepwise)
#    - edesign: matrix of experimental design
#    - cor.max: Correlation value to decide when multicollinearity is present.
#    - action: Action to perform when regulators are correlated. If "mean" (default), the average of the correlated regulators is
#              computed. If "random", one of the correlated regulators is randomly selected and the others are discarded.

#' Genome-wide Generalized Linear Models
#'
#'\code{GetGLM} fits a regression model for all the genes in the data set to identify
#' the experimental variables and potential regulators that show a significant effect on
#' the expression of each gene.
#'
#' @param GeneExpression Data frame containing gene expression data with genes in rows and
#' experimental samples in columns. Row names must be the gene IDs.
#' @param associations List where each element correspond to a different omic data type (miRNAs,
#' transcription factors, methylation, etc.). The names of the list will be the omics. Each element
#' is a data frame with 2 columns (optionally 3) describing the potential interactions between genes
#' and regulators for that omic. First column must contain the genes (or features in
#' GeneExpression object), second column must contain the regulators, and an additional column can
#' be added to describe the type of interaction (e.g., for methylation, if a CpG site is located in
#' the promoter region of the gene, in the first exon, etc.).
#   (optional).
#' @param data.omics List where each element correspond to a different omic data type (miRNAs,
#' transcription factors, methylation, etc.). The names of the list will be the omics. Each element
#' is a data matrix with omic regulators in rows and samples in columns.
#' @param edesign Data frame describing the experimental design. Rows must be the samples (columns
#' in GeneExpression) and columns must be the experimental variables to be included in the model
#' (e.g. time, treatment, etc.).
#' @param cont.var Name of the column in edesign that is to be considered as a continuous
#' explanatory variable in the regression model. NULL if edesign does not contain any continuous
#' variables. By default, "Time".
#' @param degree If cont.var is not NULL, non-linear (polynomial) relationships between the cont.var
#' and the gene expression must be studied. By default, degree = 1 and must be an integer.
#' A higher number will allow for in quadratic, cubic, etc. terms for cont.var in the model.
#' @param Res.df Number of degrees of freedom in the residuals. By default, 5. Increasing
#' Res.df will increase the power of the statistical model.
#' @param alfa Significance level. By default, 0.05.
#' @param MT.adjust Multiple testing correction method to be used within the stepwise variable
#' selection procedure. By default, "none". See the different options in ?\code{p.adjust}.
#' @param family Error distribution and link function to be used in the model (see ?\code{glm}
#' for more information). By default, negative.binomial(theta=10).
#' @param stepwise Stepwise variable selection method to be applied. It can be one of: "backward"
#' (default), "forward", "two.ways.backward" or "two.ways.forward".
#' @param interactions.exp If TRUE (default), interactions among the experimental variables
#' are included in the model.
#' @param interactions.reg If TRUE (default), interactions between regulators and experimental
#' variables are included in the model.
#' @param min.variation     ### if NULL, it computes "sd"
#' @param correlation
#' @param action
#' @param min.obs
#' @param elasticnet  0 = No ElasticNet variable selection; 1 = ElasticNet with ad-hoc penalization; 2 = ElasticNet with optimum penalization
#'
#' @return
#' @export
#'
#' @examples
GetGLM = function(GeneExpression,
                  associations,
                  data.omics,
                  edesign,
                  cont.var = "Time",
                  degree = 1,
                  Res.df = 5,
                  alfa = 0.05, MT.adjust = "none",
                  family = negative.binomial(theta=10),
                  elasticnet = 1,
                  stepwise = "backward",
                  interactions.exp = TRUE, interactions.reg = 1,
                  min.variation = 0,
                  correlation = 0.9, action = "mean",
                  min.obs = 10,
                  epsilon = 0.00001,
                  mc.cores = 1){


    # Converting matrix to data.frame
    GeneExpression = as.data.frame(GeneExpression)
    data.omics = lapply(data.omics, as.data.frame)


    # Preparing family for ElasticNet variable selection
    family2 = family$family
    family2 = strsplit(family2, "(", fixed = TRUE)[[1]][1]

    if (family2 %in% c("poisson", "quasipoisson", "Negative Binomial")) {
        family2 = "poisson"
    } else if (family2 %in% c("gaussian", "binomial")) {
        family2 = family2
    } else {
        family2 = NULL
        message(sprintf("Warning message:"))
        message(sprintf(
            "Elasticnet variable selection cannot be applied for family %s",
            family2
        ))
    }

  # Force numeric on edesign columns
  edesign[, cont.var] = apply(edesign[, cont.var, drop = F], 2, function(cvar) {
      if (is.factor(cvar)) {
          warning(sprintf("The cont.var column '%s' is a factor. The conversion to numeric type will use the underlying numeric representation.", cvar))
      }

      return(as.numeric(cvar))
  })


  # Checking that Res.df is coherent with the number of samples
  if (Res.df >= (nrow(edesign)-1)) stop("ERROR: You must decrease the Res.df so that a model can be computed.")


  # Checking that the number of samples per omic is equal to number of samples for gene expression
  for (i in 1:length(names(data.omics))){
    if(!length(colnames(data.omics[[i]]))== length(colnames(GeneExpression))){
      stop("ERROR: Samples in data.omics must be the same as in GeneExpression")
    }
  }

  ## Removing genes with too many NAs and keeping track
  genesNotNA = apply(GeneExpression, 1, function (x) sum(!is.na(x)))
  genesNotNA = names(which(genesNotNA >= min.obs))
  genesNA = setdiff(rownames(GeneExpression), genesNotNA)
  GeneExpression = GeneExpression[genesNotNA,]
  Allgenes=rownames(GeneExpression)
  nGenes = length(Allgenes)

  # Experimental groups
  ExpGroups = apply(edesign, 1, paste, collapse = "_")

  ## Experimental design matrix (with polynomial terms for cont.var, dummies for factors and interactions)
  des.mat = GenerateDesignMatrix(interactions.exp, degree, edesign, cont.var)

  # Experimental design matrix (with polynomial terms for cont.var. without dummies and interactions)
  # if (interactions.exp) {  # we compute all possible interactions between exp variables
  #
  #   if ((cont.var %in% colnames(edesign)) && (degree > 1)) {
  #     des.mat = data.frame(edesign[,-which(colnames(edesign) == cont.var), drop = FALSE])
  #     if (NCOL(des.mat) > 1) {
  #       fff = paste0("~ ", paste(sprintf("`%s`", colnames(des.mat)), collapse = "*"))
  #       fff = as.formula(fff)
  #       des.mat = model.matrix(fff, des.mat)[,-1]  # interactions without time
  #     }
  #     expcond = colnames(des.mat)
  #     # Adding time terms
  #     politerms = c(cont.var, paste(cont.var, 2:degree, sep = ""))
  #     for (i in 1:degree) {
  #       des.mat = data.frame(des.mat, edesign[,cont.var]^i)
  #     }
  #     # Adding interactions with time
  #     colnames(des.mat)[(ncol(des.mat)-degree+1):ncol(des.mat)] = politerms
  #     fff = paste0("~ ", paste(sapply(politerms, function (x) paste(expcond, x, sep = ":")), collapse = "+"))
  #     fff = as.formula(fff)
  #     des.mat = cbind(des.mat, model.matrix(fff, des.mat)[,-1])
  #     expcond = colnames(des.mat)
  #
  #   } else {  # No cont.var or degree=1
  #
  #     cont.var = NULL
  #
  #     if (NCOL(edesign) > 1) {
  #       des.mat = edesign
  #       fff = paste0("~ ", paste(sprintf("`%s`", colnames(des.mat)), collapse = "*"))
  #       fff = as.formula(fff)
  #       des.mat = model.matrix(fff, des.mat)[,-1]  # interactions between exp variables
  #     } else { # only 1 exp variable
  #       des.mat = edesign
  #     }
  #     expcond = colnames(des.mat)
  #   }
  #
  # } else {   ## NO interactions between experimental variables
  #
  #   if ((cont.var %in% colnames(edesign)) && (degree > 1)) {
  #     des.mat = data.frame(edesign[,-which(colnames(edesign) == cont.var), drop = FALSE])
  #     if (NCOL(des.mat) > 0) {
  #       # fff = paste0("~ ", paste(sprintf("`%s`", colnames(des.mat)), collapse = "*")) ## NO interactions at all!!!
  #       # fff = as.formula(fff)
  #       des.mat = model.matrix(~., des.mat)[,-1, drop = FALSE]  ## to convert factors into dummy variables
  #     }
  #     expcond = colnames(des.mat)
  #     # Adding time terms
  #     politerms = c(cont.var, paste(cont.var, 2:degree, sep = ""))
  #     for (i in 1:degree) {
  #       des.mat = data.frame(des.mat, edesign[,cont.var]^i)
  #     }
  #     colnames(des.mat)[(ncol(des.mat)-degree+1):ncol(des.mat)] = politerms
  #     expcond = colnames(des.mat)
  #
  #   } else {  # No cont.var or degree=1
  #
  #     cont.var = NULL
  #     # if (NCOL(edesign) > 1) {
  #       des.mat = edesign
  #       # fff = paste0("~ ", paste(sprintf("`%s`", colnames(des.mat)), collapse = "*"))
  #       # fff = as.formula(fff)
  #       # des.mat = model.matrix(fff, des.mat)[,-1]  # interactions between exp variables
  #       des.mat = model.matrix(~., des.mat)[,-1, drop = FALSE]  # to convert factors into dummy variables
  #     # } else { # only 1 exp variable
  #       #   des.mat = edesign
  #     # }
  #     expcond = colnames(des.mat)
  #   }
  # }


  ## Remove regulators with NA

  cat("Removing regulators with missing values...\n")

  myregNA = lapply(data.omics, rownames)
  data.omics = lapply(data.omics, na.omit)
  myregNA = lapply(1:length(data.omics), function (i) setdiff(myregNA[[i]], rownames(data.omics[[i]])))
  names(myregNA)=names(data.omics)

  cat("Number of regulators with missing values:\n")
  print(sapply(myregNA, length))
  cat("\n")


  ## Remove regulators with Low Variability
  cat("Removing regulators with low variation...\n")

  tmp = LowVariationRegu(min.variation, data.omics, ExpGroups, associations, Allgenes)
  data.omics = tmp[["data.omics"]]
  associations = tmp[["associations"]]
  myregLV = tmp[["myregLV"]]
  rm("tmp"); gc()

  # if (is.null(min.variation)){
  #
  #   #### Low variation cutoff is computed automatically
  #
  #   data.omicsMean=vector("list", length=length(data.omics))
  #
  #   for(i in 1:length(data.omics)){
  #     data.omicsMean[[i]]=t(apply(data.omics[[i]], 1, tapply, ExpGroups, mean))
  #   }
  #   names(data.omicsMean)=names(data.omics)
  #
  #   percVar=rep(10, length(data.omicsMean)) ## We fix it
  #   names(percVar)=names(data.omicsMean)
  #
  #   # Applying Low Variation filter
  #   LowVar=LowVariatFilter(data=data.omicsMean, method="sd", percVar=percVar)
  #
  #   data.omicsMean=LowVar$data  ## data.omicsMean reduced: without NA and LV
  #
  #   for (ov in names(associations)){
  #     myreg=associations[[ov]][associations[[ov]][,1] %in% Allgenes,2] # removing regulators not associated to our genes
  #     data.omicsMean[[ov]]=data.omicsMean[[ov]][intersect(myreg, rownames(data.omicsMean[[ov]])),]
  #     data.omics[[ov]]= data.omics[[ov]][rownames(data.omicsMean[[ov]]),] ## Reduced data.omics
  #   }
  #
  #
  # }
  # else {
  #
  #   #### Low variation cutoff is set by the user
  #
  #   for (ov in names(associations)){
  #     myreg=associations[[ov]][associations[[ov]][,1] %in% Allgenes,2] # removing regulators not associated to our genes
  #     data.omics[[ov]]=data.omics[[ov]][intersect(myreg, rownames(data.omics[[ov]])),]
  #   }
  #
  #
  #   if (length(min.variation) == 1) {  ## Including min.variation=0. I need a vector with omics names
  #     min.variation=rep(min.variation,length(data.omics))
  #     names(min.variation)=names(data.omics)
  #   }
  #
  #   data.omicsMean=vector("list", length = length(data.omics))  # computing mean per condition in data.omics
  #   for(i in 1:length(data.omics)){
  #     data.omicsMean[[i]] = t(apply(data.omics[[i]], 1, tapply, ExpGroups, mean))
  #   }
  #   names(data.omicsMean) = names(data.omics)
  #
  #   # Applying Low Variation filter
  #   LowVar=LowVariatFilter(data = data.omicsMean, method = "user", percVar = min.variation)
  #
  #   data.omicsMean=LowVar$data  ## data.omicsMean reduced
  #
  #   ## data.omics reduced: only mygenes and without NA and LV
  #   for (ov in names(data.omics)){
  #     data.omics[[ov]] = data.omics[[ov]][rownames(data.omicsMean[[ov]]),]
  #   }
  #
  # }
  #
  # rm("myreg"); rm("data.omicsMean"); gc()
  #
  #
  # # Regulators removed due to low variation filter
  # myregLV=LowVar$LV.reg
  # rm("LowVar"); gc()
  #
  # cat("Number of regulators with low variation:\n")
  # print(sapply(myregLV, length))
  # cat("\n")



  ### Results objects

  ## Global summary for all genes
  GlobalSummary = vector("list", length = 3)
  names(GlobalSummary) = c("GoodnessOfFit", "ReguPerGene", "GenesNOmodel")

  if (length(genesNA) > 0) {
      GlobalSummary$GenesNOmodel = data.frame("gene" = genesNA, "problem" = rep("Too many missing values", length(genesNA)))
  } else {
      GlobalSummary$GenesNOmodel = NULL
  }

  GlobalSummary$GoodnessOfFit = matrix(NA, ncol = 5, nrow = nGenes)
  rownames(GlobalSummary$GoodnessOfFit) = Allgenes
  colnames(GlobalSummary$GoodnessOfFit) = c("modelPvalue", "dfResiduals", "Rsquared", "AIC", "sigReg")

  GlobalSummary$ReguPerGene = matrix(0, ncol = 3*length(data.omics), nrow = nGenes)
  rownames(GlobalSummary$ReguPerGene) = Allgenes
  colnames(GlobalSummary$ReguPerGene) = c(paste(names(data.omics), "Ini", sep = "-"),
                                          paste(names(data.omics), "Mod", sep = "-"),
                                          paste(names(data.omics), "Sig", sep = "-"))

  ## Specific results for each gene
  ResultsPerGene=vector("list", length=length(Allgenes))
  names(ResultsPerGene) = Allgenes

  ### Computing model for each gene
  cat("Checking multicollinearity, selecting predictors and fitting model...\n")

  ## Specific results for each gene
  ResultsPerGene <- pbapply::pblapply(1:nGenes, FUN = function(i) {

    gene=Allgenes[i]

    ResultsPerGene.i = vector("list", length = 8)
    names(ResultsPerGene.i) = c("Y", "X", "coefficients", "allRegulators", "significantRegulators",
                                "GoodnessOfFit", "ReguPerGene", "GenesNOmodel")


    # Initialize global summary values
    ResultsPerGene.i$ReguPerGene <- GlobalSummary$ReguPerGene[gene, , drop = FALSE]


    # cat(paste("Gene", i, "out of", nGenes))
    # cat("\n")

    RetRegul = GetAllReg(gene=gene, associations=associations)
    RetRegul.gene = RetRegul$Results  ## RetRegul$TableGene: nr reg per omic
    ## Some of these reg will be removed, because they are not in data.omics


    # RetRegul.gene--> gene/regulator/omic/area
    RetRegul.gene=RetRegul.gene[RetRegul.gene[,"regulator"]!= "No-regulator", ,drop=FALSE] ## Remove rows with no-regulators


    ### NO INITIAL REGULATORS
    if(length(RetRegul.gene)==0){ ## En el caso de que no hayan INICIALMENTE reguladores -> Calcular modelo con variables experiment
      des.mat2 = cbind(t(GeneExpression[gene,]), des.mat)
      colnames(des.mat2)[1] = "response"
      des.mat2 = na.omit(des.mat2)

      # Removing predictors with constant values
      sdNo0 = apply(des.mat2, 2, sd)
      sdNo0 = names(sdNo0)[sdNo0 > 0]
      des.mat2 = des.mat2[,sdNo0]

      myGLM = ComputeGLM(matrix.temp = data.frame(des.mat2, check.names = FALSE),
                         alfa = alfa, stepwise = stepwise, Res.df = Res.df,
                         family = family, epsilon = epsilon, MT.adjust = MT.adjust)

      ResultsPerGene.i$X = des.mat2[, -1, drop = FALSE]
      ResultsPerGene.i$significantRegulators = NULL
      ResultsPerGene.i$allRegulators = NULL

      # GlobalSummary$ReguPerGene  # this is initially set to 0 so no need to modify it


      ### WITH INITIAL REGULATORS
    }
    else { ## There are regulators for this gene at the beginning

      ResultsPerGene.i$allRegulators = data.frame(RetRegul.gene, rep("Model",nrow(RetRegul.gene)), stringsAsFactors = FALSE)
      colnames(ResultsPerGene.i$allRegulators) = c("gene","regulator","omic","area","filter")

      ResultsPerGene.i$ReguPerGene[1, grep("-Ini", colnames(ResultsPerGene.i$ReguPerGene))] = as.numeric(RetRegul$TableGene[-1])
      # the rest of columns remain 0

      ## Identify which regulators where removed because of missing values or low variation
      res = RemovedRegulators(RetRegul.gene = ResultsPerGene.i$allRegulators,
                              myregLV=myregLV, myregNA=myregNA, data.omics=data.omics)

      if(length(res$RegulatorMatrix)==0){ ## No regulators left after the filtering to compute the model

        des.mat2 = cbind(t(GeneExpression[gene,]), des.mat)
        colnames(des.mat2)[1] = "response"
        des.mat2 = na.omit(des.mat2)

        # Removing predictors with constant values
        sdNo0 = apply(des.mat2, 2, sd)
        sdNo0 = names(sdNo0)[sdNo0 > 0]
        des.mat2 = des.mat2[,sdNo0]

        myGLM = ComputeGLM(matrix.temp = data.frame(des.mat2, check.names = FALSE),
                           alfa = alfa, stepwise = stepwise, Res.df = Res.df,
                           family = family, epsilon = epsilon, MT.adjust = MT.adjust)

        ResultsPerGene.i$X = des.mat2[, -1, drop = FALSE]
        ResultsPerGene.i$significantRegulators = NULL
        ResultsPerGene.i$allRegulators = res$SummaryPerGene

      }
      else {  ## Regulators for the model!!

        ## Multicollinearity
        res = CollinearityFilter(data = res$RegulatorMatrix, reg.table = res$SummaryPerGene,
                                 correlation = correlation, action = action)

        ResultsPerGene.i$allRegulators = res$SummaryPerGene


        ### Creating data matrix with regulators and with/without interactions
        des.mat2 = RegulatorsInteractions(interactions.reg,
                                          reguValues = res$RegulatorMatrix,
                                          des.mat, cont.var, GeneExpression,
                                          gene)

        # Removing observations with missing values
        des.mat2 = na.omit(des.mat2)


        ###  Variable selection --> Elasticnet
        tmp = ElasticNet(family2, des.mat2, epsilon, elasticnet, Res.df)

        removedCoefs = intersect(tmp[["removedCoefs"]], rownames(ResultsPerGene.i$allRegulators))
        if (length(removedCoefs) > 0)  ResultsPerGene.i$allRegulators[removedCoefs,"filter"] = "ElasticNet"

        des.mat2 = as.data.frame(tmp[["des.mat2"]])
        ResultsPerGene.i$X = des.mat2[, -1, drop = FALSE]

        rm(tmp); rm(removedCoefs); gc()

        ###################################
        if (ncol(des.mat2) > 1) {
            # Removing predictors with constant values
            sdNo0 = apply(des.mat2, 2, sd)
            sdNo0 = names(sdNo0)[sdNo0 > 0]
            des.mat2 = des.mat2[,sdNo0]

            ## Computing GLM model
            myGLM = try(ComputeGLM(matrix.temp = des.mat2,
                                   alfa = alfa, stepwise = stepwise, Res.df = Res.df,
                                   family = family, epsilon = epsilon, MT.adjust = MT.adjust), silent = TRUE)

            if (class(myGLM) == "try-error") {

                myGLM = NULL

                # GlobalSummary$GenesNOmodel = rbind(GlobalSummary$GenesNOmodel,
                #                                    data.frame("gene" = gene, "problem" = "GLM error"))
                ResultsPerGene.i$GenesNOmodel = data.frame("gene" = gene, "problem" = "GLM error")

                ## Extracting significant regulators and recovering correlated regulators
                ResultsPerGene.i$significantRegulators = NULL
                ResultsPerGene.i$allRegulators = data.frame(res$SummaryPerGene, "Sig" = NA, stringsAsFactors = FALSE)

                ## Counting original regulators in the model per omic
                contando =  ResultsPerGene.i$allRegulators[which(ResultsPerGene.i$allRegulators[,"filter"] == "Model"),]
                contando = table(contando[,"omic"])
                contando = as.numeric(contando[names(data.omics)])
                contando[is.na(contando)] = 0
                ResultsPerGene.i$ReguPerGene[gene, grep("-Mod", colnames(ResultsPerGene.i$ReguPerGene))] = contando

                ## Counting significant regulators per omic
                ResultsPerGene.i$ReguPerGene[1, grep("-Sig", colnames(ResultsPerGene.i$ReguPerGene))] = NA

            } else {

                ## Extracting significant regulators and recovering correlated regulators
                myvariables = unlist(strsplit(myGLM$SummaryStepwise$variables, ":", fixed = TRUE))
                myvariables = intersect(myvariables, rownames(ResultsPerGene.i$allRegulators))
                ResultsPerGene.i$significantRegulators = myvariables # significant regulators including "new" correlated regulators
                ResultsPerGene.i$allRegulators = data.frame(res$SummaryPerGene, "Sig" = 0, stringsAsFactors = FALSE)
                ResultsPerGene.i$allRegulators[myvariables, "Sig"] = 1
                collin.regulators = intersect(ResultsPerGene.i$significantRegulators, ResultsPerGene.i$allRegulators[,"filter"]) # "new" regulators
                if (length(collin.regulators) > 0) {  # there were correlated regulators
                    original.regulators = ResultsPerGene.i$allRegulators[ResultsPerGene.i$allRegulators[,"filter"] %in% collin.regulators,"regulator"]
                    ResultsPerGene.i$allRegulators[original.regulators, "Sig"] = 1
                    ResultsPerGene.i$significantRegulators = c(ResultsPerGene.i$significantRegulators, original.regulators)
                    if (action == "mean") {
                        ResultsPerGene.i$significantRegulators = setdiff(ResultsPerGene.i$significantRegulators, collin.regulators)
                    }
                }

                ## Counting original regulators in the model per omic
                if (length(collin.regulators) > 0) {
                    contando = ResultsPerGene.i$allRegulators
                    quitar = which(contando[,"filter"] == "MissingValue")
                    if (length(quitar) > 0) contando = contando[-quitar,]
                    quitar = which(contando[,"filter"] == "LowVariation")
                    if (length(quitar) > 0) contando = contando[-quitar,]
                    contando = contando[setdiff(rownames(contando), collin.regulators),]
                }
                else {
                    contando = ResultsPerGene.i$allRegulators[which(ResultsPerGene.i$allRegulators[,"filter"] == "Model"),]
                }
                contando = table(contando[,"omic"])
                contando = as.numeric(contando[names(data.omics)])
                contando[is.na(contando)] = 0
                ResultsPerGene.i$ReguPerGene[gene, grep("-Mod", colnames(ResultsPerGene.i$ReguPerGene))] = contando

                ## Counting significant regulators per omic
                if (length(ResultsPerGene.i$significantRegulators) > 0) {
                    contando = ResultsPerGene.i$allRegulators[ResultsPerGene.i$significantRegulators,]
                    contando = table(contando[,"omic"])
                    contando = as.numeric(contando[names(data.omics)])
                    contando[is.na(contando)] = 0
                    ResultsPerGene.i$ReguPerGene[1, grep("-Sig", colnames(ResultsPerGene.i$ReguPerGene))] = contando
                }

            }



        } else {  ## NO variables in the model because of ElasticNet selection

            myGLM = NULL

            # GlobalSummary$GenesNOmodel = rbind(GlobalSummary$GenesNOmodel,
            #                                    data.frame("gene" = gene, "problem" = "No predictors after EN"))
            ResultsPerGene.i$GenesNOmodel = data.frame("gene" = gene, "problem" = "No predictors after EN")

            ## Extracting significant regulators and recovering correlated regulators
            ResultsPerGene.i$significantRegulators = NULL
            ResultsPerGene.i$allRegulators = data.frame(ResultsPerGene.i$allRegulators, "Sig" = 0, stringsAsFactors = FALSE)

            ## Counting original regulators in the model per omic
            contando = ResultsPerGene.i$allRegulators[which(ResultsPerGene.i$allRegulators[,"filter"] == "Model"),]
            contando = table(contando[,"omic"])
            contando = as.numeric(contando[names(data.omics)])
            contando[is.na(contando)] = 0
            ResultsPerGene.i$ReguPerGene[gene, grep("-Mod", colnames(ResultsPerGene.i$ReguPerGene))] = contando

        }

        ####################################
#
#         if (interactions.reg) {  ### WITH INTERACTIONS with regulators
#
#           des.mat2 = data.frame(des.mat, res$RegulatorMatrix, check.names = FALSE)
#
#           fff = paste0("~ ",
#                        paste(sapply(colnames(res$RegulatorMatrix),
#                                     function (x) paste(expcond, sprintf("`%s`", x), sep = ":")),
#                              collapse = "+"))
#
#           fff = as.formula(fff)
#           des.mat2 = cbind(des.mat2, model.matrix(fff, des.mat2)[,-1])
#
#           des.mat2 = cbind(t(GeneExpression[gene,]), des.mat2)
#           colnames(des.mat2)[1] = "response"
#
#           colnames(des.mat2) = gsub("\`", "", colnames(des.mat2))
#
#
#           ## PUEDE OCURRIR QUE AL METER LAS INTERACCIONES TODA LA COLUMNA SEA 0. HAGO UN FILTRO PREVIO
#           sd.regulators = apply(des.mat2, 2, sd, na.rm=TRUE)
#           regulators0 = names(sd.regulators[sd.regulators==0])
#           if (length(regulators0)>0) des.mat2 = des.mat2[, setdiff(colnames(des.mat2), regulators0), drop=FALSE]
#
#         } else  {    ### WITHOUT INTERACTIONS
#
#           des.mat2 = data.frame(des.mat, res$RegulatorMatrix, check.names = FALSE)
#
#           des.mat2 = cbind(t(GeneExpression[gene,]), des.mat)
#           colnames(des.mat2)[1] = "response"
#
#         }
#
#         des.mat2 = na.omit(des.mat2)
#
#         ## It can happen that des.mat2 is not a full rank matrix --> remove dependent columns
#         des.mat3 = des.mat2[,-1, drop = FALSE]  # remove response variable
#
#         tmpPval = apply(des.mat3, 2,
#                         function (x) summary(glm(des.mat2[,1] ~ x, family = family, epsilon = epsilon))$coefficients[2,4])
#
#         tmpPval.ord = names(tmpPval)[order(tmpPval, decreasing = TRUE)]
#
#         while (Matrix::rankMatrix(des.mat3, tol = epsilon) < ncol(des.mat3)) {
#             # tmpPval = apply(des.mat3, 2,
#             #                 function (x) summary(glm(des.mat2[,1] ~ x, family = family, epsilon = epsilon))$coefficients[2,4])
#
#             # Remove the first N values with the highest p-value
#             # diff.rank = ncol(des.mat3) - des.rank
#
#             # des.mat3 = des.mat3[, - head(, diff.rank), drop = FALSE]
#             des.mat3 = des.mat3[, colnames(des.mat3) != tmpPval.ord[1], drop = FALSE]
#
#             tmpPval.ord = tail(tmpPval.ord, -1)
#         }
#
#         if (ncol(des.mat3) < (ncol(des.mat2)-1)) {
#             # print("Some explanatory variables had to be removed from the initial model to avoid multicollinearity problems.")
#             des.mat2 = data.frame(des.mat2[,1, drop = FALSE], des.mat3, check.names = FALSE)
#         }
#         rm(des.mat3); gc()
#
#         ResultsPerGene.i$X = des.mat2[,-1]
#
#         myGLM = ComputeGLM(matrix.temp = des.mat2,
#                            alfa = alfa, stepwise = stepwise, Res.df = Res.df,
#                            family = family, epsilon = epsilon, MT.adjust = MT.adjust)
#         # ResultsPerGene.i$X = des.mat2[,-1]
#         #
#         # myGLM = ComputeGLM(matrix.temp = data.frame(des.mat2, check.names = FALSE),
#         #                    alfa = alfa, stepwise = stepwise, Res.df = Res.df,
#         #                    family = family, epsilon = epsilon, MT.adjust = MT.adjust)
#
#
#         ## Extracting significant regulators and recovering correlated regulators
#         myvariables = unlist(strsplit(myGLM$SummaryStepwise$variables, ":", fixed = TRUE))
#         myvariables = intersect(myvariables, rownames(ResultsPerGene.i$allRegulators))
#         ResultsPerGene.i$significantRegulators = myvariables # significant regulators including "new" correlated regulators
#         ResultsPerGene.i$allRegulators = data.frame(res$SummaryPerGene, "Sig" = 0, stringsAsFactors = FALSE)
#         ResultsPerGene.i$allRegulators[myvariables, "Sig"] = 1
#         collin.regulators = intersect(ResultsPerGene.i$significantRegulators, ResultsPerGene.i$allRegulators[,"filter"]) # "new" regulators
#         if (length(collin.regulators) > 0) {  # there were correlated regulators
#           original.regulators = ResultsPerGene.i$allRegulators[ResultsPerGene.i$allRegulators[,"filter"] %in% collin.regulators,"regulator"]
#           ResultsPerGene.i$allRegulators[originResultsPerGene.ial.regulators, "Sig"] = 1
#           ResultsPerGene.i$significantRegulators = c(ResultsPerGene.i$significantRegulators, original.regulators)
#           if (action == "mean") {
#             ResultsPerGene.i$significantRegulators = setdiff(ResultsPerGene.i$significantRegulators, collin.regulators)
#           }
#         }
#
#         ## Counting original regulators in the model per omic
#         if (length(collin.regulators) > 0) {
#           contando = ResultsPerGene.i$allRegulators
#           quitar = which(contando[,"filter"] == "MissingValue")
#           if (length(quitar) > 0) contando = contando[-quitar,]
#           quitar = which(contando[,"filter"] == "LowVariation")
#           if (length(quitar) > 0) contando = contando[-quitar,]
#           contando = contando[setdiff(rownames(contando), collin.regulators),]
#         }
#         else {
#           contando = ResultsPerGene.i$allRegulators[which(ResultsPerGene.i$allRegulators[,"filter"] == "Model"),]
#         }
#         contando = table(contando[,"omic"])
#         contando = as.numeric(contando[names(data.omics)])
#         contando[is.na(contando)] = 0
#         GlobalSummary$ReguPerGene[gene, grep("-Mod", colnames(GlobalSummary$ReguPerGene))] = contando
#
#         ## Counting significant regulators per omic
#         if (length(ResultsPerGene.i$significantRegulators) > 0) {
#           contando = ResultsPerGene.i$allRegulators[ResultsPerGene.i$significantRegulators,]
#           contando = table(contando[,"omic"])
#           contando = as.numeric(contando[names(data.omics)])
#           contando[is.na(contando)] = 0
#           GlobalSummary$ReguPerGene[gene, grep("-Sig", colnames(GlobalSummary$ReguPerGene))] = contando
#         }


      } ## Close "else" --> None regulators from begining



      if (is.null(myGLM)) {

          ResultsPerGene.i$Y = des.mat2[,1]
          ResultsPerGene.i$coefficients = NULL

          ResultsPerGene.i$GoodnessOfFit = c(NA, NA, NA, NA, NA)

      } else {
          ResultsPerGene.i$Y = data.frame("y" = myGLM$GLMfinal$y, "fitted.y" = myGLM$GLMfinal$fitted.values,
                                             "residuals" = residuals(myGLM$GLMfinal))

          ResultsPerGene.i$coefficients = summary(myGLM$GLMfinal)$coefficients[,c(1,4), drop = FALSE]
          if (nrow(ResultsPerGene.i$coefficients) > 0) {
              colnames(ResultsPerGene.i$coefficients) = c("coefficient", "p-value")
          }

          ResultsPerGene.i$GoodnessOfFit = c(myGLM$SummaryStepwise$"p.value",
                                                 summary(myGLM$GLMfinal)$df.residual,
                                                 myGLM$SummaryStepwise$"R.squared",
                                                 summary(myGLM$GLMfinal)$aic,
                                                 length(ResultsPerGene.i$significantRegulators))
      }
    }

    #   ResultsPerGene.i$Y = data.frame("y" = myGLM$GLMfinal$y, "fitted.y" = myGLM$GLMfinal$fitted.values,
    #                                      "residuals" = residuals(myGLM$GLMfinal))
    #
    #   ResultsPerGene.i$coefficients = summary(myGLM$GLMfinal)$coefficients[,c("Estimate", "Pr(>|t|)"), drop = FALSE]
    #   if (nrow(ResultsPerGene.i$coefficients) > 0) {
    #     colnames(ResultsPerGene.i$coefficients) = c("coefficient", "p-value")
    #   }
    #
    #   ResultsPerGene.i$GoodnessOfFit <- c(myGLM$SummaryStepwise$"p.value",
    #                                            summary(myGLM$GLMfinal)$df.residual,
    #                                            myGLM$SummaryStepwise$"R.squared",
    #                                            summary(myGLM$GLMfinal)$aic,
    #                                            length(ResultsPerGene.i$significantRegulators))
    # }

    return(ResultsPerGene.i)

    }, cl = mc.cores)


  # Restore GoodnessOfFit to the proper data frame using a simply lapply (not multicore)
  globalValues <- c("GoodnessOfFit", "ReguPerGene", "GenesNOmodel")

  ResultsPerGene <- lapply(seq_along(ResultsPerGene), function(i) {

      for (gValue in globalValues) {
          if (exists(gValue, ResultsPerGene[[i]]) && ! is.null(ResultsPerGene[[i]][[gValue]])) {
              # GenesNOmodel should be rbinded to the global data frame.
              # TODO: skip data frame checking for something more reliable?
              if (! is.data.frame(ResultsPerGene[[i]][[gValue]])) {
                  GlobalSummary[[gValue]][Allgenes[i], ] <<- ResultsPerGene[[i]][[gValue]]
              } else {
                  GlobalSummary[[gValue]] <<- rbind(GlobalSummary[[gValue]], ResultsPerGene[[i]][[gValue]])
              }
          }
      }

      # Remove the global values keys
      return(ResultsPerGene[[i]][- grep(paste0(globalValues, collapse = "|"),
                                        names(ResultsPerGene[[i]]))])
  })

  names(ResultsPerGene) = Allgenes

  myarguments = list(edesign = edesign, degree = degree, Res.df = Res.df, alfa = alfa, family = family,
                     stepwise = stepwise, interactions.exp = interactions.exp, interactions.reg = interactions.reg,
                     min.variation = min.variation, correlation = correlation, action = action,
                     MT.adjust = MT.adjust, min.obs = min.obs, cont.var = cont.var, epsilon = epsilon,
                     dataOmics = data.omics)

  return(list("ResultsPerGene" = ResultsPerGene, "GlobalSummary" = GlobalSummary, "arguments" = myarguments))

}






# Get all regulators ------------------------------------------------------


## Auxiliar function to paste different areas

#' Title
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
uniquePaste = function (x) {
  x = unique(x)

  if (length(x) > 1) {
    x = paste(x, collapse = ";")
  }

  return (x)
}




## INPUT
## gene: we are seeking information
## associations: List containing as many elements as original association files between genes and each omic (association$miRNA, association$DNase, etc.)
## data.omics: For removing regulators with no-information

#' Title
#'
#' @param gene
#' @param associations
#'
#' @return
#' @export
#'
#' @examples
GetAllReg=function(gene, associations){
  Reg.matrix=NULL
  NrReg=NULL
  myomic=names(associations)

  for(ov in myomic){

    colnames(associations[[ov]])[1]="gene"

    ## Regulator with area

    if(ncol(associations[[ov]])>2){

      myregulators=associations[[ov]][associations[[ov]]$gene==gene, ,drop=FALSE] ## "regulator" with Area--> Matrix

      if(nrow(myregulators)==0){
        ov.nr=nrow(myregulators) ## regulators nr per omic --> It could be 0
        myregulators=c("No-regulator","")
        myregulators=t(as.matrix(myregulators))

      } else {
        myregulators=aggregate(myregulators,by=list(myregulators[,2]),uniquePaste) ## Agrupado por "region", me devuelve Area separado por comas
        myregulators=myregulators[,-c(1:2)]
        myregulators[,2]=sapply(myregulators[,2], paste, collapse = ";")
        ov.nr=nrow(myregulators)
      }

      Reg.matrix.temp = cbind(gene,myregulators,ov)
      Reg.matrix.temp = Reg.matrix.temp[,c(1,2,4,3),drop=FALSE] ## Dejo el mismo orden que tenía
      colnames(Reg.matrix.temp)=c("gene","regulator","omic","area")
      NrReg.temp=ov.nr
      Reg.matrix=rbind(Reg.matrix,Reg.matrix.temp)
      colnames(Reg.matrix)=c("gene","regulator","omic","area")
      NrReg=c(NrReg, NrReg.temp)

      ## Regulators without area

    } else {

      myregulators=associations[[ov]][associations[[ov]]$gene==gene,2, drop=FALSE] ## "regulator" --> No tiene Area
      myregulators=unique(myregulators) ## Could it be repeated??
      ov.nr=nrow(myregulators) ## regulators nr per omic --> It could be 0

      if(nrow(myregulators)==0){
        myregulators=c("No-regulator")
        myregulators=t(as.matrix(myregulators))
      }

      Reg.matrix.temp=cbind(gene,myregulators,ov,"")  ## Lo tengo que dejar igual que las otras omicas
      NrReg.temp=ov.nr
      colnames(Reg.matrix.temp)=c("gene","regulator","omic","area")
      Reg.matrix.temp = t(apply(Reg.matrix.temp,1, as.character))  ## adding this here to avoid factors
      Reg.matrix=rbind(Reg.matrix,Reg.matrix.temp)
      colnames(Reg.matrix)=c("gene","regulator","omic","area")
      NrReg=c(NrReg, NrReg.temp)
    }
  }

  Results2=c(gene,NrReg)
  myomic=paste(myomic,"Ini",sep="-")
  names(Results2)=c("gene",myomic)

  Results=vector("list", length=2)
  Results[[1]]=Reg.matrix
  Results[[2]]=Results2
  names(Results)=c("Results","TableGene")

  return(Results)

}






# Low variation filtering -------------------------------------------------

## data: For computing variability. It could be mean values between replicas
## method: One of "sd","range", "IQrange" or "user"
## percVar: percentage of variation defined by the user

#' Title
#'
#' @param data
#' @param method
#' @param percVar
#'
#' @return
#' @export
#'
#' @examples
LowVariatFilter=function(data, method, percVar){

  SummaryRes=vector("list", length=length(data))
  names(SummaryRes)=names(data)
  LV.reg=vector("list", length=length(data))
  names(LV.reg)=names(data)

  if(method=="sd"){

    for (ov in names(data)) {

      met=apply(data[[ov]], 1, sd, na.rm=TRUE)  ## Compute standard deviation between conditions
      maxMet=max(met)*(percVar[[ov]]/100)  ## Compute minimum variation allowed
      myreg=met[met>maxMet]  # Regulators to be kept
      LV.reg[[ov]]=names(met[met<=maxMet]) ## Keep names of removed regulators
      data[[ov]]=data[[ov]][names(myreg), ,drop=FALSE]
    }

  }

  if(method=="range"){

    for (ov in names(data)){

      met=apply(data[[ov]],1,function(x) max(x, na.rm=TRUE)-min(x, na.rm=TRUE) )
      maxMet=summary(met)["Max."]*(percVar[[ov]]/100)
      myreg=met[met>maxMet]
      LV.reg[[ov]]=names(met[met<=maxMet]) ## Keep names of removed regulators
      data[[ov]]=data[[ov]][names(myreg), , drop=FALSE]
    }

  }

  if(method=="IQrange"){

    for (ov in names(data)){

      met=apply(data[[ov]],1,function(x) BiocGenerics::IQR(x, na.rm=TRUE) )
      maxMet=summary(met)["Max."]*(percVar[[ov]]/100)
      myreg=met[met>maxMet]
      LV.reg[[ov]]=names(met[met<=maxMet]) ## Keep names of removed regulators
      data[[ov]]=data[[ov]][names(myreg), , drop=FALSE]
    }

  }

  if(method=="user"){

    for (ov in names(data)){

      if (min(dim(data[[ov]])) > 0) {
        met = apply(data[[ov]], 1, function(x) max(x, na.rm=TRUE)-min(x, na.rm=TRUE) )
        maxMet = percVar[[ov]] ## We don't consider the max, just the percentage defined by the user
        myreg=met[met>maxMet]
        LV.reg[[ov]]=names(met[met<=maxMet]) ## Keep names of removed regulators
        data[[ov]]=data[[ov]][names(myreg), , drop=FALSE]
      }

    }


  }

  results=vector("list", length=2)
  results[[1]]=data
  results[[2]]=LV.reg
  names(results)=c("data", "LV.reg")

  return(results)
}



# Obtain which regulators have been removed and why -----------------------
## Input:
## RetRegul.gene: Initial regulator matrix
## myregLV: regulators removed for low variability. List by omics
## myregNA: regulators removed for NA. List by omics

#' Title
#'
#' @param RetRegul.gene
#' @param myregLV
#' @param myregNA
#' @param data.omics
#'
#' @return
#' @export
#'
#' @examples
RemovedRegulators = function(RetRegul.gene, myregLV, myregNA, data.omics){
  RegulatorsValue=NULL
  RetRegul.geneNEW = NULL

  rownames(RetRegul.gene) = RetRegul.gene[,"regulator"]
  myregini = RetRegul.gene[,"regulator"]  ## In our case, one regulator cannot belong to 2 omics
  mygene = RetRegul.gene[1,"gene"]

  for(ov in unique(RetRegul.gene[,"omic"])){

    ## remove regulators not in data.omics
    regmodel = intersect(RetRegul.gene[,"regulator"], rownames(data.omics[[ov]]))

    if (length(regmodel) > 0) {
        RegulatorsValue = cbind(RegulatorsValue, t(data.omics[[ov]][regmodel, , drop=FALSE]))
        RetRegul.geneNEW = rbind(RetRegul.geneNEW, RetRegul.gene[regmodel, , drop=FALSE])
    }

    ## NA

    if(length(intersect(myregini, myregNA[[ov]]))>0){
      RetRegul.geneNEW = rbind(RetRegul.geneNEW, data.frame(gene = mygene, regulator = intersect(myregini, myregNA[[ov]]),
                                                            omic = ov, area=RetRegul.gene[intersect(myregini, myregNA[[ov]]),"area"],
                                                            filter = "MissingValue", stringsAsFactors = FALSE))
    }

    ## LV

    if (length(intersect(myregini, myregLV[[ov]]))>0){
      RetRegul.geneNEW = rbind(RetRegul.geneNEW, data.frame(gene = mygene, regulator = intersect(myregini, myregLV[[ov]]),
                                                            omic = ov, area=RetRegul.gene[intersect(myregini, myregLV[[ov]]),"area"],
                                                            filter = "LowVariation", stringsAsFactors = FALSE))
    }
    # browser()

    if (!is.null(dim(RetRegul.geneNEW))) {
        ## RetRegul.geneNEW=as.data.frame(RetRegul.geneNEW) ## Si es una matriz con una fila,al hacer el apply me lo convierte en vector
        RetRegul.geneNEW=apply(RetRegul.geneNEW, c(1,2), as.character)
        ##RetRegul.geneNEW[,"filter"] = as.character(RetRegul.geneNEW[,"filter"])
    }
  }

  return(list("SummaryPerGene" = RetRegul.geneNEW, "RegulatorMatrix" = RegulatorsValue))

}




# Checking multi-collinearity ---------------------------------------------


#' Title
#'
#' @param data
#' @param reg.table
#' @param correlation
#' @param action
#'
#' @return
#' @export
#'
#' @examples
CollinearityFilter = function(data, reg.table, correlation = 0.8, action = c("random", "mean")) {

  ## data = Regulator data matrix for all omics where missing values and regulators with low variation have been filtered out
  #         (regulators must be in columns)
  ## reg.table = Table with "gene", "regulator", "omic", "area", filter" where omics with no regulators have been removed

  row.names(reg.table) = reg.table[,"regulator"]

  resultado = list(RegulatorMatrix = data, SummaryPerGene = reg.table)

  for (omic in unique(reg.table[,"omic"])) {

    # Initial regulators
    myreg = reg.table[which(reg.table[,"omic"] == omic), ,drop=FALSE]
    myreg = as.character(myreg[which(myreg[,"filter"] == "Model"),"regulator"])

    if (length(myreg) > 1) {  # if there is more than one regulator for this omic:

      # Checking if correlation is higher than threshold value
      mycor = data.frame(t(combn(myreg,2)), as.numeric(as.dist(cor(data[,myreg]))), stringsAsFactors = FALSE)
      mycor = mycor[mycor[,3] >= correlation,]

      if (nrow(mycor) == 1) {  ### only 2 regulators are correlated in this omic
        correlacionados = unlist(mycor[,1:2])
        resultado = CorreAction(correlacionados, data = resultado$RegulatorMatrix, reg.table = resultado$SummaryPerGene,
                                action = action, nombre = "mc1", omic = omic)
      }

      if (nrow(mycor) >= 2) {   ### more than 2 regulators might be correlated in this omic

        mygraph = graph.data.frame(mycor, directed=F)
        mycomponents = clusters(mygraph)

        for (i in 1:mycomponents$no) {
          correlacionados = names(mycomponents$membership[mycomponents$membership == i])
          resultado = CorreAction(correlacionados, data = resultado$RegulatorMatrix, reg.table = resultado$SummaryPerGene,
                                  action = action, nombre = paste("mc", i, sep = ""), omic = omic)
        }
      }

    }

  }

  rownames(resultado$SummaryPerGene) = resultado$SummaryPerGene[,"regulator"]
  return(resultado)
}




# Action to perform when having correlated regulators ---------------------

#' Title
#'
#' @param correlacionados
#' @param data
#' @param reg.table
#' @param omic
#' @param action
#' @param nombre
#'
#' @return
#' @export
#'
#' @examples
CorreAction = function (correlacionados, data, reg.table, omic, action = c("mean", "random"), nombre = "mc1") {

  regus = colnames(data)

  if (action == "random") {
    dejar = sample(correlacionados, 1)  # correlated regulator to keep
    quitar = setdiff(correlacionados,dejar) # correlated regulator to remove
    regus = setdiff(regus, quitar)  # all regulators to keep
    data = data[,regus]  # excluding filtered regulator from data
    reg.table[quitar,"filter"] = dejar    # update info in regulators table
    reg.table = apply(reg.table, 2, as.character)
    rownames(reg.table) = reg.table[,"regulator"]
  }

  if (action == "mean") {
    nuevo = rowMeans(data[,correlacionados])  # mean of the correlated regulators

    # Collapse the areas of all the regulators
    myarea = unique(reg.table[correlacionados, "area"])
    if (length(grep(";", myarea)) > 0) {
      myarea = unique(unlist(sapply(myarea, strsplit, ";")))
    }
    myarea = paste(myarea, collapse = ";")

    regus = setdiff(regus, correlacionados) # remove correlated regulators
    data = data[,regus]   # new data without correlated regulators
    data = cbind(data, nuevo)  # add new "mean" regulator
    nombre = paste(omic, nombre, sep = "_")
    colnames(data)[ncol(data)] = nombre
    reg.table[correlacionados,"filter"] = rep(nombre, length(correlacionados)) # update info in regulators table
    reg.table = rbind(reg.table,
                      data.frame("gene" = reg.table[1,1], "regulator" = nombre, "omic" = omic, "area" = myarea, "filter" = "Model",
                                 stringsAsFactors = FALSE))  # add a new row to regulators table with new regulators
    reg.table = apply(reg.table, 2, as.character)
    rownames(reg.table) = reg.table[,"regulator"]
  }

  return(list(RegulatorMatrix = data, SummaryPerGene = reg.table))

}












# Plot GLM results --------------------------------------------------------

#' Title
#'
#' @param GLMoutput
#' @param gene
#' @param regulator
#' @param reguValues
#' @param plotPerOmic
#' @param fittedModel
#' @param replicates
#' @param gene.col
#' @param regu.col
#' @param xlab
#' @param cont.var
#' @param cond2plot
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
plotGLM = function (GLMoutput, gene, regulator = NULL, reguValues = NULL, plotPerOmic = TRUE,
                    fittedModel = FALSE, replicates = FALSE, gene.col = 4, regu.col = NULL,
                    xlab = "", cont.var = "Time", cond2plot = NULL,...) {

  # Colors for omics
  omic.col = c("red3", "darkorange1", "purple2", "royalblue3", "green4")
  names(omic.col) = c("Methyl", "ChIP", "miRNA", "DNase", "TF")

  if (is.null(regu.col)) { any.col = 2 } else { any.col = regu.col }

  # Changing margin
  par(mar = c(5,3,4,3)+0.1)

  # Replicates
  if (!is.null(cont.var) && (cont.var %in% colnames(GLMoutput$arguments$edesign))) {  # we have continuous variable
    if (!is.null(cond2plot) && (cond2plot %in% colnames(GLMoutput$arguments$edesign))) { # cont.var + cond2plot

      disseny = GLMoutput$arguments$edesign[,c(cond2plot, cont.var)]
      disseny = disseny[order(disseny[,1], disseny[,2]),]
      myreplicates = apply(disseny, 1, paste, collapse = "_")
      condi1 = sort(unique(GLMoutput$arguments$edesign[,cond2plot]))
      condi2 = sort(unique(GLMoutput$arguments$edesign[,cont.var]))
      allTreatments = apply(data.frame(rep(condi1, each = length(condi2)), rep(condi2, length(condi1))), 1, paste, collapse = "_")

    } else {  # only continuous variable
      disseny = GLMoutput$arguments$edesign[,cont.var, drop = FALSE]
      disseny = disseny[order(disseny[,1]), , drop = FALSE]
      myreplicates = disseny[,1]
      condi1 = NULL
      condi2 = sort(unique(GLMoutput$arguments$edesign[,cont.var]))
      allTreatments = condi2
    }

  } else {   # no cont.var

    if (!is.null(cond2plot) && (cond2plot %in% colnames(GLMoutput$arguments$edesign))) { # only cond2plot

      disseny = GLMoutput$arguments$edesign[, cond2plot, drop = FALSE]
      disseny = disseny[order(disseny[,1]), , drop = FALSE]
      myreplicates = disseny[,1]
      condi1 = NULL
      condi2 = sort(unique(GLMoutput$arguments$edesign[,cond2plot]))
      allTreatments = condi2

    } else {  # nothing! ERROR

      stop("plotGLM() requires at least either cont.var or cond2plot")

    }
  }

  # Cast myreplicates to character
  myreplicates = as.character(myreplicates)
  names(myreplicates) = rownames(GLMoutput$arguments$edesign)

  # Error values
  # errorValues = errorValuesRegu = NULL

  getErrorValues = function(realValues, repsInfo) {
    # Disable it
    if (! replicates)
      return(NULL)

    out_values = tapply(realValues, repsInfo, function(reps) sd(reps)/sqrt(length(reps)))

    return(out_values)
  }

  if (is.null(regulator)) {  ### Plot all regulators for the given gene

    GLMgene = GLMoutput$ResultsPerGene[[gene]]

    if (is.null(GLMgene)) {
      stop(paste("No GLM was obtained for gene", gene))
    }

    if (is.null(GLMgene$significantRegulators)) { ## No significant regulators -> plot only gene

      cat("No significant regulators were found for this gene.\n")
      ######## dibujar gen !!!  TO DO

    } else {  ## Significant regulators:

      # Considering multicollinearity
      SigReg = GLMgene$allRegulators
      SigReg = SigReg[SigReg$Sig == 1, c("regulator", "omic", "area", "filter")]
      # myMC = grep("_mc", SigReg[,"filter"])
      # if (length(myMC) > 0) {
      #   losMC = SigReg[myMC,]
      #   SigReg = SigReg[-myMC,]
      # }

      SigReg = SigReg[GLMgene$significantRegulators,,drop = FALSE]

      cat(paste(nrow(SigReg), "significant regulators are to be plotted for gene", gene)); cat("\n")

      # Gene values
      geneValues = GLMgene$Y$y  ### modificar si se pinta fitted
      # myreplicates = myreplicates[rownames(GLMgene$X)]
      geneValues = tapply(geneValues, myreplicates, mean)
      # geneValues = geneValues[unique(myreplicates)]
      geneValues = geneValues[allTreatments]
      names(geneValues) = allTreatments

      errorValues = getErrorValues(GLMgene$Y$y, myreplicates)

      # X values
      x.points = 1:length(allTreatments)
      numLines = 0    ########### esto es si hay cont.var
      cont.var = GLMoutput$arguments$cont.var
      # condiciones = GLMoutput$arguments$edesign[,setdiff(colnames(GLMoutput$arguments$edesign), cont.var),drop = FALSE]
      # if (ncol(condiciones) > 0) {
      #   for (cc in 1:length(condiciones)) {
      #     numLines = numLines + length(unique(condiciones[,cc])) - 1
      #   }
      # }

      if (!is.null(condi1)) numLines = length(condi1) - 1

      # eje = aggregate(GLMoutput$arguments$edesign[,cont.var],
      #         by = list(myreplicates), unique); print(eje)
      # rownames(eje) = eje[,1]
      # eje = eje[unique(myreplicates),2]
      eje = allTreatments


      if (plotPerOmic) { ## All regulators from the same omic in the same plot

        myomics = unique(SigReg$omic)

        for (oo in myomics) {

          SigRegOmic = SigReg[SigReg$omic == oo,]

          # reguValues = GLMgene$X
          omicValues = t(GLMoutput$arguments$dataOmics[[oo]])
          omicValues = omicValues[rownames(GLMgene$X),]
          reguValues = omicValues[, colnames(omicValues) == SigRegOmic$regulator[1]]

          errorValuesRegu = getErrorValues(reguValues, myreplicates)

          myreplicates = myreplicates[rownames(GLMgene$X)]
          reguValues = tapply(reguValues, myreplicates, mean)
          reguValues = reguValues[allTreatments]
          names(reguValues) = allTreatments

          mycol = omic.col[oo]
          if (is.na(mycol)) mycol = any.col

          if (nrow(SigRegOmic) == 1) {
            leftlab = SigRegOmic$regulator[1]
          } else { leftlab = oo }

          yleftlim = range(omicValues[,SigRegOmic$regulator], na.rm = TRUE)

          if (! is.null(errorValuesRegu)) {
            yleftlim = range(
              apply(omicValues[, SigRegOmic$regulator, drop = FALSE], 2, function(x) {
                errorInd = getErrorValues(x, myreplicates)
                meanValues = tapply(x, myreplicates, mean)
                meanValues = meanValues[unique(myreplicates)]

                return(c(meanValues - errorInd, meanValues + errorInd))
              }))
          }

          plotGeneRegu(x.points = x.points, geneValues = geneValues, reguValues = reguValues,
                       col = c(gene.col, mycol), yleftlim = yleftlim,
                       xlab = xlab, condi1 = condi1, condi2 = condi2,
                       yylab = c(gene, leftlab), pch = c(16,16),
                       main = oo, numLines = numLines, x.names = eje,
                       geneErrorValues = errorValues, reguErrorValues = errorValuesRegu)

          if (nrow(SigRegOmic) > 1) {
            for (i in 2:nrow(SigRegOmic)) {

              # reguValues = GLMgene$X
              omicValues = t(GLMoutput$arguments$dataOmics[[oo]])
              omicValues = omicValues[rownames(GLMgene$X),]
              reguValues = omicValues[, colnames(omicValues) == SigRegOmic$regulator[i]]
              # reguValues = reguValues[, colnames(reguValues) == SigRegOmic$regulator[i]]

              errorValuesRegu = getErrorValues(reguValues, myreplicates)

              reguValues = tapply(reguValues, myreplicates, mean)
              reguValues = reguValues[allTreatments]
              names(reguValues) = allTreatments

              lines(x.points, reguValues, type = "o", lwd = 2, pch = i, col = mycol, lty = i)

              if (! is.null(errorValuesRegu)) {
                arrows(x.points, reguValues - errorValuesRegu, x.points, reguValues + errorValuesRegu,
                       code = 3, length = 0.02, angle = 90, col = mycol)
              }
            }
          }

        }

      } else {  ## Each regulator in a separate plot

        for (rr in SigReg$regulator) {

          # reguValues = GLMgene$X
          oo = GLMoutput$ResultsPerGene[[gene]]$allRegulators[rr,"omic"]
          omicValues = t(GLMoutput$arguments$dataOmics[[oo]])
          # omicValues = omicValues[rownames(GLMgene$X),]
          reguValues = omicValues[, colnames(omicValues) == rr]
          # reguValues = reguValues[, colnames(reguValues) == rr]

          errorValuesRegu = getErrorValues(reguValues, myreplicates)

          # myreplicates = myreplicates[rownames(GLMgene$X)]
          reguValues = tapply(reguValues, myreplicates, mean)
          # reguValues = reguValues[unique(myreplicates)]
          reguValues = reguValues[allTreatments]
          names(reguValues) = allTreatments

          mycol = omic.col[SigReg[rr, "omic"]]
          if (is.na(mycol)) mycol = any.col

          plotGeneRegu(x.points = x.points, geneValues = geneValues, reguValues = reguValues,
                       col = c(gene.col, mycol),
                       xlab = xlab,
                       yylab = c(gene, rr), pch = c(16,16),
                       main = paste(as.character(SigReg[rr, c("omic", "area")]), collapse = " "),
                       numLines = numLines, x.names = eje, condi1 = condi1, condi2 = condi2,
                       geneErrorValues = errorValues, reguErrorValues = errorValuesRegu)

        }

      }

      return(GLMgene$allRegulators[GLMgene$significantRegulators, -6])
    }

  }


  if (is.null(gene)) {  ### Plot all genes regulated by the regulator

    SigniReguGene = GetPairsGeneRegulator(genes = NULL, getGLMoutput = GLMoutput)
    SigniReguGene = SigniReguGene[SigniReguGene[,"regulator"] == regulator,]
    myomic = SigniReguGene[1,"omic"]

    if (nrow(SigniReguGene) > 0) {  # When there are genes regulated by this regulator

      if (is.null(reguValues)) {  # User does not provide reguValues
          # i=1
          reguValues = GLMoutput$arguments$dataOmics[[myomic]][regulator,]
          # reguValues = GLMoutput$PerGene[[SigniReguGene[i,"gene"]]]$X
          # reguValues = reguValues[, colnames(reguValues) == regulator]
          # if (length(reguValues) == 0) {
          #   while ((i <= nrow(SigniReguGene)) && (length(reguValues) == 0)) {
          #     i = i+1
          #     reguValues = GLMoutput$ResultsPerGene[[SigniReguGene[i,"gene"]]]$X
          #     reguValues = reguValues[, colnames(reguValues) == regulator]
          #   }
          # }
      }

      numGenes = length(SigniReguGene$gene)
      cat(paste(numGenes, "genes are regulated by", regulator)); cat("\n")

      if (length(reguValues) > 0) {  # reguValues are available (recovered or given by user)

        errorValuesRegu = getErrorValues(reguValues, myreplicates)

        # myreplicates = myreplicates[rownames(geneResults$X)] ## OJO, lo quito porque no esta GeneResults, hace falta realmente??
        reguValues = tapply(reguValues, myreplicates, mean)
        reguValues = reguValues[allTreatments]
        names(reguValues) = allTreatments

        lapply(1:numGenes, function (i) {

          geneValues = GLMoutput$ResultsPerGene[[SigniReguGene[i,"gene"]]]$Y$y  ### modificar si se pinta fitted
          geneValues = tapply(geneValues, myreplicates, mean)
          geneValues = geneValues[allTreatments]
          names(geneValues) = allTreatments

          errorValues = getErrorValues(GLMoutput$ResultsPerGene[[SigniReguGene[i,"gene"]]]$Y$y, myreplicates)

          x.points = 1:length(allTreatments)
          numLines = 0    ########### esto es si hay cont.var
          # cont.var = GLMoutput$arguments$cont.var
          # condiciones = GLMoutput$arguments$edesign[,setdiff(colnames(GLMoutput$arguments$edesign), cont.var),drop = FALSE]
          # if (ncol(condiciones) > 0) {
          #   for (cc in 1:ncol(condiciones)) {
          #     numLines = numLines + length(unique(condiciones[,cc])) - 1
          #   }
          # }
          # eje = aggregate(GLMoutput$arguments$edesign[,cont.var],
          #                 by = list(myreplicates), unique)
          # rownames(eje) = eje[,1]
          # eje = eje[unique(myreplicates),2]

          if (!is.null(condi1)) numLines = length(condi1) - 1

          eje = allTreatments

          plotGeneRegu(x.points = x.points, geneValues = geneValues, reguValues = reguValues,
                       col = c(gene.col, any.col),
                       xlab = xlab, condi1 = condi1, condi2 = condi2,
                       yylab = c(SigniReguGene[i,"gene"], regulator), pch = c(16,16),
                       main = paste(as.character(SigniReguGene[1,c("omic", "area")]), collapse = " "),
                       numLines = numLines, x.names = eje,
                       geneErrorValues = errorValues, reguErrorValues = errorValuesRegu)
        })
      } else { cat("Regulator values could not be recovered from GLMoutput. Please provide them in reguValues argument to generate the plot.\n") }

      return(SigniReguGene$gene)

    } else { cat(paste("There are no genes significantly regulated by", regulator)); cat("\n") }

  }


  if (!is.null(gene) && !is.null(regulator)) {  ### Plot only the given gene and the given regulator

    geneResults = GLMoutput$ResultsPerGene[[gene]]

    if (is.null(geneResults)) {
      stop(paste("No GLM was obtained for gene", gene))
    } else {

      if (is.null(reguValues)) {  # User does not provide reguValues
          reguValues = geneResults$allRegulators[regulator, "omic"] # omic
          reguValues = as.numeric(GLMoutput$arguments$dataOmics[[reguValues]][regulator,]) # regulator values
          # reguValues = reguValues[, colnames(reguValues) == regulator]
      }

      if (length(reguValues) > 0) {  # reguValues are available (recovered or given by user)

        errorValuesRegu = getErrorValues(reguValues, myreplicates)

        # myreplicates = myreplicates[rownames(geneResults$X)]
        reguValues = tapply(reguValues, myreplicates, mean, na.rm = TRUE)
        reguValues = reguValues[allTreatments]
        names(reguValues) = allTreatments

        geneValues = GLMoutput$ResultsPerGene[[gene]]$Y$y  ### modificar si se pinta fitted
        geneValues = tapply(geneValues, myreplicates, mean)
        geneValues = geneValues[allTreatments]
        names(geneValues) = allTreatments

        errorValues = getErrorValues(GLMoutput$ResultsPerGene[[gene]]$Y$y, myreplicates)

        x.points = 1:length(allTreatments)
        numLines = 0    ########### esto es si hay cont.var
        # cont.var = GLMoutput$arguments$cont.var
        # condiciones = GLMoutput$arguments$edesign[,setdiff(colnames(GLMoutput$arguments$edesign), cont.var),drop = FALSE]
        # if (ncol(condiciones) > 0) {
        #   for (cc in 1:length(condiciones)) {
        #     numLines = numLines + length(unique(condiciones[,cc])) - 1
        #   }
        # }
        # eje = aggregate(GLMoutput$arguments$edesign[,cont.var],
        #                 by = list(myreplicates), unique)
        # rownames(eje) = eje[,1]
        # eje = eje[unique(myreplicates),2]

        if (!is.null(condi1)) numLines = length(condi1) - 1

        eje = allTreatments

        plotGeneRegu(x.points = x.points, geneValues = geneValues, reguValues = reguValues,
                     col = c(gene.col, any.col),
                     xlab = xlab, condi1 = condi1, condi2 = condi2,
                     yylab = c(gene, regulator), pch = c(16,16),
                     main = paste(as.character(geneResults$allRegulators[regulator, c("omic", "area")]),
                                  collapse = " "),
                     numLines = numLines, x.names = eje,
                     geneErrorValues = errorValues, reguErrorValues = errorValuesRegu)

      } else {

        cat("Regulator values could not be recovered from GLMoutput.\n")

        regulator = geneResults$allRegulators[regulator,"filter"]

        if (regulator %in% rownames(geneResults$allRegulators)) {

          cat(paste(regulator, "values will be plotted instead.")); cat("\n")
          cat(paste(regulator, "summarizes information from the following correlated regulators:")); cat("\n")
          cat(geneResults$allRegulators[geneResults$allRegulators[,"filter"] == regulator,"regulator"]); cat("\n")

          reguValues = geneResults$X
          reguValues = reguValues[, colnames(reguValues) == regulator]

          errorValuesRegu = getErrorValues(reguValues, myreplicates)

          reguValues = tapply(reguValues, myreplicates, mean)
          reguValues = reguValues[allTreatments]
          names(reguValues) = allTreatments

          geneValues = GLMoutput$ResultsPerGene[[gene]]$Y$y  ### modificar si se pinta fitted
          geneValues = tapply(geneValues, myreplicates, mean)
          geneValues = geneValues[allTreatments]
          names(geneValues) = allTreatments

          errorValues = getErrorValues(GLMoutput$ResultsPerGene[[gene]]$Y$y , myreplicates)

          x.points = 1:length(unique(myreplicates))
          numLines = 0    ########### esto es si hay cont.var
          # cont.var = GLMoutput$arguments$cont.var
          # condiciones = GLMoutput$arguments$edesign[,setdiff(colnames(GLMoutput$arguments$edesign), cont.var),drop = FALSE]
          # if (ncol(condiciones) > 0) {
          #   for (cc in 1:length(condiciones)) {
          #     numLines = numLines + length(unique(condiciones[,cc])) - 1
          #   }
          # }
          # eje = aggregate(GLMoutput$arguments$edesign[,cont.var],
          #                 by = list(myreplicates), unique)
          # rownames(eje) = eje[,1]
          # eje = eje[unique(myreplicates),2]

          if (!is.null(condi1)) numLines = length(condi1) - 1

          eje = allTreatments

          plotGeneRegu(x.points = x.points, geneValues = geneValues, reguValues = reguValues,
                       col = c(gene.col, any.col),
                       xlab = xlab, condi1 = condi1, condi2 = condi2,
                       yylab = c(gene, regulator), pch = c(16,16),
                       main = paste(as.character(geneResults$allRegulators[regulator, c("omic", "area")]),
                                    collapse = " "),
                       numLines = numLines, x.names = eje,
                       geneErrorValues = errorValues, reguErrorValues = errorValuesRegu)
        } else {
          cat("The selected regulator was not declared as significant by the GLM.\n")
          cat("Please select another regulator or provide the regulator values.\n")
        }

      }

    }

  }


}





# Function to obtain all significant pairs gene-regulator per omic --------

# For all genes
#' Title
#'
#' @param genes
#' @param getGLMoutput
#' @param correSense
#'
#' @return
#' @export
#'
#' @examples
GetPairsGeneRegulator = function (genes = NULL, getGLMoutput, correSense = FALSE) {

  if (is.null(genes)) genes = rownames(getGLMoutput$GlobalSummary$ReguPerGene)

  myresults = do.call("rbind", lapply(genes, GetPairs1GeneAllReg, getGLMoutput, correSense))

  #   colnames(myresults) = c("gene", "regulator", "omic", "area")
  return(myresults)
}


# For only 1 gene
#' Title
#'
#' @param gene
#' @param getGLMoutput
#' @param correSense
#'
#' @return
#' @export
#'
#' @examples
GetPairs1GeneAllReg = function (gene, getGLMoutput, correSense = FALSE) {

  reguSignif = getGLMoutput$ResultsPerGene[[gene]]$significantRegulators

  if (is.null(reguSignif)) {  # NO significant regulators
    return (NULL)

  } else {  # Significant regulators

    reguSignif = getGLMoutput$ResultsPerGene[[gene]]$allRegulators[reguSignif,]
    reguSignif = reguSignif[,c("gene", "regulator", "omic", "area")]
    return (reguSignif)
  }
}






# Plot 1 gene versus 1 regulator ------------------------------------------

#' Title
#'
#' @param x.points
#' @param geneValues
#' @param reguValues
#' @param geneErrorValues
#' @param reguErrorValues
#' @param col
#' @param xlab
#' @param yylab
#' @param pch
#' @param main
#' @param numLines
#' @param x.names
#' @param yleftlim
#' @param yrightlim
#'
#' @return
#' @export
#'
#' @examples
plotGeneRegu = function (x.points, geneValues, reguValues, geneErrorValues, reguErrorValues, col = c(1,2),
                         xlab = "", yylab = c("right", "left"), pch = c(16,17), main = "",
                         numLines = 0, x.names = NULL, yleftlim, yrightlim,
                         condi1 = NULL, condi2) {

  # Adjust the axis to include the error value
  if (missing(yrightlim)) {
    if (! missing(geneErrorValues) && ! is.null(geneErrorValues)) {
      yrightlim = range(c(geneValues - geneErrorValues, geneValues + geneErrorValues), na.rm = TRUE)
    } else {
      yrightlim = range(geneValues, na.rm = TRUE)
    }
  }

  if (missing(yleftlim)) {
   if (! missing(reguErrorValues) && ! is.null(reguErrorValues)) {
    yleftlim = range(c(reguValues - reguErrorValues, reguValues + reguErrorValues), na.rm = TRUE)
   } else {
    yleftlim = range(reguValues, na.rm = TRUE)
   }
  }

  plot.y2(x = x.points, yright = geneValues, yleft = reguValues, yleftlim = yleftlim,
          col = col, xlab = xlab, yylab = yylab, pch = pch, main = main, yrightlim = yrightlim,
          yrightErrorValues = geneErrorValues, yleftErrorValues = reguErrorValues)

  if (numLines > 0) {
      # aqui = length(x.points) / (numLines + 1)

      if (!is.null(condi1)) {
          aqui = length(condi2)
          for (i in 1:numLines) {
              abline(v = (2*aqui + 1)/2, lty = 2, col = 1)
              aqui = aqui + length(condi2)
          }
    }
  }

  if (!is.null(x.names)) {
    axis(side=1, at = x.points, labels = x.names, cex.axis = 0.8, las=2)
  }
}








# Plot Y2 -----------------------------------------------------------------

# By Ajay Shah (taken from [R] Plot 2 time series with different y axes (left and right),
# in https://stat.ethz.ch/pipermail/r-help/2004-March/047775.html)

# Modified by: Sonia Tarazona

### PARAMETERS (default):
# x: data to be drawn on X-axis
# yright: data to be drawn on Y right axis
# yleft: data to be drawn on Y left axis
# yrightlim (range(yright, na.rm = TRUE)): ylim for rigth Y-axis
# yleftlim (range(yleft, na.rm = TRUE)): ylim for left Y-axis
# xlab (NULL): Label for X-axis
# yylab (c("","")): Labels for right and left Y-axis
# pch (c(1,2)): Type of symbol for rigth and left data
# col (c(1,2)): Color for rigth and left data
# linky (TRUE): If TRUE, points are connected by lines.
# smooth (0): Friedman's super smoothing
# lwds (1): Line width for smoothed line
# length (10): Number of tick-marks to be drawn on axis
# ...: Other graphical parameters to be added by user (such as main, font, etc.)
###


#' Title
#'
#' @param x
#' @param yright
#' @param yleft
#' @param yrightlim
#' @param yleftlim
#' @param xlim
#' @param xlab
#' @param yylab
#' @param lwd
#' @param pch
#' @param col
#' @param type
#' @param linky
#' @param smooth
#' @param bg
#' @param lwds
#' @param length
#' @param ...
#' @param x2
#' @param yright2
#' @param yleft2
#' @param col2
#' @param yrightErrorValues
#' @param yleftErrorValues
#'
#' @return
#' @export
#'
#' @examples
plot.y2 <- function(x, yright, yleft, yrightlim = range(yright, na.rm = TRUE),
                    yleftlim = range(yleft, na.rm = TRUE),
                    xlim = range(x, na.rm = TRUE),
                    xlab = NULL, yylab = c("",""), lwd = c(2,2),
                    pch = c(1,2), col = c(1,2), type = c("o","o"),
                    linky = TRUE, smooth = 0, bg = c("white","white"),
                    lwds = 1, length = 10, ...,
                    x2 = NULL, yright2 = NULL, yleft2 = NULL, col2 = c(3,4),
                    yrightErrorValues, yleftErrorValues
                    )
{
  #par(mar = c(5,2,4,2), oma = c(0,3,0,3))

  ## Plotting RIGHT axis data

  plot(x, yright, axes = FALSE, ylab = "", xlab = xlab, ylim = yrightlim,
       xlim = xlim, pch = pch[1], type = type[1], lwd = lwd[1],
       col = col[1], ...)

  axis(4, pretty(yrightlim, length), col = 1, col.axis = 1)

  if (is.null(yright2) == FALSE) {
    points(x2, yright2, type = type[1], pch = pch[1], lwd = lwd[1], col = col2[1], ...)
  }

  #if (linky) lines(x, yright, col = col[1], ...)

  if (smooth != 0) lines(supsmu(x, yright, span = smooth), col = col[1], lwd = lwds, ...)

  if(yylab[1]=="") {
    mtext(deparse(substitute(yright)), side = 4, outer = FALSE, line = 2,
          col = col[1], cex = 0.9,...)
  } else {
    mtext(yylab[1], side = 4, outer = FALSE, line = 2, col = col[1], cex = 0.9,...)
  }

  # Plot arrows showing standard error
  if (!missing(yrightErrorValues) && ! is.null(yrightErrorValues)) {
    arrows(x, yright - yrightErrorValues, x, yright + yrightErrorValues,
           code = 3, length = 0.02, angle = 90, col = col[1])
  }


  par(new = T)

  ## Plotting LEFT axis data
  plot(x, yleft, axes = FALSE, ylab = "" , xlab = xlab, ylim = yleftlim,
       xlim = xlim, bg = bg[1],
       pch = pch[2], type = type[2], lwd = lwd[2], col = col[2], ...)

  box()

  axis(2, pretty(yleftlim, length), col = 1, col.axis = 1)

  if (is.null(yleft2) == FALSE) {
    points(x2, yleft2, type = type[2], pch = pch[2], bg = bg[2],
           lwd = lwd[2], col = col2[2], ...)
  }


  #if (linky) lines(x, yleft, col = col[2], ...)

  if (smooth != 0) lines(supsmu(x, yleft, span = smooth), col = col[2], lwd=lwds, ...)

  if(yylab[2] == "") {
    mtext(deparse(substitute(yleft)), side = 2, outer = FALSE, line = 2, col = col[2], cex = 0.9, ...)
  } else {
    mtext(yylab[2], side = 2, outer = FALSE, line = 2, col = col[2], cex = 0.9, ...)
  }

  if (!missing(yleftErrorValues) && ! is.null(yleftErrorValues)) {
    arrows(x, yleft - yleftErrorValues, x, yleft + yleftErrorValues,
           code = 3, length = 0.02, angle = 90, col = col[2])
  }

  ## X-axis
  ##  axis(1, at = pretty(xlim, length))  ## Comment last line


}
