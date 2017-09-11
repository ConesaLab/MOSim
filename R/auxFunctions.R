#' @include ComputeGLM_function.R
NULL

#########################################
#### Auxiliar functions
#########################################

## By Sonia
## 6-Jun-2017


# Generating design matrix with interactions between experimental variables ------------------------------
## We only consider interactions of order 2

GenerateDesignMatrix = function (interactions.exp, degree, edesign, cont.var) {

  if ((cont.var %in% colnames(edesign)) && (degree > 1)) {  # polynomial terms
    politerms = paste0(cont.var, 2:degree)
    poliEquat = NULL
    for (i in 2:degree) {
      poliEquat = cbind( poliEquat, edesign[,cont.var]^i)
    }
    colnames(poliEquat) = politerms

  } else { politerms = NULL }


  if (interactions.exp) {   # we compute all possible interactions between exp variables

    if (NCOL(edesign) > 1) {
      pares = combn(colnames(edesign), 2) # all possible pairs
      pares = apply(pares, 2, function (x) paste(sprintf("`%s`", x), collapse="*")) # interactions of order 2
      fff = paste0("~ ", paste(pares, collapse = "+")) # complete formula
      fff = as.formula(fff)
      desmat = model.matrix(fff, edesign)[,-1]  # design matrix with all predictors and all interactions of order 2

      # polynomial terms?
      if (!is.null(politerms)) {
        expcond = setdiff(colnames(desmat), cont.var) # remove cont.var
        # remove interactions to avoid having int > 2
        if (length(grep(":", expcond, fixed = TRUE)) > 0) expcond = expcond[-grep(":", expcond, fixed = TRUE)]
        desmat = as.data.frame(cbind(desmat, poliEquat))
        fff = paste0("~ ", paste(sapply(politerms, function (x) paste(expcond, x, sep = ":")), collapse = "+"))
        fff = as.formula(fff)
        desmat = as.data.frame(cbind(desmat, model.matrix(fff, desmat)[,-1, drop = FALSE]))
      }


    } else {  # only 1 experimental covariate
      desmat = model.matrix(~., data = edesign)[, -1, drop = FALSE]
      if (!is.null(politerms)) {
        expcond = setdiff(colnames(desmat), cont.var) # remove cont.var
        desmat = as.data.frame(cbind(desmat, poliEquat))
        if (length(expcond) > 0) {
          fff = paste0("~ ", paste(sapply(politerms, function (x) paste(expcond, x, sep = ":")), collapse = "+"))
          fff = as.formula(fff)
          desmat = as.data.frame(cbind(desmat, model.matrix(fff, desmat)[,-1]))
        }
      }
    }

  } else {  # no interactions
    desmat = model.matrix(~., data = edesign)[, -1, drop = FALSE]
    if (!is.null(politerms)) desmat = as.data.frame(cbind(desmat, poliEquat))
  }

  return(desmat)

}



# Removing regulators with low variation ----------------------------------

LowVariationRegu = function(min.variation, data.omics, ExpGroups, associations, Allgenes) {

  if (is.null(min.variation)){

    #### Low variation cutoff is computed automatically
    data.omicsMean = vector("list", length=length(data.omics))

    for(i in 1:length(data.omics)){
      data.omicsMean[[i]]=t(apply(data.omics[[i]], 1, tapply, ExpGroups, mean))
    }
    names(data.omicsMean)=names(data.omics)

    percVar=rep(10, length(data.omicsMean)) ## We fix it
    names(percVar)=names(data.omicsMean)

    # Applying Low Variation filter
    LowVar=LowVariatFilter(data=data.omicsMean, method="sd", percVar=percVar)

    # data.omicsMean reduced: without NA and LV
    data.omicsMean=LowVar$data

    for (ov in names(associations)){
      myreg=associations[[ov]][associations[[ov]][,1] %in% Allgenes,2] # removing regulators not associated to our genes
      data.omicsMean[[ov]]=data.omicsMean[[ov]][intersect(myreg, rownames(data.omicsMean[[ov]])),]
      data.omics[[ov]]= data.omics[[ov]][rownames(data.omicsMean[[ov]]),] ## Reduced data.omics
    }


  }
  else {

    #### Low variation cutoff is set by the user

    # removing regulators not associated to our genes
    for (ov in names(associations)){
      myreg=associations[[ov]][associations[[ov]][,1] %in% Allgenes,2]
      data.omics[[ov]]=data.omics[[ov]][intersect(myreg, rownames(data.omics[[ov]])),]
    }

    # Creating vector for min.variation
    if (length(min.variation) == 1) {  ## Including min.variation = 0. I need a vector with omics names
      min.variation=rep(min.variation,length(data.omics))
      names(min.variation)=names(data.omics)
    }

    # computing mean per condition in data.omics
    data.omicsMean=vector("list", length = length(data.omics))
    for(i in 1:length(data.omics)){
      data.omicsMean[[i]] = t(apply(data.omics[[i]], 1, tapply, ExpGroups, mean))
    }
    names(data.omicsMean) = names(data.omics)

    # Applying Low Variation filter
    LowVar=LowVariatFilter(data = data.omicsMean, method = "user", percVar = min.variation)

    data.omicsMean=LowVar$data  ## data.omicsMean reduced

    ## data.omics reduced: only mygenes and without NA and LV
    for (ov in names(data.omics)){
      data.omics[[ov]] = data.omics[[ov]][rownames(data.omicsMean[[ov]]),]
    }

  }

  rm("myreg"); rm("data.omicsMean"); gc()

  # Regulators removed due to low variation filter
  myregLV=LowVar$LV.reg
  rm("LowVar"); gc()

  cat("Number of regulators with low variation:\n")
  print(sapply(myregLV, length))
  cat("\n")

  return(list("myregLV" = myregLV, "data.omics" = data.omics, "associations" = associations))
}





# Low variation filtering -------------------------------------------------
## data: For computing variability. It could be mean values between replicates
## method: One of "sd","range", "IQrange" or "user"
## percVar: percentage of variation defined by the user

LowVariatFilter=function(data, method, percVar){

  SummaryRes = LV.reg = vector("list", length=length(data))
  names(SummaryRes) = names(LV.reg) = names(data)

  if (method=="sd") {

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







# Adding interations with regulatores -------------------------------------

RegulatorsInteractions = function (interactions.reg, reguValues, des.mat, cont.var, GeneExpression, gene) {

  # Adding regulators
  des.mat2 = data.frame(des.mat, reguValues, check.names = FALSE)

  ### WITH INTERACTIONS with regulators

  if (interactions.reg > 0) {
    expcond = colnames(des.mat)

    if (interactions.reg == 1) {  # Max order of interaction = 2 & Interactions with cont.var not allowed
      if (length(grep(":", expcond, fixed = TRUE)))  expcond = expcond[-grep(":", expcond, fixed = TRUE)]
      if (length(grep(cont.var, expcond)))  expcond = expcond[-grep(cont.var, expcond)]
    }

    if (interactions.reg == 2) { # Max order of interaction = 2
      if (length(grep(":", expcond, fixed = TRUE)))  expcond = expcond[-grep(":", expcond, fixed = TRUE)]
    }

    fff = paste0("~ ",
                 paste(sapply(colnames(reguValues),
                              function (x) paste(expcond, sprintf("`%s`", x), sep = ":")),
                       collapse = "+"))

    fff = as.formula(fff)
    des.mat2 = cbind(des.mat2, model.matrix(fff, des.mat2)[,-1])

    des.mat2 = cbind(t(GeneExpression[gene,]), des.mat2)
    colnames(des.mat2)[1] = "response"

    colnames(des.mat2) = gsub("\`", "", colnames(des.mat2))


    ## PUEDE OCURRIR QUE AL METER LAS INTERACCIONES TODA LA COLUMNA SEA 0. HAGO UN FILTRO PREVIO
    sd.regulators = apply(des.mat2[,-1, drop = FALSE], 2, sd, na.rm=TRUE)
    regulators0 = names(sd.regulators[sd.regulators==0])
    if (length(regulators0)>0) des.mat2 = des.mat2[, setdiff(colnames(des.mat2), regulators0), drop=FALSE]


  } else  {    ### WITHOUT INTERACTIONS

    des.mat2 = cbind(t(GeneExpression[gene,]), des.mat)
    colnames(des.mat2)[1] = "response"

  }

  return(des.mat2)

}






# ElasticNet variable selection -------------------------------------------

ElasticNet = function (family2, des.mat2, epsilon, elasticnet, Res.df) {

  if (elasticnet > 0) {  # Variable selection with ElasticNet

    if (!is.null(family2)) {

      if (nrow(des.mat2) > 200) { mynfolds =  ceiling(nrow(des.mat2)/20) } else { mynfolds = nrow(des.mat2) }  # for CV

      cvEN = cv.glmnet(x = as.matrix(des.mat2[,-1]),
                       y = des.mat2[,1],
                       nfolds = mynfolds, alpha = 0.5, standardize = FALSE, thres = epsilon,
                       family = family2, grouped = FALSE)

      if (elasticnet == 1) {  #  "ad-hoc" penalization parameter
        myNZ = nrow(des.mat2) - 1 - Res.df
        rangeNZ = range(cvEN$nzero)

        if (myNZ < rangeNZ[1]) {
          # No ElasticNet selection is performed
          removedCoefs = NULL
          myS = NULL
        } else if ( myNZ > rangeNZ[2] ) {
          myNZ = rangeNZ[2] # take the max of non-zero coefficients
          myS = cvEN$lambda[cvEN$nzero == myNZ]
        } else { # myNZ is within the computed range of non-zero values
          myS = cvEN$lambda[cvEN$nzero == myNZ]
          if (length(myS) == 0) {
            tmp = abs(cvEN$nzero - myNZ)
            myS = cvEN$lambda[tmp == min(tmp)]
          }
        }
      }

      if (elasticnet == 2) {  # optimum penalization parameter
        myS = cvEN$lambda.min
      }

      if (length(myS) > 1) {  # more than 1 values for lambda
        myCVerror = cvEN$cvm[sapply(myS, function (i) which(cvEN$lambda == i))]
        myS = myS[which.min(myCVerror)]
      }

      if (length(myS) == 1) {
        mycoef = colnames(des.mat2)[which(as.logical(coef(cvEN, s = myS) != 0))] # selected coefficients
        mycoef = unique(c(colnames(des.mat2)[1], mycoef))

        removedCoefs = setdiff(colnames(des.mat2), mycoef)
        removedCoefs = gsub(".*[:?](.*)$", "\\1", removedCoefs)

        des.mat2 = des.mat2[, mycoef, drop = FALSE]
      }

    } else { removedCoefs = NULL }

  } else { removedCoefs = NULL }



  return(list("des.mat2" = des.mat2, "removedCoefs" = removedCoefs))
}

