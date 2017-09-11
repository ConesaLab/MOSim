#' @import igraph MASS glmnet
NULL

##################################
## ComputeGLM:
## Obtain final GLM and Summary with main results
## By Monica: 16/12/15
## Modified by Sonia: 7/07/16
###################################



## INPUT
# - matrix.temp: The final matrix for glm
# - Res.df: residuals degree of freedom. By default:3
# - alfa: significance level used for variable selection in the stepwise regression. By default: 0.05
# - stepwise: Can be either "backward", "forward", "two.ways.backward" or "two.ways.forward"
# - edesign: edesign (necessary for the summary results)
## epsilon: the same as fixed internally

## OUTPUT
## GLMfinal: Final GLM obtained after applying different stepwise
## SummaryStepwise: List with 3 elements "p-value", "R-squared", "variables"


#' Title
#'
#' @param matrix.temp
#' @param alfa
#' @param stepwise
#' @param Res.df
#' @param family
#' @param epsilon
#' @param MT.adjust
#' @param iter.max
#'
#' @return
#' @export
#'
#' @examples
ComputeGLM = function(matrix.temp, alfa = 0.05, stepwise = "two.ways.backward",
                      Res.df = 4, family = negative.binomial(theta=10), epsilon = 0.00001,
                      MT.adjust = "fdr", iter.max = 100) {

  dim.reg = ncol(matrix.temp)-1
  df.total = nrow(matrix.temp)-(Res.df+1)

  y = matrix.temp[,1]
  x = matrix.temp[,-1, drop = FALSE]

  ## NOT ENOUGH DF: VARIABLE SELECTION

  if (df.total< dim.reg){

    if(stepwise=="backward" | stepwise=="two.ways.backward"){

      glmPlot = PreviousBackward(y=y, d=x, alfa=alfa, family=family, epsilon=epsilon, stepwise=stepwise,
                                 gdl.max = df.total, MT.adjust = MT.adjust, iter.max = iter.max)   ## multiple testing method!!!

    } else {

      glmPlot = StepwiseProcedure(stepwise=stepwise, y=y, x=x, family=family, alfa=alfa, epsilon=epsilon,
                                  MT.adjust = MT.adjust, iter.max = iter.max)    ## multiple testing method!!!

    }

    ## ENOUGH DF: Apply the stepwise procedure choosen by the user if it is necessary

  } else {
    glm1 = glm(y ~ ., data = x, family = family, epsilon=epsilon)
    resTemp=summary(glm1)

    mypvals = resTemp$coefficients[,4]
    mypvals = p.adjust(mypvals, method = MT.adjust)

    if(all(mypvals <= alfa)) {   ## multiple testing method!!!

      glmPlot=glm1

    } else {

      glmPlot = StepwiseProcedure(stepwise=stepwise, y=y, x=x, family=family, alfa=alfa, epsilon=epsilon,
                                  MT.adjust = MT.adjust, iter.max = iter.max)           ## multiple testing method!!!

    }
  }

  SummaryStepwise = ResultsTable(glm=glmPlot, family=family, epsilon=epsilon)


  return(list("GLMfinal" = glmPlot, "SummaryStepwise" = SummaryStepwise))

}



## PreviousBackward

# When we don't have enough df and we want to apply a backward methodology,
# it is necessary to make a previous selection

# INPUT
## y: y values (response)
## d: x values
## family: By default is negative.binomial(theta=10)
## alfa
## epsilon
## gdl.max: degree of freedom maximum allowed

## OUTPUT
## Final GLM

#' Title
#'
#' @param y
#' @param d
#' @param alfa
#' @param family
#' @param epsilon
#' @param stepwise
#' @param gdl.max
#' @param MT.adjust
#' @param iter.max
#'
#' @return
#' @export
#'
#' @examples
PreviousBackward = function(y, d, alfa, family, epsilon, stepwise, gdl.max, MT.adjust = "fdr", iter.max = 100) {
  pval <- NULL
  noms = NULL
  design <- NULL
  j = 1
  d=as.data.frame(d, check.names = FALSE)

  for (i in 1:ncol(d)) {  # one model per variable
    sub <- d[, i,drop=FALSE]
    sub <- as.data.frame(sub, check.names = FALSE)
    glm1 <- glm(y ~ ., data = sub, family = family, epsilon=epsilon)
    result <- summary(glm1)$coefficients[, 4]
    pval[i] <- result[j + 1]
    noms[i] = names(result)[j + 1]
  }

  noms=gsub("\`","",noms )  ## Remove quotation mark
  names(pval)=noms
  pval = p.adjust(pval, method = MT.adjust) ## multiple testing method!!!
  pval=sort(pval,decreasing=FALSE) ## with sorting NAs are lost --> at this step pval has not the same variables than d
  pval.ini=pval[1:gdl.max] ## p-values of variables to initially enter the model
  d.ini=d[,names(pval.ini)]  ## Initial data
  d.res=d[,names(pval[(gdl.max+1):length(pval)])] ## The rest of the var. We want to maintain the order
  varNA = setdiff(colnames(d), c(colnames(d.ini), colnames(d.res))) # variables with pvalue=NA (singularities in the GLM)
  if (length(varNA) > 0) d.res = cbind(d.res, d[,varNA])

  tmax=ncol(d)-gdl.max  # number of variables that must be OUT of the model (only for "backward")

  if(stepwise=="backward"){

      glm2=stepback.gdl(y = y, d.ini = d.ini, d.res = d.res, alfa = alfa, family = family, epsilon = epsilon,
                        gdl.max = gdl.max, tmax = tmax, MT.adjust = MT.adjust)

  }

  if(stepwise=="two.ways.backward"){

    glm2=two.ways.stepback.gdl(y = y, d.ini = d.ini, d.res = d.res, alfa = alfa, family = family, epsilon = epsilon,
                               gdl.max = gdl.max, MT.adjust = MT.adjust, iter.max = 100)

  }

  return(glm2)

}






### stepback.gdl

# INPUT
## y: y values (response)
## d.ini: the x variables with lesser p-value (the number is determined by the degrees of freedom)
## d.res: the rest of variables (to complete the initial number)
## family: By default is negative.binomial(theta=10)
## alfa
## epsilon
## gdl.max
## tmax: Obtained previously (tmax=ncol(d)-gdl.max)

## OUTPUT
## GLM after applying the stepwise methodology

#' Title
#'
#' @param y
#' @param d.ini
#' @param d.res
#' @param alfa
#' @param family
#' @param epsilon
#' @param gdl.max
#' @param tmax
#' @param MT.adjust
#'
#' @return
#' @export
#'
#' @examples
stepback.gdl = function (y, d.ini, d.res, alfa, family, epsilon, gdl.max, tmax, MT.adjust = "fdr"){

  # d.orig = cbind(d.ini, d.res)

  t=1 ## We introduce the var one per one, lesser p-value

  d.ini=as.data.frame(d.ini, check.names = FALSE)  # variables with significant p-value in individual GLMs
  lm1 <- glm(y ~ ., data = d.ini, family=family, epsilon=epsilon)
  result <- summary(lm1)$coefficients[,4]  ## p-valores
  result = p.adjust(result, method = MT.adjust)  ## multiple testing method
  names(result)=gsub("\`","",names(result))  ## Remove quotation mark
  if (length(result[-1]) < ncol(d.ini)) { d.res = cbind(d.res, d.ini[,setdiff(colnames(d.ini), names(result)[-1])]) }
  d.ini=d.ini[,names(result)[-1], drop = FALSE]

  max <- max(result[-1], na.rm = TRUE)

  while (max > alfa) {   # there are variables with non-significant p-value in the global GLM
    varout <- names(result[-1])[result[-1] == max][1]
    varout = gsub("\`","",varout)
    pos <- findPosition(matrix = d.ini, vari = varout)
    d.ini <- d.ini[, -pos, drop = FALSE]

    while (t < tmax) {
      d.ini=cbind(d.ini, d.res[,t, drop=FALSE])
      t=t+1

      lm1 <- glm(y ~ ., data = d.ini, family=family, epsilon=epsilon)

      result <- summary(lm1)$coefficients[,4]  ## p-valores
      result = p.adjust(result, method = MT.adjust)   ## multiple testing method
      names(result)=gsub("\`","",names(result))
      max <- max(result[-1], na.rm = TRUE)

      if(max>alfa){
        varout <- names(result[-1])[result[-1] == max][1]
        varout=gsub("\`","",varout)
        pos <- findPosition(matrix = d.ini, vari = varout)
        d.ini <- d.ini[, -pos, drop = FALSE]

      } else {
        t=tmax
      }
    }

    lm1 <- glm(y ~ ., data = d.ini, family=family, epsilon=epsilon)
    result <- summary(lm1)$coefficients[,4]   ## p-valores
    result = p.adjust(result, method = MT.adjust) ## multiple testing method
    names(result)=gsub("\`","",names(result))
    max <- max(result[-1], na.rm = TRUE)
    if (length(result[-1]) == 1) {
      max <- result[-1]
      if (max > alfa) {
        max = 0
        lm1 <- glm(y ~ 1,  family=family, epsilon=epsilon)
      }
    }
  }
  return(lm1)
}




### TWO WAYS STEPBACK- NOT ENOUGH GDL

# INPUT
## y: y values (response)
## d.ini: the x variables with lesser p-value (the number is determined by the degrees of freedom)
## d.res: the rest of variables (to complete the initial number)
## gdl.max: degree of freedom max
## family: By default is negative.binomial(theta=10)
## alfa
## epsilon

## OUTPUT
## GLM after applying the stepwise methodology

#' Title
#'
#' @param y
#' @param d.ini
#' @param d.res
#' @param alfa
#' @param family
#' @param epsilon
#' @param gdl.max
#' @param MT.adjust
#' @param iter.max
#'
#' @return
#' @export
#'
#' @examples
two.ways.stepback.gdl = function (y, d.ini, d.res, alfa , family, epsilon, gdl.max, MT.adjust = "fdr", iter.max = 100) {

  d.ini=as.data.frame(d.ini, check.names = FALSE)  # Data for variables initially in the model
  OUT=d.res  # Data for the variables previously discarded

  lm1 <- glm(y ~ ., data = d.ini, family = family, epsilon=epsilon)
  result <- summary(lm1)$coefficients[, 4]
  result = p.adjust(result, method = MT.adjust)   ## multiple testing method
  names(result)=gsub("\`","",names(result))  ## Remove quotation mark
  max <- max(result[-1], na.rm = TRUE)
  if (length(result[-1]) < ncol(d.ini)) { OUT = data.frame(OUT,
                                                           d.ini[,setdiff(colnames(d.ini), names(result)[-1]), drop = FALSE],
                                                           check.names = FALSE) }
  d.ini <- d.ini[, names(result)[-1], drop = FALSE]  # we need to do this to have the same order and because we may have NAs in the model due to singularities

  num.iter = 1

  while ((max > alfa) && (num.iter <= iter.max)) {

    num.iter = num.iter + 1

    varout <- names(result)[result == max]
    varout=gsub("\`","",varout)
    pos <-  findPosition(matrix = d.ini, vari = varout)
    OUT <- as.data.frame(cbind(OUT, d.ini[, pos,drop=FALSE]), check.names = FALSE)

    d.ini <- d.ini[, -pos, drop = FALSE]

    j = ncol(d.ini) + 1
    pval <- NULL
    for (i in 1:ncol(OUT)) {
      sub <- cbind(d.ini, OUT[, i, drop=FALSE])
      sub <- as.data.frame(sub, check.names = FALSE)
      lm2 <- glm(y ~ ., data = sub, family = family, epsilon=epsilon)
      result <- summary(lm2)$coefficients[, 4]
      names(result)=gsub("\`","",names(result))
      pval[i] <- result[j + 1]
    }

    pval = p.adjust(pval, method = MT.adjust)  ## multiple testing method
    min <- min(pval, na.rm = TRUE)

    while ((min <= alfa) && (ncol(d.ini) < gdl.max)) {

      pos <- which(pval == min)
      pos=na.omit(pos)
      if ((ncol(d.ini) + length(pos)) > gdl.max) { set.seed(123); pos = sample(pos, size = gdl.max - ncol(d.ini)) }
      d.ini <- cbind(d.ini, OUT[, pos, drop=FALSE])
      d.ini <- as.data.frame(d.ini, check.names = FALSE)

      #       if (ncol(OUT) == 2) {
      #         max <- max(pval, na.rm = TRUE)
      #         b <- pval == max
      #         c <- c(1:length(pval))
      #         last <- c[b]
      #         lastname <- colnames(OUT)[last]
      #       }

      OUT <- OUT[, -pos, drop = FALSE]

      #       if (is.null(dim(OUT))) {
      #         OUT <- as.data.frame(OUT)
      #         colnames(OUT) <- lastname
      #       }

      j = ncol(d.ini) + 1
      pval <- NULL
      for (i in 1:ncol(OUT)) {
        sub <- cbind(d.ini, OUT[, i])
        sub <- as.data.frame(sub, check.names = FALSE)
        lm2 <- glm(y ~ ., data = sub, family = family, epsilon=epsilon)
        result <- summary(lm2)$coefficients[, 4]
        names(result)=gsub("\`","",names(result))
        pval[i] <- result[j + 1]
      }

      pval = p.adjust(pval, method = MT.adjust)  ## multiple testing method
      min <- min(pval, na.rm = TRUE)
      if (ncol(OUT) == 1) {
        if (min <= alfa) {
          d.ini <- cbind(d.ini, OUT[, 1])
          d.ini <- as.data.frame(d.ini, check.names = FALSE)
          colnames(d.ini)[j] <- colnames(OUT)[1]
        }
        min = 1
      }
    }

    if (ncol(d.ini) > 0) {
      lm1 <- glm(y ~ ., data = d.ini, family = family, epsilon=epsilon)
      result <- summary(lm1)$coefficients[, 4]
      result = p.adjust(result, method = MT.adjust) ## multiple testing method
      names(result)=gsub("\`","",names(result))
      max <- max(result[-1], na.rm = TRUE)
    } else { result = c(1,1) }

    if (length(result[-1]) == 1) {
      max <- result[-1]
      if (max > alfa) {
        max = 0
        lm1 <- glm(y ~ 1, family = family, epsilon=epsilon)
      }
    }
  }

  return(lm1)
}





#### StepwiseProcedure

## INPUT
## stepwise:1,2,3 or 4 (different stepwise methodology)
## y: y values (response)
## d: x values (predictors)
## family: By default is negative.binomial(theta=10)
## alfa
## epsilon

## OUTPUT
## GLM after applying the stepwise methodology

#' Title
#'
#' @param stepwise
#' @param y
#' @param x
#' @param family
#' @param alfa
#' @param epsilon
#' @param MT.adjust
#' @param iter.max
#'
#' @return
#' @export
#'
#' @examples
StepwiseProcedure = function(stepwise, y, x, family, alfa, epsilon, MT.adjust = "fdr", iter.max = 100){

  if(stepwise=="forward") {

    glm2 = stepforMOD(y = y, d = x, family = family, alfa = alfa, epsilon = epsilon, MT.adjust = MT.adjust)

    glmPlot=glm2

  }

  if(stepwise=="two.ways.forward") {

    glm3=two.ways.stepforMOD(y=y,d=x,family=family,alfa=alfa, epsilon=epsilon, MT.adjust = MT.adjust, iter.max = 100)

    glmPlot=glm3

  }

  if(stepwise=="two.ways.backward") {

    glm4=two.ways.stepbackMOD(y=y,d=x,family=family,alfa=alfa, epsilon=epsilon, MT.adjust = MT.adjust, iter.max = 100)

    glmPlot=glm4

  }

  if(stepwise=="backward") {

    glm5=stepbackMOD(y=y,d=x,family=family,alfa=alfa, epsilon=epsilon, MT.adjust = MT.adjust)

    glmPlot=glm5

  }

  return(glmPlot)

}





### Method: forward
# If you have enough degree of freedom or you are using "forward" methods,
# you can use the classical procedure

## INPUT
## y: y values (response)
## d: x values (predictors)
## alfa
## family: By default is negative.binomial(theta=10)
## epsilon

## OUTPUT
## GLM after applying the stepwise methodology

#' Title
#'
#' @param y
#' @param d
#' @param alfa
#' @param family
#' @param epsilon
#' @param MT.adjust
#'
#' @return
#' @export
#'
#' @examples
stepforMOD = function (y, d, alfa, family, epsilon, MT.adjust = "fdr") {
  pval <- NULL
  design <- NULL
  j = 1
  d=as.data.frame(d, check.names = FALSE)

  for (i in 1:ncol(d)) {
    sub <- cbind(design, d[, i])
    sub <- as.data.frame(sub, check.names = FALSE)
    lm2 <- glm(y ~ ., data = sub, family = family, epsilon=epsilon)
    result <- summary(lm2)
    pval[i] <- result$coefficients[, 4][j + 1]
  }
  ## multiple testing method!!!
  pval = p.adjust(pval, method = MT.adjust)
  min <- min(pval,na.rm=TRUE)
  while (min < alfa) {
    b <- pval == min
    c <- c(1:length(pval))
    pos <- c[b]
    pos <- pos[!is.na(pos)][1]
    design <- cbind(design, d[, pos])
    design <- as.data.frame(design, check.names = FALSE)
    colnames(design)[j] <- colnames(d)[pos]
    j = j + 1

    if (ncol(d) == 2) { lastname <- colnames(d)[!b] }
    d <- as.data.frame(d[, -pos], check.names = FALSE)
    if(ncol(d) == 1) {colnames(d) = lastname}

    pval <- NULL
    if (ncol(d) != 0) {
      for (i in 1:ncol(d)) {
        sub <- cbind(design, d[, i])
        sub <- as.data.frame(sub, check.names = FALSE)
        lm2 <- glm(y ~ ., data = sub, family = family, epsilon=epsilon)
        result <- summary(lm2)
        pval[i] <- result$coefficients[, 4][j + 1] ## la var que ya tenía metida no se va
      }

      pval = p.adjust(pval, method = MT.adjust)  ## multiple testing method
      min <- min(pval, na.rm = TRUE)
    }
    else min <- 1 ## Asi paro el bucle
  }
  if (is.null(design)) { ## No tengo nada. Hago el modelo solo para el intercept
    lm1 <- glm(y ~ 1, family = family, epsilon=epsilon)
  }
  else {
    lm1 <- glm(y ~ ., data = design, family = family, epsilon=epsilon)
  }
  return(lm1)
}







#### Method: Backward

## INPUT
## y: y values (response)
## d: x values (predictors)
## alfa
## family: By default is negative.binomial(theta=10)
## epsilon

## OUTPUT
## GLM after applying the stepwise methodology


#' Title
#'
#' @param y
#' @param d
#' @param alfa
#' @param family
#' @param epsilon
#' @param MT.adjust
#'
#' @return
#' @export
#'
#' @examples
stepbackMOD = function (y, d, alfa, family, epsilon, MT.adjust = "fdr"){

  d=as.data.frame(d, check.names = FALSE)
  lm1 <- glm(y ~ ., data = d, family=family, epsilon=epsilon)
  result <- summary(lm1)$coefficients[,4]
  result = p.adjust(result, method = MT.adjust)  ## multiple testing method
  names(result)=gsub("\`","",names(result))  ## Remove quotation mark
  max <- max(result[-1], na.rm = TRUE)
  while (max > alfa) {
    varout <- names(result)[result == max][1]
    varout=gsub("\`","",varout)
    pos <- findPosition(matrix = d, vari = varout)
    d <- d[, -pos]
    if (length(result[-1]) == 2) {
      min <- min(result[-1], na.rm = TRUE)
      lastname <- names(result[-1])[result[-1] == min][1]
    }
    if (is.null(dim(d))) {
      d <- as.data.frame(d, check.names = FALSE)
      colnames(d) <- lastname
    }

    if (ncol(d) > 0) {    ### IF included by STC
      lm1 <- glm(y ~ ., data = d, family = family, epsilon=epsilon)
      result <- summary(lm1)$coefficients[, 4]
      result = p.adjust(result, method = MT.adjust) ## multiple testing method
      names(result)=gsub("\`","",names(result))
      max <- max(result[-1], na.rm = TRUE)
    } else { result = c(1,1) }

    if (length(result[-1]) == 1) {
      max <- result[-1]
      if (max > alfa) {
        max = 0
        lm1 <- glm(y ~ 1,  family=family, epsilon=epsilon)
      }
    }
  }
  return(lm1)
}





### Method: two.ways.backward

# INPUT
## y: y values (response)
## d: x values (predictors)
## family: By default is negative.binomial(theta=10)
## alfa
## epsilon

## OUTPUT
## GLM after applying the stepwise methodology

#' Title
#'
#' @param y
#' @param d
#' @param alfa
#' @param family
#' @param epsilon
#' @param MT.adjust
#' @param iter.max
#'
#' @return
#' @export
#'
#' @examples
two.ways.stepbackMOD=function (y, d, alfa , family, epsilon, MT.adjust = "fdr", iter.max = 100) {

  d=as.data.frame(d, check.names = FALSE)  # explanatory variables (X)
  OUT <- NULL  # variables to go out of the model
  lm1 <- glm(y ~ ., data = d, family = family, epsilon=epsilon)  # model with ALL variables
  result <- summary(lm1)$coefficients[, 4]  # p-values
  result = p.adjust(result, method = MT.adjust)  # adjusted p-values
  names(result)=gsub("\`","",names(result))  ## Remove quotation mark
  max <- max(result[-1], na.rm = TRUE)  # highest p-value
  if (length(result[-1]) < ncol(d)) { OUT = cbind(OUT, as.matrix(d[,setdiff(colnames(d), names(result)[-1])])) }
  d <- d[, names(result)[-1], drop = FALSE]

  num.iter = 1

  while ((max > alfa) && (num.iter <= iter.max)) {   # There are non-significant variables that must go out

    num.iter = num.iter + 1

    varout <- names(result)[result == max]  # variables to go out
    varout=gsub("\`","",varout)
    pos <-  findPosition(matrix = d, vari = varout)
    saliendo = colnames(d)[pos]
    if (is.null(OUT)) {
      x = 0
      OUT <- as.data.frame(d[, pos], check.names = FALSE)   # Data of variables to go out
    } else {
      x <- ncol(OUT)
      OUT <- as.data.frame(cbind(OUT, d[, pos]), check.names = FALSE)   # Data of variables to go out
    }

    x <- (x+1):(x+length(pos))
    colnames(OUT)[x] <- saliendo

    d <- d[, -pos, drop = FALSE]  # New X without  less significant va

    j = ncol(d) + 1  # Now we will try to include variables that went out

    pval <- NULL  # Computing p-values for OUT variables (candidates to re-enter in the model)
    for (i in 1:ncol(OUT)) {
      sub <- cbind(d, OUT[, i])
      sub <- as.data.frame(sub, check.names = FALSE)
      lm2 <- glm(y ~ ., data = sub, family = family, epsilon=epsilon)
      result <- summary(lm2)$coefficients[, 4]
      names(result)=gsub("\`","",names(result))
      pval[i] <- result[j + 1]   # p-value of the new candidate variable
    }
    ## multiple testing method!!!
    pval = p.adjust(pval, method = MT.adjust)

    min <- min(pval, na.rm = TRUE)  # lowest p-value

    while (min <= alfa) {  ## Possible variables to re-enter the model
      pos <- which(pval == min)
      entrando = colnames(OUT)[pos]

      if (setequal(entrando, saliendo)) {
        min = 1    ## to avoid that the same variables keep going in and out

      } else {  ## what is going to enter is different from what has just exited

        d <- cbind(d, OUT[, pos])  # new significant variables enter into the model
        d <- as.data.frame(d)
        colnames(d)[j:(j+length(pos)-1)] <- colnames(OUT)[pos]  # Adding to X the new significant variables

        OUT <- OUT[, -pos, drop = FALSE]  # Removing added variables from OUT set

        if (ncol(OUT) > 0) {
          j = ncol(d) + 1
          pval <- NULL   # p-values for variables still OUT of the model
          for (i in 1:ncol(OUT)) {
            sub <- cbind(d, OUT[, i])
            sub <- as.data.frame(sub, check.names = FALSE)
            lm2 <- glm(y ~ ., data = sub, family = family, epsilon=epsilon)
            result <- summary(lm2)
            names(result)=gsub("\`","",names(result))
            pval[i] <- result$coefficients[, 4][j + 1]
          }
          ## multiple testing method!!!
          pval = p.adjust(pval, method = MT.adjust)

          min <- min(pval, na.rm = TRUE)
          if (ncol(OUT) == 1) {
            if (min <= alfa) {
              d <- cbind(d, OUT[, 1])
              d <- as.data.frame(d, check.names = FALSE)
              colnames(d)[j] <- colnames(OUT)[1]
            }
            min = 1
          }
        } else { min = 1 }

      }

    }

    # No more OUT variables can be added to the model -> New model
    if (ncol(d) > 0) {
      lm1 <- glm(y ~ ., data = d, family = family, epsilon=epsilon)
      result <- summary(lm1)$coefficients[, 4]
      result = p.adjust(result, method = MT.adjust)  ## multiple testing method!!!
      names(result)=gsub("\`","",names(result))
      max <- max(result[-1], na.rm = TRUE)
      if (length(unique(result[-1])) == 1) {
        # max <- result[-1]
        if (max > alfa) {
          max = 0
          lm1 <- glm(y ~ 1, family = family, epsilon=epsilon)
        }
      }
    } else {
      max = 0
      lm1 <- glm(y ~ 1, family = family, epsilon=epsilon)
    }


  }

  return(lm1)
}







#### Method: two.ways.forward

## INPUT
## y: y values (response)
## d: x values (predictors)
## alfa
## family: By default is negative.binomial(theta=10)
## epsilon

## OUTPUT
## GLM after applying the stepwise methodology

#' Title
#'
#' @param y
#' @param d
#' @param alfa
#' @param family
#' @param epsilon
#' @param MT.adjust
#' @param iter.max
#'
#' @return
#' @export
#'
#' @examples
two.ways.stepforMOD = function (y, d, alfa, family, epsilon, MT.adjust = "fdr", iter.max = 100){

  d=as.data.frame(d, check.names = FALSE)
  pval <- NULL
  design <- NULL
  j = 1

  for (i in 1:ncol(d)) {
    sub <- cbind(design, d[, i])
    sub <- as.data.frame(sub, check.names = FALSE)
    lm2 <- glm(y ~ ., data = sub, family = family, epsilon=epsilon)
    result <- summary(lm2)
    pval[i] <- result$coefficients[, 4][j + 1]
  }

  pval = p.adjust(pval, method = MT.adjust)  ## multiple testing method
  min <- min(pval, na.rm = TRUE)  ## Aquí está mi punto de partida

  num.iter = 1

  while ((min <= alfa) && (num.iter <= iter.max)) {

    num.iter = num.iter + 1

    b <- pval == min
    c <- c(1:length(pval))
    pos <- c[b]
    pos <- pos[!is.na(pos)][1]
    design <- cbind(design, d[, pos])
    design <- as.data.frame(design, check.names = FALSE)
    colnames(design)[j] <- colnames(d)[pos]

    d <- d[, -pos, drop = FALSE]

    result2 <- summary(glm(y ~ ., data = design, family = family, epsilon=epsilon))$coefficients[,4]
    result2 = p.adjust(result2, method = MT.adjust)  ## multiple testing method
    max <- max(result2[-1], na.rm = TRUE)

    while (max > alfa) {
      varout <- names(result2)[result2 == max]
      varout=gsub("\`","",varout)
      pos <-  findPosition(matrix = design, vari = varout)
      d <- as.data.frame(cbind(d, design[, pos]), check.names = FALSE)
      x <- ncol(d)
      colnames(d)[x] <- colnames(design)[pos]
      if (ncol(design) == 2) {
        min <- min(result2[-1], na.rm = TRUE)
        lastname <- names(result2)[result2 == min][1]
      }
      design <- design[, -pos]
      if (is.null(dim(design))) {
        design <- as.data.frame(design, check.names = FALSE)
        colnames(design) <- lastname
      }
      result2 <- summary(glm(y ~ ., data = design, family = family, epsilon=epsilon))$coefficients[,4]
      ## multiple testing method!!!
      result2 = p.adjust(result2, method = MT.adjust)
      max <- max(result2[-1], na.rm = TRUE)
    }
    j = ncol(design) + 1
    pval <- NULL
    for (i in 1:ncol(d)) {
      sub <- cbind(design, d[, i])
      sub <- as.data.frame(sub, check.names = FALSE)
      lm2 <- glm(y ~ ., data = sub , family = family, epsilon=epsilon)
      result <- summary(lm2)
      pval[i] <- result$coefficients[, 4][j + 1]
    }
    ## multiple testing method!!!
    pval = p.adjust(pval, method = MT.adjust)
    min <- min(pval, na.rm = TRUE)
    if (ncol(d) == 1) {
      if (min <= alfa) {
        design <- cbind(design, d[, 1])
        design <- as.data.frame(design, check.names = FALSE)
        colnames(design)[j] <- colnames(d)[1]
      }
      min = 1
    }
  }

  if (is.null(design)) {
    lm1 <- glm(y ~ 1, family = family, epsilon=epsilon)

  } else {
    lm1 <- glm(y ~ ., data = design, family = family, epsilon=epsilon)
  }
  return(lm1)
}







# FUNCTION FOR OBTAINING A SUMMARY ----------------------------------------

## INPUT
## glm: final glm obtained in previous steps
## alfa
## family
## epsilon
## edesign

## OUTPUT
## list: "p-value", "R-squared", "variables"

#' Title
#'
#' @param glm
#' @param family
#' @param epsilon
#'
#' @return
#' @export
#'
#' @examples
ResultsTable = function(glm, family, epsilon){

  y=glm$y

  # Computing p-values

  model.glm.0=glm(y~1, family=family, epsilon=epsilon)

  if(family$family=="gaussian")
  {
    test=anova(model.glm.0, glm, test="F")
    p.value = test[6][2,1]
  }
  else
  {
    test = anova(model.glm.0, glm, test="Chisq")
    p.value = test[5][2,1]
  }

  # Computing goodness of fit:

  resFin=summary(glm)
  bondad = (resFin$null.deviance-resFin$deviance)/resFin$null.deviance  ### R-SQUARED
  rownames(resFin$coefficients)=gsub("`","",rownames(resFin$coefficients)) ## Remove quotation marks
  vars.in=rownames(resFin$coefficients)   #### VARIABLES

  output=list(p.value,bondad,vars.in)
  names(output)=c("p.value","R.squared","variables")

  return(output)

}





# maSigPro::position ------------------------------------------------------

# ## old function from maSigPro
# findPosition = function (matrix, vari) {
#   a <- colnames(matrix)
#   b <- a == vari
#   c <- c(1:length(a))
#   d <- c[b]
#   return(d)
# }



#' Title
#'
#' @param matrix
#' @param vari
#'
#' @return
#' @export
#'
#' @examples
findPosition = function (matrix, vari) {
  d = NULL
  a <- colnames(matrix)
  for (i in 1:length(vari)) {
    d = c(d, which(a == vari[i]))
  }
  return(d)
}

