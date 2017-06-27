##' Summary functions for 'seqICP' objects.
##'
##' @title summary function
##' @param object object of class 'seqICP'.
##' @param show.past 'TRUE' if lagged variables should also be shown.
##' @param ... additional arguments affecting the summary produced.
##'
##' @author Niklas Pfister and Jonas Peters
##'
##' @export

summary.seqICP <- function(object, show.past=TRUE, ...){
  stopifnot(inherits(object, "seqICP"))

  ns <- length(object$parent.set)
  cat(paste("\n Invariant Linear Causal Regression ", "at level ", object$alpha, sep = ""))
  if (object$modelReject){
    cat(paste("\n Model has been rejected at the chosen level ", object$alpha,
              ", that is no subset of variables leads to invariance across time. This can be for example due to  presence of \n (a) non-linearities or \n (b) hidden variables or \n (c) interventions on the target variable. \n", sep = ""))
    if (object$stopIfEmpty){
      cat("\n In this run, option 'stopIfEmpty' was set to TRUE so not all sets of variables have maybe been tested; rerun with option set to FALSE to get definite answer whether model is rejected")
    }
  }
  else {
    ## iid model
    if(object$model=="iid"){
      wh <- object$parent.set+1
      varnames <- character(object$n.var)
      varnames[1] <- "intercept"
      varnames[2:(object$n.var+1)] <- sapply(1:object$n.var,function(i) paste("X",i,sep=""))
      if (ns > 0){
        plural <- c("Variables "," show")
        if(ns == 1){
          plural <- c("Variable "," shows")
        }
        cat(paste("\n ",plural[1], paste(varnames[wh],collapse=", "), plural[2], " a significant causal effect", sep = ""))
      }
      else{
        cat("\n No variable shows a significant causal effect")
      }
      if(object$stopIfEmpty&length(object$S[[length(object$S)]])<object$n.var){
        pvals <- rep(NA,nrow(object$coefficients))
        cat("\n\n Not all sets were considered since stopIfEmpty was TRUE.")
      }
      else{
        pvals <- c(NA,object$p.values)
      }
    }
    ## ar model
    else{
      wh <- object$parent.set+1
      parent.ind <- which(sapply(object$S,function(y) identical(y,object$parent.set)))
      p <- object$test.results$p[parent.ind]
      varnames <- character(nrow(object$coefficients))
      varnames[1] <- "intercept"
      varnames[2:(object$n.var+1)] <- sapply(1:object$n.var,function(i) paste("X",i,"[t]",sep=""))
      if(p>0){
        for(k in 1:p){
          varnames[(object$n.var+(k-1)*(object$n.var+1)+2)] <- paste("Y0[t-",k,"]",sep="")
          varnames[(object$n.var+(k-1)*(object$n.var+1)+3):(object$n.var+1+k*(object$n.var+1))] <- sapply(1:object$n.var,function(i) paste("X",i,"[t-",k,"]",sep=""))
        }
      }
      if (ns > 0){
        plural <- c("Variables "," show")
        if(ns == 1){
          plural <- c("Variable "," shows")
        }
        cat(paste("\n ",plural[1], paste(varnames[wh],collapse=", "), plural[2], " a significant causal effect", sep = ""))
      }
      else{
        cat("\n No variable shows a significant causal effect")
      }
      if(object$stopIfEmpty&length(object$S[[length(object$S)]])<object$n.var){
        pvals <- rep(nrow(object$coefficients))
        cat("\n\n Not all sets were considered since stopIfEmpty was TRUE.")
      }
      else{
        pvals <- c(NA,object$p.values,rep(NA,(object$n.var+1)*p))
      }
    }
    # print table of coefficients and p-values
    result.table <- cbind(object$coefficients, pvals)
    rownames(result.table) <- varnames
    colnames(result.table) <- c("coefficient", "lower bound", "upper bound", " p-value")
    if(!show.past){
      result.table <- result.table[1:(length(object$p.values)+1),]
    }
    cat("\n \n ")
    printCoefmat(result.table, digits = 3, signif.stars = TRUE, 
                 P.values = TRUE, has.Pvalue = TRUE, tst.ind = 1, 
                 zap.ind = 2, eps.Pvalue = 10^(-9))
  }
  cat("\n")
  if(!object$pknown & object$model=="ar"){
    cat("Confidence intervals are only approximate due to post selection (pknown is FALSE).\n")
  }
  cat("\n")
}
