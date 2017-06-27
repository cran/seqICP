##' Summary functions for 'seqICPnl' objects.
##'
##' @title summary function
##' @param object object of class 'seqICPnl'.
##' @param show.past 'TRUE' if lagged variables should also be shown.
##' @param ... additional arguments affecting the summary produced.
##'
##' @author Niklas Pfister and Jonas Peters
##'
##' @export

summary.seqICPnl <- function(object, show.past=TRUE, ...){
  stopifnot(inherits(object, "seqICPnl"))

  ns <- length(object$parent.set)
  cat(paste("\n Non-linear Invariant Causal Regression ", "at level ", object$alpha, sep = ""))
  if (object$modelReject){
    cat(paste("\n Model has been rejected at the chosen level ", object$alpha,
              ", that is no subset of variables leads to invariance across time. This can be for example due to  presence of \n (a) wrong function class or \n (b) hidden variables or \n (c) interventions on the target variable. \n", sep = ""))
    if (object$stopIfEmpty){
      cat("\n In this run, option 'stopIfEmpty' was set to TRUE so not all sets of variables have maybe been tested; rerun with option set to FALSE to get definite answer whether model is rejected")
    }
  }
  else {
    wh <- object$parent.set
    varnames <- sapply(1:object$n.var,function(i) paste("X",i,sep=""))
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
      pvals <- rep(NA,length(object$p.values))
      cat("\n\n Not all sets were considered since stopIfEmpty was TRUE.")
    }
    else{
      pvals <- c(object$p.values)
    }
    # print table of p-values
    symp <- symnum(pvals, corr = FALSE,
                   cutpoints = c(0,  .001,.01,.05, .1, 1),
                   symbols = c("***","**","*","."," "))
    result.table <- noquote(cbind(p.vals=format.pval(pvals), sig=symp))
    rownames(result.table) <- varnames
    colnames(result.table) <- c("p-value","")
    cat("\n")
    print(result.table)
    cat("--- \n")
    cat("Signif. codes: ")
    cat(attr(symp, "legend"))
  }
  cat("\n \n")
}
