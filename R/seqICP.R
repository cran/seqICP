##' Estimates the causal parents S of the target variable Y using invariant causal prediction and fits a linear model of the form \cr Y = a X^S + N.
##'
##' The function can be applied to two types of models\cr
##' (1) a linear model (model="iid")\cr Y_i = a X_i^S + N_i\cr with iid noise N_i and \cr
##' (2) a linear autoregressive model (model="ar")\cr Y_t = a_0 X_t^S + ... + a_p (Y_(t-p),X_(t-p)) + N_t\cr with iid noise N_t.
##'
##' For both models the invariant prediction procedure is applied
##' using the hypothesis test specified by the \code{test} parameter
##' to determine whether a candidate model is invariant. For further
##' details see the references.
##' @title Sequential Invariant Causal Prediction
##' @param X matrix of predictor variables. Each column corresponds to
##'   one predictor variable.
##' @param Y vector of target variable, with length(Y)=nrow(X).
##' @param test string specifying the hypothesis test used to test for
##'   invariance of a parent set S (i.e. the null hypothesis
##'   H0_S). The following tests are available: "decoupled",
##'   "combined", "trend", "variance", "block.mean", "block.variance",
##'   "block.decoupled", "smooth.mean", "smooth.variance",
##'   "smooth.decoupled" and "hsic".
##' @param par.test parameters specifying hypothesis test. The
##'   following parameters are available: \code{grid},
##'   \code{complements}, \code{link}, \code{alpha}, \code{B} and
##'   \code{permutation}. The parameter \code{grid} is an increasing
##'   vector of gridpoints used to construct enviornments for change
##'   point based tests. If the parameter \code{complements} is 'TRUE'
##'   each environment is compared against its complement if it is
##'   'FALSE' all environments are compared pairwise. The parameter
##'   \code{link} specifies how to compare the pairwise test
##'   statistics, generally this is either max or sum. The parameter
##'   \code{alpha} is a numeric value in (0,1) indicting the
##'   significance level of the hypothesis test. The parameter
##'   \code{B} is an integer and specifies the number of Monte-Carlo
##'   samples (or permutations) used in the approximation of the null
##'   distribution. If the parameter \code{permutation} is 'TRUE' a
##'   permuatation based approach is used to approximate the null
##'   distribution, if it is 'FALSE' the scaled residuals approach is
##'   used.
##' @param model string specifying the underlying model class. Either
##'   "iid" if Y consists of independent observations or "ar" if Y has
##'   a linear time dependence structure.
##' @param par.model parameters specifying model. The following
##'   parameters are available: \code{pknown}, \code{p} and
##'   \code{max.p}. If \code{pknown} is 'FALSE' the number of lags will be
##'   determined by comparing all fits up to \code{max.p} lags using
##'   the AIC criterion. If \code{pknown} is 'TRUE' the procedure will fit
##'   \code{p} lags.
##' @param max.parents integer specifying the maximum size for
##'   admissible parents. Reducing this below the number of predictor
##'   variables saves computational time but means that the confidence
##'   intervals lose their coverage property.
##' @param stopIfEmpty if ‘TRUE’, the procedure will stop computing
##'   confidence intervals if the empty set has been accepted (and
##'   hence no variable can have a signicificant causal
##'   effect). Setting to ‘TRUE’ will save computational time in these
##'   cases, but means that the confidence intervals lose their
##'   coverage properties for values different to 0.
##' @param silent If 'FALSE', the procedure will output progress
##'   notifications consisting of the currently computed set S
##'   together with the p-value resulting from the null hypothesis
##'   H0_S
##' 
##' @return object of class 'seqICP' consisting of the following
##'   elements
##'
##' \item{parent.set}{vector of the estimated causal parents.}
##' 
##' \item{test.results}{matrix containing the result from each
##' individual test as rows.}
##'
##' \item{S}{list of all the sets that were
##' tested. The position within the list corresponds to the index in
##' the first column of the test.results matrix.}
##'
##' \item{p.values}{p-value for being not included in the set of true
##' causal parents. (If a p-value is smaller than alpha, the
##' corresponding variable is a member of parent.set.)}
##'
##' \item{coefficients}{vector of coefficients resulting from a
##' regression based on the estimated parent set.}
##'
##' \item{stopIfEmpty}{a boolean value indicating whether computations
##' stop as soon as intersection of accepted sets is empty.}
##' 
##' \item{modelReject}{a boolean value indicating if the whole model was
##' rejected (the p-value of the best fitting model is too low).}
##'
##' \item{pknown}{a boolean value indicating whether the number of
##' lags in the model was known. Only relevant if model was set to
##' "ar".}
##'
##' \item{alpha}{significance level at which the hypothesis tests were performed.}
##'
##' \item{n.var}{number of predictor variables.}
##'
##' \item{model}{either "iid" or "ar" depending on which model was selected.}
##' 
##' @export
##'
##' @import stats utils
##'
##' @author Niklas Pfister and Jonas Peters
##'
##' @references
##' Pfister, N., P. Bühlmann and J. Peters (2017).
##' Invariant Causal Prediction for Sequential Data. ArXiv e-prints (1706.08058).
##'
##' Peters, J., P. Bühlmann, and N. Meinshausen (2016).
##' Causal inference using invariant prediction: identification and confidence intervals.
##' Journal of the Royal Statistical Society, Series B (with discussion) 78 (5), 947–1012.
##'
##' @seealso The function \code{\link{seqICP.s}} allows to perform
##'   hypothesis test for individual sets S. For non-linear
##'   models the functions \code{\link{seqICPnl}} and
##'   \code{\link{seqICPnl.s}} can be used.
##'
##' @examples
##' set.seed(1)
##' 
##' # environment 1
##' na <- 140
##' X1a <- 0.3*rnorm(na)
##' X3a <- X1a + 0.2*rnorm(na)
##' Ya <- -.7*X1a + .6*X3a + 0.1*rnorm(na)
##' X2a <- -0.5*Ya + 0.5*X3a + 0.1*rnorm(na)
##' 
##' # environment 2
##' nb <- 80
##' X1b <- 0.3*rnorm(nb)
##' X3b <- 0.5*rnorm(nb)
##' Yb <- -.7*X1b + .6*X3b + 0.1*rnorm(nb)
##' X2b <- -0.5*Yb + 0.5*X3b + 0.1*rnorm(nb)
##' 
##' # combine environments
##' X1 <- c(X1a,X1b)
##' X2 <- c(X2a,X2b)
##' X3 <- c(X3a,X3b)
##' Y <- c(Ya,Yb)
##' Xmatrix <- cbind(X1, X2, X3)
##' 
##' # Y follows the same structural assignment in both environments
##' # a and b (cf. the lines Ya <- ... and Yb <- ...).
##' # The direct causes of Y are X1 and X3.
##' # A linear model considers X1, X2 and X3 as significant.
##' # All these variables are helpful for the prediction of Y.
##' summary(lm(Y~Xmatrix))
##' 
##' # apply seqICP to the same setting
##' seqICP.result <- seqICP(X = Xmatrix, Y,
##' par.test = list(grid = seq(0, na + nb, (na + nb)/10), complements = FALSE, link = sum,
##' alpha = 0.05, B =100), max.parents = 4, stopIfEmpty=FALSE, silent=FALSE)
##' summary(seqICP.result)
##' # seqICP is able to infer that X1 and X3 are causes of Y


seqICP <- function(X, Y,
                 test="decoupled",
                 par.test=list(grid=c(0,round(nrow(X)/2),nrow(X)),
                               complements=FALSE,
                               link=sum,
                               alpha=0.05,
                               B=100,
                               permutation=FALSE),
                 model="iid",
                 par.model=list(pknown=FALSE, p=0, max.p=10),
                 max.parents=ncol(X),
                 stopIfEmpty=TRUE,
                 silent=TRUE)
{

  # preprocess X and Y
  X <- unname(as.matrix(X))
  Y <- unname(as.matrix(Y))
  

  # add defaults if missing in parameter lists
  if(missing(max.parents)){
    max.parents <- ncol(X)
  }
  if(missing(stopIfEmpty)){
    stopIfEmpty <- TRUE
  }
  if(missing(silent)){
    silent <- TRUE
  }
  if(!exists("grid",par.test)){
    par.test$grid <- c(0,round(nrow(X)/2),nrow(X))
  }
  if(!exists("complements",par.test)){
    par.test$complements <- FALSE
  }
  if(!exists("link",par.test)){
    par.test$link <- sum
  }
  if(!exists("alpha",par.test)){
    par.test$alpha <- 0.05
  }
  if(!exists("B",par.test)){
    par.test$B <- 100
  }
  if(!exists("permutation",par.test)){
    par.test$permutation <- FALSE
  }
  if(!exists("pknown",par.model)){
    par.model$pknown <- FALSE
  }
  if(!exists("p",par.model)){
    par.model$p <- 0
  }
  if(!exists("max.p",par.model)){
    par.model$max.p <- 10
  }
  
  # read out parameters from par.test
  alpha <- par.test$alpha
  
  # number of variables
  d <- ncol(X)
  # number of observations
  n <- nrow(X)

  ###
  # Compute parent set
  ###

  # initialize variables
  parent.set <- 1:d
  len <- 1
  ind <- 1
  S <- list(numeric(0))
  # check wheter model is invariant over the empty set S={}
  if(!silent){
    message("Currently fitting set S = {}")
  }
  res <- seqICP.s(X, Y, S[[1]], test, par.test, model, par.model)
  result <- data.frame(ind=1,p.value=res$p.value,test.stat=res$test.stat,crit.value=res$crit.value, p=res$p)
  if(!silent){
    message(paste("p-value:",toString(round(res$p.value,digits=2))))
  }
  model.fits <- list(res$model.fit)
  # compute intersection of rejected sets
  if(res$p.value>alpha){
    parent.set <- numeric(0)
  }
  empty.int <- length(parent.set)==0
  # iterate over all sets until intersection is empty or max.parents is reached
  while(len<=max.parents & len<=d & !(empty.int&stopIfEmpty)){
    S <- append(S,combn(d,len,simplify=FALSE))
    while(ind<length(S) & !(empty.int&stopIfEmpty)){
      ind <- ind+1
      # use the candiate set S[[ind]]
      if(!silent){
        message(paste("Currently fitting set S = {",toString(S[[ind]]),"}",sep=""))
      }
      res <- seqICP.s(X, Y, S[[ind]], test, par.test, model, par.model)
      result <- rbind(result,data.frame(ind,p.value=res$p.value,test.stat=res$test.stat,crit.value=res$crit.value, p=res$p))
      if(!silent){
        message(paste("p-value:",toString(round(res$p.value,digits=2))))
      }
      model.fits[[ind]] <- res$model.fit
      # compute intersection of rejected sets
      if(res$p.value>alpha){
        parent.set <- intersect(parent.set,S[[ind]])
      }
      empty.int <- length(parent.set)==0
    }
    len <- len+1
  }
  # check if all tests where rejected (model is rejected)
  if(sum(result$p.value<=alpha)==length(result$p.value)){
    modelReject <- TRUE
    empty.int <- TRUE
    parent.set <- numeric(0)
  }
  else{
    modelReject <- FALSE
    if(empty.int){
      parent.set <- numeric(0)
    }
  }

  
  ###
  # Compute coefficients and confidence intervals for parent set
  ###
  
  # compute index of parent set
  parent.ind <- which(sapply(S,function(x) identical(x,parent.set)))
  parents <- S[[parent.ind]]
  accepted.ind <- result$ind[result$p.value>alpha]  
  # iid model
  if(model=="iid"){
    coef.mat <- matrix(0,d+1,3)
    # coefficients
    coef.tmp <- coef(model.fits[[parent.ind]])
    coef.mat[c(1,parents+1),1] <- coef.tmp
    # confidence intervals
    coef.mat[parents+1,2:3] <- NA
    coef.mat[1,2:3] <- confint(model.fits[[1]],1,1-alpha)
    for(ind in accepted.ind){
      fit.tmp <- model.fits[[ind]]
      int.int <- confint(fit.tmp,1,1-alpha)
      coef.mat[1,2:3] <- c(min(c(coef.mat[1,2],int.int[1])),max(c(coef.mat[1,3],int.int[2])))
      for(j in intersect(parents,S[[ind]])){
        k <- which(j==S[[ind]])
        interval <- confint(fit.tmp,k+1,1-alpha)
        coef.mat[j+1,2:3] <- c(min(c(coef.mat[j+1,2],interval[1]),na.rm=TRUE),max(c(coef.mat[j+1,3],interval[2]),na.rm=TRUE))
      }
    }
    colnames(coef.mat) <- c("coefficients","lower bound","upper bound")
  }
  else if(model=="ar"){
    p <- result$p[parent.ind]
    if(p>0){
      past.ind <- (d+2):((d+1)*(p+1))
    }
    else{
      past.ind <- numeric(0)
    }
    coef.mat <- matrix(0,(d+1)*(p+1),3)
    # coefficients
    coef.tmp <- coef(model.fits[[parent.ind]])
    coef.mat[c(1,parents+1,past.ind),1] <- coef.tmp
    # confidence intervals
    coef.mat[parents+1,2:3] <- NA
    coef.mat[1,2:3] <- confint(model.fits[[1]],1,1-alpha)
    coef.mat[past.ind,2:3] <- confint(model.fits[[parent.ind]],past.ind-(d-length(parents)),1-alpha)
    for(ind in accepted.ind){
      fit.tmp <- model.fits[[ind]]
      # intercept
      int.int <- confint(fit.tmp,1,1-alpha)
      coef.mat[1,2:3] <- c(min(c(coef.mat[1,2],int.int[1])),max(c(coef.mat[1,3],int.int[2])))
      # past
      int.past <- confint(fit.tmp,past.ind-(d-length(parents)),1-alpha)
      int.past[is.na(int.past)] <- 0
      coef.mat[c(rep(FALSE,nrow(coef.mat)-length(past.ind)),int.past[,1]<coef.mat[past.ind,2]),2] <- int.past[int.past[,1]<coef.mat[past.ind,2],1]
      coef.mat[c(rep(FALSE,nrow(coef.mat)-length(past.ind)),int.past[,2]>coef.mat[past.ind,3]),3] <- int.past[int.past[,2]>coef.mat[past.ind,3],2]
      for(j in intersect(parents,S[[ind]])){
        k <- which(j==S[[ind]])
        interval <- confint(fit.tmp,k+1,1-alpha)
        coef.mat[j+1,2:3] <- c(min(c(coef.mat[j+1,2],interval[1]),na.rm=TRUE),max(c(coef.mat[j+1,3],interval[2]),na.rm=TRUE))
      }
    }
    colnames(coef.mat) <- c("coefficients","lower bound","upper bound")
  }
  else{
    stop(paste("The model",model,"is unknown."))
  }
    
    
  
  ###
  # Compute p-values for all variables (only if stopIfEmpty==FALSE)
  ###

  if(!stopIfEmpty|length(S[[length(S)]])==d){
    p.value <- numeric(d)
    j.vec <- vector("logical",length(S))
    for(j in 1:d){
      for(i in 1:length(S)){
        j.vec[i] <- j %in% S[[i]]
      }
      pval1 <- result$p.value[j.vec]
      pval2 <- result$p.value[!j.vec]
      if(max(pval1)>alpha){
        p.value[j] <- max(pval2)
      }
      else{
        p.value[j] <- 1
      }
    }
  }
  else{
    p.value="not computed since stopIfEmpty==TRUE"
  } 
  

  icp.result <- list(parent.set=parent.set,
                      test.results=result,
                      S=S,
                      p.values=p.value,
                      coefficients=coef.mat,
                      stopIfEmpty=stopIfEmpty,
                      modelReject=modelReject,
                      pknown=par.model$pknown,
                      alpha=alpha,
                      n.var=d,
                      model=model)
  
  class(icp.result) <- "seqICP"
  
  return(icp.result)

}
