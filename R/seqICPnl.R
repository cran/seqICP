##' Estimates the causal parents S of the target variable Y using
##' invariant causal prediction and fits a general model of the form
##' \cr Y = f(X^S) + N.
##'
##' The function can be applied to models of the form\cr Y_i =
##' f(X_i^S) + N_i\cr with iid noise N_i and f is from a specific
##' function class, which the regression procedure given by the
##' parameter \code{regression.fun} should be able to approximate.\cr
##'
##' The invariant prediction procedure is applied using the hypothesis
##' test specified by the \code{test} parameter to determine whether a
##' candidate model is invariant. For further details see the
##' references.
##' @title Non-linear Invariant Causal Prediction
##' @param X matrix of predictor variables. Each column corresponds to
##'   one predictor variable.
##' @param Y vector of target variable, with length(Y)=nrow(X).
##' @param test string specifying the hypothesis test used to test for
##'   invariance of a parent set S (i.e. the null hypothesis
##'   H0_S). The following tests are available: "block.mean",
##'   "block.variance", "block.decoupled", "smooth.mean",
##'   "smooth.variance", "smooth.decoupled" and "hsic".
##' @param par.test parameters specifying hypothesis test. The
##'   following parameters are available: \code{grid}, \code{complements},
##'   \code{link}, \code{alpha} and \code{B}. The parameter \code{grid} is an increasing
##'   vector of gridpoints used to construct enviornments for change
##'   point based tests. If the parameter \code{complements} is 'TRUE' each
##'   environment is compared against its complement if it is 'FALSE'
##'   all environments are compared pairwise. The parameter \code{link}
##'   specifies how to compare the pairwise test statistics, generally
##'   this is either max or sum. The parameter \code{alpha} is a numeric
##'   value in (0,1) indicting the significance level of the
##'   hypothesis test. The parameter \code{B} is an integer and specifies
##'   the number of Monte-Carlo samples used in the approximation of
##'   the null distribution.
##' @param regression.fun regression function used to fit the function
##'   f. This should be a function which takes the argument (X,Y) and
##'   outputs the predicted values f(Y).
##' @param max.parents integer specifying the maximum size for
##'   admissible parents. Reducing this below the number of predictor
##'   variables saves computational time but means that the confidence
##'   intervals lose their coverage property.
##' @param stopIfEmpty if ‘TRUE’, the procedure will stop computing
##'   confidence intervals if the empty set has been accepted (and
##'   hence no variable can have a signicificant causal
##'   effect). Setting to‘TRUE’ will save computational time in these
##'   cases, but means that the confidence intervals lose their
##'   coverage properties for values different to 0.
##' @param silent If 'FALSE', the procedure will output progress
##'   notifications consisting of the currently computed set S
##'   together with the p-value resulting from the null hypothesis
##'   H0_S
##' 
##' @return object of class 'seqICPnl' consisting of the following
##'   elements
##'
##' \item{parent.set}{vector of the estimated causal parents.}
##' 
##' \item{test.results}{matrix containing results from each individual
##' hypothesis test H0_S as rows. The first column ind links the rows in the
##' \code{test.results} matrix to the position in the list of the variable
##' \code{S}.}
##'
##' \item{S}{list of all the sets that were
##' tested. The position within the list corresponds to the index in
##' the first column of the test.results matrix.}
##'
##' \item{p.values}{p-value for being not included in the set of true
##' causal parents. (If a p-value is smaller than alpha, the
##' corresponding variable is a member of parent.set.)}
##'
##' \item{stopIfEmpty}{a boolean value indicating whether computations
##' stop as soon as intersection of accepted sets is empty.}
##' 
##' \item{modelReject}{a boolean value indicating if the whole model was
##' rejected (the p-value of the best fitting model is too low).}
##'
##' \item{alpha}{significance level at which the hypothesis tests were performed.}
##'
##' \item{n.var}{number of predictor variables.}
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
##' @seealso The function \code{\link{seqICPnl.s}} can be used to
##'   perform the individual hypothesis tests
##'   H0_S. \code{\link{seqICP}} and \code{\link{seqICP.s}} are the
##'   corresponding functions for classical sequential invariant
##'   causal prediction.
##' 
##' @examples
##' set.seed(2)
##' 
##' # environment 1
##' na <- 120
##' X1a <- 0.3*rnorm(na)
##' X3a <- X1a + 0.2*rnorm(na)
##' Ya <- 2*X1a^2 + 0.6*sin(X3a) + 0.1*rnorm(na)
##' X2a <- -0.5*Ya + 0.5*X3a + 0.1*rnorm(na)
##' 
##' # environment 2
##' nb <- 80
##' X1b <- 2*rnorm(nb)
##' X3b <- rnorm(nb)
##' Yb <- 2*X1b^2 + 0.6*sin(X3b) + 0.1*rnorm(nb)
##' X2b <- -0.5*Yb + 0.8*rnorm(nb)
##' 
##' # combine environments
##' X1 <- c(X1a,X1b)
##' X2 <- c(X2a,X2b)
##' X3 <- c(X3a,X3b)
##' Y <- c(Ya,Yb)
##' Xmatrix <- cbind(X1, X2, X3)
##'
##' # use GAM as regression function
##' GAM <- function(X,Y){
##'   d <- ncol(X)
##'   if(d>1){
##'     formula <- "Y~1"
##'     names <- c("Y")
##'     for(i in 1:(d-1)){
##'       formula <- paste(formula,"+s(X",toString(i),")",sep="")
##'       names <- c(names,paste("X",toString(i),sep=""))
##'     }
##'     data <- data.frame(cbind(Y,X[,-1,drop=FALSE]))
##'       colnames(data) <- names
##'     fit <- fitted.values(mgcv::gam(as.formula(formula),data=data))
##'   } else{
##'     fit <- rep(mean(Y),nrow(X))
##'   }
##'   return(fit)
##' }
##' 
##' # Y follows the same structural assignment in both environments
##' # a and b (cf. the lines Ya <- ... and Yb <- ...).
##' # The direct causes of Y are X1 and X3.
##' # A GAM model fit considers X1, X2 and X3 as significant.
##' # All these variables are helpful for the prediction of Y.
##' summary(mgcv::gam(Y~s(X1)+s(X2)+s(X3)))
##'
##' # apply seqICP to the same setting
##' seqICPnl.result <- seqICPnl(X = Xmatrix, Y, test="block.variance",
##' par.test = list(grid = seq(0, na + nb, (na + nb)/10), complements = FALSE, link = sum,
##' alpha = 0.05, B =100), regression.fun = GAM,  max.parents = 4, stopIfEmpty=FALSE, silent=FALSE)
##' summary(seqICPnl.result)
##' # seqICPnl is able to infer that X1 and X3 are causes of Y

seqICPnl <- function(X, Y,
                 test="block.variance",
                 par.test=list(grid=c(0,round(nrow(X)/2),nrow(X)),
                               complements=FALSE,
                               link=sum,
                               alpha=0.05,
                               B=100),
                 regression.fun=function(X,Y) fitted.values(lm.fit(X,Y)),
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
  res <- seqICPnl.s(X, Y, S[[1]], test, par.test, regression.fun)
  result <- data.frame(ind=1,p.value=res$p.value,test.stat=res$test.stat,crit.value=res$crit.value, p=res$p)
  if(!silent){
    message(paste("p-value:",toString(round(res$p.value,digits=2))))
  }
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
      res <- seqICPnl.s(X, Y, S[[ind]], test, par.test, regression.fun)
      result <- rbind(result,data.frame(ind,p.value=res$p.value,test.stat=res$test.stat,crit.value=res$crit.value, p=res$p))
      if(!silent){
        message(paste("p-value:",toString(round(res$p.value,digits=2))))
      }
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
                     stopIfEmpty=stopIfEmpty,
                     modelReject=modelReject,
                     alpha=alpha,
                     n.var=d)
  
  class(icp.result) <- "seqICPnl"
  
  return(icp.result)

}
