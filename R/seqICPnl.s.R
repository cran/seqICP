##' Tests whether the conditional distribution of Y given X^S is
##' invariant across time, by allowing for arbitrary non-linear
##' additive dependence models.
##'
##' The function can be applied to models of the form\cr Y_i =
##' f(X_i^S) + N_i\cr with iid noise N_i and f is from a specific
##' function class, which the regression procedure given by the
##' parameter \code{regression.fun} should be able to approximate.\cr
##'
##' For both models the hypothesis test specified by the \code{test}
##' parameter specifies the hypothesis test used to test whether the
##' set S leads to an invariant model. For futher details see the
##' references.
##' @title Non-linear Invariant Causal Prediction for an individual
##'   set S
##' @param X matrix of predictor variables. Each column corresponds to
##'   one predictor variable.
##' @param Y vector of target variable, with length(Y)=nrow(X).
##' @param S vector containing the indicies of predictors to be tested
##' @param test string specifying the hypothesis test used to test for
##'   invariance of a parent set S (i.e. the null hypothesis
##'   H0_S). The following tests are available: "block.mean",
##'   "block.variance", "block.decoupled", "smooth.mean",
##'   "smooth.variance", "smooth.decoupled" and "hsic".
##' @param par.test parameters specifying hypothesis test. The
##'   following parameters are available: \code{grid},
##'   \code{complements}, \code{link}, \code{alpha} and \code{B}. The
##'   parameter \code{grid} is an increasing vector of gridpoints used
##'   to construct enviornments for change point based tests. If the
##'   parameter \code{complements} is 'TRUE' each environment is
##'   compared against its complement if it is 'FALSE' all
##'   environments are compared pairwise. The parameter \code{link}
##'   specifies how to compare the pairwise test statistics, generally
##'   this is either max or sum. The parameter \code{alpha} is a
##'   numeric value in (0,1) indicting the significance level of the
##'   hypothesis test. The parameter \code{B} is an integer and
##'   specifies the number of Monte-Carlo samples used in the
##'   approximation of the null distribution.
##' @param regression.fun regression function used to fit the function
##'   f. This should be a function which takes the argument (X,Y) and
##'   outputs the predicted values f(Y).
##' 
##' @return list containing the following elements
##'
##' \item{test.stat}{value of the test statistic.}
##'
##' \item{crit.value}{critical value computed using a Monte-Carlo
##' simulation of the null distribution.}
##'
##' \item{p.value}{p-value.}
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
##' @seealso To estimate the set of causal parents use the function
##'   \code{\link{seqICPnl}}. For linear models use the corresponding
##'   functions \code{\link{seqICP}} and \code{\link{seqICP.s}}.
##' 
##' @examples
##' set.seed(1)
##'
##' # environment 1
##' na <- 130
##' X1a <- rnorm(na,0,0.1)
##' Ya <- 5*X1a+rnorm(na,0,0.5)
##' X2a <- Ya+rnorm(na,0,0.1)
##'
##' # environment 2
##' nb <- 70
##' X1b <- rnorm(nb,-1,1)
##' Yb <- X1b^2+rnorm(nb,0,0.5)
##' X2b <- rnorm(nb,0,0.1)
##' 
##' # combine environments
##' X1 <- c(X1a,X1b)
##' X2 <- c(X2a,X2b)
##' Y <- c(Ya,Yb)
##' Xmatrix <- cbind(X1, X2)
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
##' # apply seqICPnl.s to all possible sets using the regression
##' # function GAM - only the true parent set S=1 is
##' # invariant in this example
##' seqICPnl.s(Xmatrix, Y, S=numeric(), par.test=list(grid=c(0,50,100,150,200)), regression.fun=GAM)
##' seqICPnl.s(Xmatrix, Y, S=1, par.test=list(grid=c(0,50,100,150,200)), regression.fun=GAM)
##' seqICPnl.s(Xmatrix, Y, S=2, par.test=list(grid=c(0,50,100,150,200)), regression.fun=GAM)
##' seqICPnl.s(Xmatrix, Y, S=c(1,2), par.test=list(grid=c(0,50,100,150,200)), regression.fun=GAM)

seqICPnl.s <- function(X, Y, S,
                   test="block.variance",
                   par.test=list(grid=c(0,round(nrow(X)/2),nrow(X)),
                                 complements=FALSE,
                                 link=sum,
                                 alpha=0.05,
                                 B=100),
                   regression.fun=function(X,Y) fitted.values(lm.fit(X,Y)))
{

  # add defaults if missing in parameter lists
  if(!exists("grid",par.test)){
    par.test$grid <- seq(0,nrow(X),length.out=3)
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
  grid <- par.test$grid
  complements <- par.test$complements
  link <- par.test$link
  alpha <- par.test$alpha
  B <- par.test$B
  
  # number of observations
  n <- nrow(X)
  # convert vectors to matrices
  Y <- as.matrix(Y)
  
  ###
  # Adjust X, Y based on model
  ###
  
  # compute new design matrix
  X <- cbind(matrix(1,n,1),X[,S,drop=FALSE])
  # number of variables
  d <- ncol(X)

  
  ###
  # Compute block environments if necessary
  ###
  if(test=="block.mean"||test=="block.variance"||test=="block.decoupled"){
    # Comparisons against the complement or general comparisons
    if(complements){
      # Compare environments with complement environments
      # number of grid points
      L <- length(grid)
      # check whether grid was specified
      if(L<=1){
        stop("not enough grid points given, specify correct grid")
      }
      # if first grid point is 1 change it to 0
      if(grid[1]==1){
        grid[1] <- 0
      }
      # total number of environments
      L.tot <- choose(L,2)
      # compute all combinations of elements on the grid
      t <- combn(L,2,simplify=TRUE)
      # initialize variables
      Z <- vector("list",2*L.tot)
      env <- vector("list",2*L.tot)
      fullenv <- L-1
      for(i in 1:L.tot){
        cp1 <- grid[t[1,i]]
        cp2 <- grid[t[2,i]]
        env[[i]] <- (cp1+1):cp2
        if(cp1==0 & cp2!=n){
          env[[i+L.tot]] <- (cp2+1):n
        }
        else if(cp2==n & cp1!=0){
          env[[i+L.tot]] <- 1:cp1
        }
        else if(cp1==0 & cp2==n){
          env[[i+L.tot]] <- 1:n
        }
        else{
          env[[i+L.tot]] <- c(1:cp1,(cp2+1):n)
        }
      }
      # write grid to CP
      CP <- grid
      # compute relevant comparisons
      pairwise <- cbind((1:L.tot)[-fullenv],(1:L.tot)[-fullenv]+L.tot)
      # set L.tot to the full number of environments
      L.tot <- 2*L.tot
    }
    else{
      # Compare environements against all non-intersecting environments
      # number of grid points
      L <- length(grid)
      # check whether grid was specified
      if(L<=1){
        stop("not enough grid points given, specify correct grid")
      }
      # if first grid point is 1 change it to 0
      if(grid[1]==1){
        grid[1] <- 0
      }
      # total number of environments
      L.tot <- choose(L,2)
      # compute all combinations of elements on the grid
      t <- combn(L,2,simplify=TRUE)
      # initialize variables
      Z <- vector("list",L.tot)
      env <- vector("list",L.tot)
      for(i in 1:L.tot){
        cp1 <- grid[t[1,i]]
        cp2 <- grid[t[2,i]]
        env[[i]] <- (cp1+1):cp2
      }
      # write grid to CP
      CP <- grid
      # compute relevant comparisons
      index <- t[1,]>=t[2,1]
      pairwise <- rbind(rep(1,sum(index)),(1:L.tot)[index])
      for(i in 2:L.tot){
        index <- t[1,]>=t[2,i]
        pairwise <- cbind(pairwise,rbind(rep(i,sum(index)),(1:L.tot)[index]))
      }
      pairwise <- t(cbind(pairwise,rbind(pairwise[2,],pairwise[1,])))
    }
  }
  
  ###
  # Define function to compute test statistic (depends on test)
  ###
  if(test=="hsic"){
    test.stat.fun <- function(R.scaled){
      test.stat <- dHSIC::dhsic(1:n,R.scaled)$dHSIC
      return(test.stat)
    }
    # size of output
    num <- 1
  }
  else if(test=="block.mean"){
    test.stat.fun <- function(R.scaled){
      mean.env <- vector("numeric",L.tot)
      for(i in 1:L.tot){
        mean.env[i] <- mean(R.scaled[env[[i]],,drop=FALSE])
      }
      # go over all pairs of non-intersecting environments
      statfun <- function(i,j){
        T <- abs(mean.env[i]-mean.env[j])
        return(T)
      }
      T.pairs <- mapply(statfun,pairwise[,1],pairwise[,2])
      # combine the pairwise test statistics with the link function
      test.stat <- link(T.pairs)
      return(test.stat)
    }
    # size of output
    num <- 1
  }
  else if(test=="block.variance"){
    test.stat.fun <- function(R.scaled){
      var.env <- vector("numeric",L.tot)
      for(i in 1:L.tot){
        var.env[i] <- var(R.scaled[env[[i]],,drop=FALSE])
      }
      # go over all pairs of non-intersecting environments
      statfun <- function(i,j){
        T <- abs(var.env[i]/var.env[j]-1)
        return(T)
      }
      T.pairs <- mapply(statfun,pairwise[,1],pairwise[,2])
      # combine the pairwise test statistics with the link function
      test.stat <- link(T.pairs)
      return(test.stat)
    }
    # size of output
    num <- 1
  }
  else if(test=="block.decoupled"){
    test.stat.fun <- function(R.scaled){
      var.env <- vector("numeric",L.tot)
      mean.env <- vector("numeric",L.tot)
      for(i in 1:L.tot){
        var.env[i] <- var(R.scaled[env[[i]],,drop=FALSE])
        mean.env[i] <- mean(R.scaled[env[[i]],,drop=FALSE])
      }
      # go over all pairs of non-intersecting environments
      statfun <- function(i,j){
        T1 <- abs(mean.env[i]-mean.env[j])
        T2 <- abs(var.env[i]/var.env[j]-1)
        T <- c(T1,T2)
        return(T)
      }
      T.pairs <- mapply(statfun,pairwise[,1],pairwise[,2])
      # combine the pairwise test statistics with the link function
      test.stat <- rbind(link(T.pairs[1,]),link(T.pairs[2,]))
      return(test.stat)
    }
    # size of output
    num <- 2
  }
  else if(test=="smooth.variance"){
    test.stat.fun <- function(R.scaled){
      time <- 1:n
      fit <- mgcv::gam(R.scaled^2~s(time))
      y <- fitted.values(fit)
      test.stat <- sum(y^2)
      return(test.stat)
    }
    # size of output
    num <- 1
  }
  else if(test=="smooth.mean"){
    test.stat.fun <- function(R.scaled){
      time <- 1:n
      fit <- mgcv::gam(R.scaled~s(time))
      y <- fitted.values(fit)
      test.stat <- sum(y^2)
      return(test.stat)
    }
    # size of output
    num <- 1
  }
  else if(test=="smooth.decoupled"){
    test.stat.fun <- function(R.scaled){
      time <- 1:n
      fit.m <- mgcv::gam(R.scaled~s(time))
      fit.v <- mgcv::gam(R.scaled^2~s(time))
      y.m <- fitted.values(fit.m)
      y.v <- fitted.values(fit.v)
      test.stat <- c(sum(y.m^2),sum(y.v^2))
      return(test.stat)
    }
    # size of output
    num <- 2
  }
  else{
    stop("given test unknown")
  }

  
  
  ###  
  # Compute test statistic
  ###
  Ypred <- regression.fun(X,Y)
  R <- Y-Ypred
  R.scaled <- R/sqrt(sum(R^2))
  test.stat <- test.stat.fun(R.scaled)

  ## ###
  ## # Diagnostic plots
  ## ###
  ## if(verbose){
  ##   ### plot the cross correlation in the after ar fit
  ##   fit <- lm(Y~-1+X)
  ##   res <- residuals(fit)
  ##   par(mfrow=c(1,2))
  ##   ccf(Y[,1],res)
  ##   acf(res)
  ##   par(mfrow=c(1,1))
  ##   #readline(prompt="Press [enter] to continue")
    
  ##   ### plot the scaled residuals versus time
  ##   time <- 1:n
  ##   Rs <- list(R.scaled,p)
  ##   #save(Rs,file="Rs.Rda")
  ##   plot(R.scaled^2)
  ##   #lines(fitted.values(smooth.spline(time,R.scaled^2)),col="red")
  ##   #lines(lokerns(time,R.scaled^2,hetero=TRUE)$y,col="red")
  ##   lines(fitted.values(gam(R.scaled^2~s(time))),col="red")
  ## }
  
  ###
  # Use permutation of the scaled residuals to approximate null distribution
  ###
  boot.fun <- function(i){
    return(test.stat.fun(matrix(R.scaled[sample(1:n)],n,1)))
  }
  test.stat.boot <- matrix(sapply(1:B,boot.fun),num,B)
  crit.value <- apply(test.stat.boot,1,function(x) quantile(x,1-alpha))
  p.value <- (rowSums(test.stat.boot>=matrix(rep(test.stat,B),num,B))+1)/(B+1)
  
  # perform Bonferroni adjustment incase method uses multiple testing
  min.ind <- which.min(p.value)
  test.stat <- test.stat[min.ind]
  crit.value <- crit.value[min.ind]
  p.value <- min(num*p.value[min.ind],1)
  
  return(list(test.stat=as.numeric(test.stat),
              crit.value=crit.value,
              p.value=as.numeric(p.value)))
         
}
  
