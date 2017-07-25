##' Tests whether the conditional distribution of Y given X^S is
##' invariant across time, by assuming a linear dependence model.
##'
##' The function can be applied to two types of models\cr
##' (1) a linear model (model="iid")\cr Y_i = a X_i^S + N_i\cr with iid noise N_i and \cr
##' (2) a linear autoregressive model (model="ar")\cr Y_t = a_0 X_t^S + ... + a_p (Y_(t-p),X_(t-p)) + N_t\cr with iid noise N_t.
##'
##' For both models the hypothesis test specified by the \code{test}
##' parameter is used to test whether the set S leads to an invariant
##' model. For futher details see the references.
##' @title Sequential Invariant Causal Prediction for an individual
##'   set S
##' @param X matrix of predictor variables. Each column corresponds to
##'   one predictor variable.
##' @param Y vector of target variable, with length(Y)=nrow(X).
##' @param S vector containing the indicies of predictors to be tested
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
##' \item{p}{number of lags that were used.}
##'
##' \item{model.fit}{'lm' object of linear model fit.}
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
##'   \code{\link{seqICP}}. For non-linear models use the
##'   corresponding functions \code{\link{seqICPnl}} and
##'   \code{\link{seqICPnl.s}}.
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
##' Yb <- 5*X1b+rnorm(nb,0,0.5)
##' X2b <- rnorm(nb,0,0.1)
##'
##' # combine environments
##' X1 <- c(X1a,X1b)
##' X2 <- c(X2a,X2b)
##' Y <- c(Ya,Yb)
##' Xmatrix <- cbind(X1, X2)
##'
##' # apply seqICP.s to all possible sets - only the true parent set S=1
##' # is invariant in this example
##' seqICP.s(Xmatrix, Y, S=numeric(), par.test=list(grid=c(0,50,100,150,200)))
##' seqICP.s(Xmatrix, Y, S=1, par.test=list(grid=c(0,50,100,150,200)))
##' seqICP.s(Xmatrix, Y, S=2, par.test=list(grid=c(0,50,100,150,200)))
##' seqICP.s(Xmatrix, Y, S=c(1,2), par.test=list(grid=c(0,50,100,150,200)))

seqICP.s <- function(X, Y, S,
                   test="decoupled",
                   par.test=list(grid=c(0,round(nrow(X)/2),nrow(X)),
                                 complements=FALSE,
                                 link=sum,
                                 alpha=0.05,
                                 B=100,
                                 permutation=FALSE),
                   model="iid",
                   par.model=list(pknown=FALSE, p=0, max.p=10))
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
  grid <- par.test$grid
  complements <- par.test$complements
  link <- par.test$link
  alpha <- par.test$alpha
  B <- par.test$B
  permutation <- par.test$permutation
  
  # number of observations
  n <- nrow(X)
  # convert vectors to matrices
  Y <- as.matrix(Y)
  
  ###
  # Adjust X, Y based on model
  ###
  
  # iid model
  if(model=="iid"){
    # compute new design matrix
    X <- cbind(matrix(1,n,1),X[,S,drop=FALSE])
    # number of variables
    d <- ncol(X)
    p <- NA
  }
  # ar model
  else if(model=="ar"){
    d <- ncol(X)
    # check wether p is known
    if(par.model$pknown){
      p <- par.model$p
    }
    else{
      # use AIC/BIC to choose p (check for all p upto max.p)
      max.p <- par.model$max.p
      #AICc <- numeric(max.p+1)
      AIC <- numeric(max.p+1)
      #BIC <- numeric(max.p+1)
      Xtmp <- cbind(matrix(1,n,1),X[,S,drop=FALSE])
      # fit model and compute residuals
      res <- lm.fit(x=Xtmp,y=Y)$residuals
      ntmp <- n
      res.var <- sum(res^2)/n
      # AIC/BIC
      ktmp <- length(S)+1
      #AICc[1] <- n*(log(2*pi*res.var)+1)+2*(ntmp*ktmp)/(ntmp-ktmp-1)
      AIC[1] <- 2*length(S)-2*sum(log(1/sqrt(2*pi*res.var))-res^2/(2*res.var))
      #BIC[1] <- log(n)*length(S)-2*sum(log(1/sqrt(2*pi*res.var))-res^2/(2*res.var))
      for(k in 1:max.p){
        ntmp <- n-k
        Xnew <- matrix(0,ntmp,(d+1)*k+length(S))
        if(length(S)>0){
          for(i in 1:length(S)){
            Xnew[,i] <- X[(k+1):n,S[i]]
          }
        }
        for(i in 1:k){
          Xnew[,(length(S)+(d+1)*(i-1)+1)] <- Y[(k+1-i):(n-i)]
          Xnew[,(length(S)+(d+1)*(i-1)+2):(length(S)+(d+1)*i)] <- X[(k+1-i):(n-i),]
        }
        Xtmp <- cbind(matrix(1,n-k,1),Xnew)
        Ytmp <- Y[(k+1):n,,drop=FALSE]
        # fit model and compute residuals
        res <- lm.fit(x=Xtmp,y=Ytmp)$residuals
        #df <- n-ncol(Xtmp)-1
        res.var <- sum(res^2)/ntmp
        # AIC/BIC
        ktmp <- (d+1)*k+length(S)+1
        #AICc[k+1] <- n*(log(2*pi*res.var)+1)+2*(ntmp*ktmp)/(ntmp-ktmp-1)
        AIC[k+1] <- 2*(length(S)+k*(d+1))-2*sum(log(1/sqrt(2*pi*res.var))-res^2/(2*res.var))
        #BIC[k+1] <- log(n-k)*(length(S)+k*(d+1))-2*sum(log(1/sqrt(2*pi*res.var))-res^2/(2*res.var))
      }
      #p <- which.min(AICc)-1
      p <- which.min(AIC)-1
      #p <- which.min(BIC)-1
      #message(paste("minimal AIC for p=",toString(p),sep=""))
    }
    # compute new design matrix X
    if(p>0){
      Xnew <- matrix(0,n-p,(d+1)*p+length(S))
      if(length(S)>0){
        for(i in 1:length(S)){
          Xnew[,i] <- X[(p+1):n,S[i]]
        }
      }
      for(i in 1:p){
        Xnew[,(length(S)+(d+1)*(i-1)+1)] <- Y[(p+1-i):(n-i)]
        Xnew[,(length(S)+(d+1)*(i-1)+2):(length(S)+(d+1)*i)] <- X[(p+1-i):(n-i),]
      }
      X <- cbind(rep(1,n-p),Xnew)
    }
    else{
      X <- cbind(rep(1,n),X[,S,drop=FALSE])
    }
    # compute new response Y
    Y <- Y[(p+1):n,,drop=FALSE]
    # adjust number of variables
    d <- ncol(X)
    # adjust sample size
    n <- n-p
  }
  else{
    stop("given model unknown")
  }


  ###
  # Compute regression matricies
  ###
  if(test=="trend"||test=="combined"||test=="decoupled"||test==
       "variance"||test=="block.decoupled"||test=="block.mean"||test==
       "block.variance"){
    # Comparisons against the complement or general comparisons
    if(complements){
      # Compare environments with complement environments
      # check whether grid is adjusted to p if model=="ar"
      if(model=="ar"& grid[length(grid)]>n){
        warning(paste("grid points larger than ",n ,
                      " are removed, to account for lag p=",p,sep=""))
        grid <- grid[grid<=n]
        if(grid[length(grid)]!=n){
          grid <- c(grid,n)
        }
      }
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
        X.comp <- X[env[[i+L.tot]],,drop=FALSE]
        X.env <- X[env[[i]],,drop=FALSE]
        Z[[i]] <- solve(t(X.env)%*%X.env,t(X.env))
        Z[[i+L.tot]] <- solve(t(X.comp)%*%X.comp,t(X.comp))
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
      # check whether grid is adjusted to p if model=="ar"
      if(model=="ar"& grid[length(grid)]>n){
        warning(paste("grid points larger than ",n ,
                      " are removed, to account for lag p=",p,sep=""))
        grid <- grid[grid<=n]
        if(grid[length(grid)]!=n){
          grid <- c(grid,n)
        }
      }
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
        X.env <- X[env[[i]],,drop=FALSE]
        Z[[i]] <- solve(t(X.env)%*%X.env,t(X.env))
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
  if(test=="combined"){
    test.stat.fun <- function(R.scaled){
      R.env <- vector("list",L.tot)
      X.env <- vector("list",L.tot)
      var.env <- vector("numeric",L.tot)
      gamma.env <- vector("list",L.tot)
      size <- vector("numeric",L.tot)
      for(i in 1:L.tot){
        R.env[[i]] <- R.scaled[env[[i]],,drop=FALSE]
        X.env[[i]] <- X[env[[i]],,drop=FALSE]
        gamma.env[[i]] <- Z[[i]]%*%R.env[[i]]
        resid <- R.env[[i]]-X.env[[i]]%*%gamma.env[[i]]
        size[i] <- length(env[[i]])
        var.env[i] <- sum(resid^2)/(size[i]-d)
      }
      # go over all pairs of non-intersecting environments
      statfun <- function(i,j){
        T <- abs(sum((R.env[[i]]-X.env[[i]]%*%gamma.env[[j]])^2)/(var.env[j]*size[i])-1)
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
  else if(test=="decoupled"){
    test.stat.fun <- function(R.scaled){
      var.env <- vector("numeric",L.tot)
      gamma.env <- vector("list",L.tot)
      for(i in 1:L.tot){
        R.env <- R.scaled[env[[i]],,drop=FALSE]
        X.env <- X[env[[i]],,drop=FALSE]
        gamma.env[[i]] <- Z[[i]]%*%R.env
        resid <- R.env-X.env%*%gamma.env[[i]]
        var.env[i] <- sum(resid^2)/(length(env[[i]])-d)
      }
      # go over all pairs of non-intersecting environments
      statfun <- function(i,j){
        T1 <- sum((gamma.env[[i]]-gamma.env[[j]])^2)
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
  else if(test=="variance"){
    test.stat.fun <- function(R.scaled){
      var.env <- vector("numeric",L.tot)
      gamma.env <- vector("list",L.tot)
      for(i in 1:L.tot){
        R.env <- R.scaled[env[[i]],,drop=FALSE]
        X.env <- X[env[[i]],,drop=FALSE]
        gamma.env[[i]] <- Z[[i]]%*%R.env
        resid <- R.env-X.env%*%gamma.env[[i]]
        var.env[i] <- sum(resid^2)/(length(env[[i]])-d)
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
  else if(test=="trend"){
    test.stat.fun <- function(R.scaled){
      gamma.env <- vector("list",L.tot)
      for(i in 1:L.tot){
        R.env <- R.scaled[env[[i]],,drop=FALSE]
        gamma.env[[i]] <- Z[[i]]%*%R.env
      }
      # go over all pairs of non-intersecting environments
      statfun <- function(i,j){
        T <- sum((gamma.env[[i]]-gamma.env[[j]])^2)
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
  else if(test=="hsic"){
    test.stat.fun <- function(R.scaled){
      test.stat <- dHSIC::dhsic(1:n,R.scaled)$dHSIC
      return(test.stat)
    }
    # size of output
    num <- 1
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
        var.env[i] <- mean(R.scaled[env[[i]],,drop=FALSE])
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
  else{
    stop("given test unknown")
  }

  
  ###
  # Fit linear model
  ###
  model.fit <- lm(Y~-1+X)
  
  ###  
  # Compute test statistic
  ###
  M <- diag(n)-X%*%solve(t(X)%*%X,t(X))
  R <- M%*%Y
  R.scaled <- R/sqrt(sum(R^2))
  test.stat <- test.stat.fun(R.scaled)

  ###
  # Diagnostic plots
  ###
  
  ### plot the cross correlation in the after ar fit
  ## fit <- lm(Y~-1+X)
  ## res <- residuals(fit)
  ## par(mfrow=c(1,2))
  ## ccf(Y[,1],res)
  ## acf(res)
  ## par(mfrow=c(1,1))
  #readline(prompt="Press [enter] to continue")
  
  ### plot the scaled residuals versus time
  ## time <- 1:n
  ## Rs <- list(R.scaled,p)
  ## #save(Rs,file="Rs.Rda")
  ## graphics::plot(R.scaled^2)
  ## #graphics::lines(fitted.values(smooth.spline(time,R.scaled^2)),col="red")
  ## #graphics::lines(lokern::glkerns(time,R.scaled^2,x.out=time)$est,col="red")
  ## graphics::lines(fitted.values(mgcv::gam(R.scaled^2~s(time))),col="red")
  ## #readline(prompt="Press [enter] to continue")
  
  
  ###
  # Use parametric bootstrap of the scaled residuals to approximate the null distribution
  ###
  if(permutation){
    boot.fun <- function(i){
      R.scaled.perm <- matrix(R.scaled[sample(1:n)],n,1)
      return(test.stat.fun(R.scaled.perm))
    }
  }
  else{
    boot.fun <- function(i){
      tmp <- M%*%matrix(rnorm(n),n,1)
      R.scaled.boot <- tmp/sqrt(sum(tmp^2))
      return(test.stat.fun(R.scaled.boot))
    }
  }
  test.stat.boot <- matrix(sapply(1:B,boot.fun),num,B)
  crit.value <- apply(test.stat.boot,1,function(x) quantile(x,1-alpha))
  p.value <- (rowSums(test.stat.boot>=matrix(rep(test.stat,B),num,B))+1)/(B+1)

  # perform Bonferroni adjustment if required
  min.ind <- which.min(p.value)
  test.stat <- test.stat[min.ind]
  crit.value <- crit.value[min.ind]
  p.value <- min(num*p.value[min.ind],1)

  return(list(test.stat=as.numeric(test.stat),
              crit.value=crit.value,
              p.value=as.numeric(p.value),
              p=p,
              model.fit=model.fit))
  
}
