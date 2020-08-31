data_cf<-function()
{
  # Generate data 
  # parameters for data generating
  p <- 10 # number of total covariates
  pt <- 4 # number of covariates affecting treatment effects
  py <- 4 # number of covariates affecting outcomes but not treatment effects
  asym <- .5 # whether treatment effects are distributed asymmetrically across treated and control
  n <- 3000 # total size of the dataset
  propens <- .5 #treatment probability
  sig = .01
  treatsize <- .5 # treatment effect size
  levsize <- 1
  # draw W
  w <- rbinom(n, 1, propens)
  w_ <- 1-w
  # draw X
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  # generate treatment effects as function of X
  if (p<pt+py) print("error: p>=pt+py required")
  tau <- 0
  for (iii in 1:pt) {
    tau <- tau + treatsize*pmax(X[,iii],array(0,n))*(2*X[,iii])
  }
  # generate average value of outcomes
  mu <- treatsize*rowSums(X[,1:pt])+levsize*rowSums(X[,(pt+1):(pt+py)])
  
  # generate outcomes as function of treatment status, mu, tau, and noise
  y <- mu + asym*w*tau + (asym-1)*(1-w)*tau + rnorm(n,0,sig)
  y_ <- mu + asym*w_*tau + (asym-1)*(1-w_)*tau + rnorm(n,0,sig)
  
  # create formulas for estimation
  f <- ""
  nextx <- ""
  if (p>1) {
    for (ii in 1:(p-1)) {
      nextx <- paste("x",ii, sep="")
      if (ii==1) {name <- nextx}
      if (ii>1) {name <- c(name, nextx)}
      f <- paste(f, nextx, "+", sep="")
    }
    f <- paste(f, "x", ii+1, sep="")
  } else if (p==1) {
    f <- "x1"
  }
  
  for (ii in 1:p) {
    nextx <- paste("x",ii, sep="")
    if (ii==1) {name <- nextx}
    if (ii>1) {name <- c(name, nextx)}
  }
  nameall <- c( name,  "y", "w", "tau_true")
  tau_true <- (1-2*w)*(y_ - y)
  ntr <- round(.9*n)
  ntest <- n - ntr 
  dataTrain <- data.frame(X[1:ntr,], y[1:ntr], w[1:ntr], tau_true[1:ntr])
  dataTest <- data.frame(X[(ntr+1):n,], y[(ntr+1):n], w[(ntr+1):n], tau_true[(ntr+1):n])
  names(dataTrain)=nameall
  names(dataTest)=nameall
  data_cf<-as.list(NULL)
  data_cf[[1]]<-dataTrain
  data_cf[[2]]<-dataTest
  data_cf[[3]]<-p
  data_cf[[4]]<-f
  return (data_cf)
}