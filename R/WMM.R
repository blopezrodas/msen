###### Optimization Function
msen.WMM <- function(X, w = NULL, Rfun = "optimize") {
  if(is.vector(X))
    X <- matrix(X,ncol=1)
  if(is.data.frame(X))
    X <- as.matrix(X)
  if(any(is.na(X)))
    stop('No NAs allowed.')
  
  d <- ncol(X)
  n <- nrow(X)
  
  if(is.null(w))
    w <- rep(1,n)
  
  mu   <- apply(X, 2, stats::weighted.mean, w = w)
  Var  <- stats::cov.wt(x = X, wt = w, cor = FALSE, center = TRUE, method = "ML")$cov
  Kurt <- max(sum(w/sum(w)*stats::mahalanobis(x = X, center = mu, cov = Var, inverted = FALSE)^2), d*(d+2))
  
  # variance factor
  v <- function(theta) {
    theta*exp(theta)*expint::gammainc(a = 0, x = theta)
  }
  
  # kurtosis factor
  k <- function(theta) {
    (1 - theta*exp(theta)*expint::gammainc(a = 0, x = theta))/(theta*exp(2*theta)*expint::gammainc(a = 0, x = theta)^2)
  }
  
  # objective function
  ktheta <- Kurt/(d*(d+2))
  if(ktheta == 1) {
    warning("The best MSEN distribution is the normal one; for approximation, theta is fixed to 100.")
    theta <- 100
  }
  
  if(ktheta > 1) {
    if(Rfun == "optimize") {
      fobj <- function(theta, ktheta) {
        L2dist <- abs(k(theta = theta) - ktheta)^2
        return(L2dist)
      }
      res   <- stats::optimize(f = fobj, ktheta = ktheta, lower = 0, upper = 100) # , control=list(maxit = maxit, reltol = reltol, trace = 1)
      theta <- res$minimum
    }
    if(Rfun == "optim") {
      fobj <- function(theta, ktheta){
        theta.new <- exp(theta)
        L2dist <- abs(k(theta = theta.new) - ktheta)^2
        return(L2dist)
      }
      res   <- stats::optim(par = log(1), fn = fobj, ktheta = ktheta, method = "Nelder-Mead") # , control=list(maxit = maxit, reltol = reltol, trace = 1)
      theta <- exp(res$par)
    }
  }
  
  return(
    list(
      d     = d,
      mu    = mu,
      Var   = Var,
      Kurt  = Kurt,
      Sigma = Var/v(theta = theta),
      theta = theta
    )
  )
  
}