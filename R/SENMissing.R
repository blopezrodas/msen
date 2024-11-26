library(stats) 
library(mvtnorm) 
library(mnormt) 
library(expint)
library(MASS)
library(matrixcalc)


SEN_m<-function(X){
formula = "indirect"
reltol = 10^-4
max_iter= 200
epsilon=0.0001
check     <- 0
iter <- 0
loglik    <- NULL

n <- nrow(X) # number of observations
d <- ncol(X) # number of variables/dimensions

 # number of observed values per observation
  mv<- m <- is.na(X)
  ov <- !mv
  do<-apply(ov,1,sum)

# definition for the E-step

w     <- numeric(n)
delta <- numeric(n)

# Distribution

dens <- numeric(n)


X_tilde     <- X
sigma_tilde <- array(NA, dim = c(d, d, n))

Xcom=X[do==d,]
Xinc=X[!do==d,]

temp  <- msen.WMM(X = Xcom)
mu    <- temp$mu
Sigma <- temp$Sigma
theta <- temp$theta
if(!is.positive.definite(Sigma))Sigma=diag(d)
#### E step

while(iter < max_iter & getall(loglik) > epsilon){
  iter=iter+1
for (i in 1:n) {
  xi <- X[i, ]
  m <- is.na(xi)
  o <- !m
  
  mu_o <- mu[ o]
  sigma_oo <- as.matrix(Sigma[o, o])
  # sigma_oo <- round(sigma_oo, digits = 16)
  delta <- as.vector(xi[o]-mu_o) %*% solve(sigma_oo) %*% as.vector(xi[o]-mu_o)

  num <- expint::gammainc(a = (do[i]/2 + 2), x = (delta/2 + theta))
  den <- (delta/2 + theta)*expint::gammainc(a = (do[i]/2 + 1), x = (delta/2 + theta))
 
  
   w[i] <- num/den
  
  
  if (any(m)) {
  
    
    mu_m <- mu[ m]
    mu_o <- mu[ o]
    
    sigma_mo <- Sigma[m, o]
    sigma_om <- Sigma[o, m]
    sigma_mm <- Sigma[m, m]
    sigma_oo_inv <- mnormt::pd.solve(Sigma[o, o])
    # sigma_oo_inv <- solve(sigma[o, o, g])
    
    x_ig_tilde <- mu_m + sigma_mo %*% sigma_oo_inv %*% (xi[o] - mu_o)
    X_tilde[i, m] <- x_ig_tilde
    
    sigma_tilde[o, o, i] <- tcrossprod((xi[o] - mu_o))
    sigma_tilde[o, m, i] <- tcrossprod((xi[o] - mu_o), (x_ig_tilde - mu_m))
    sigma_tilde[m, o, i] <- t(sigma_tilde[o, m, i])
    M <- sigma_mm - sigma_mo %*% sigma_oo_inv %*% sigma_om
    sigma_tilde[m, m, i] <- tcrossprod((x_ig_tilde - mu_m)) + M / w[i]
    sigma_tilde[,,i]=w[i] * sigma_tilde[,,i]
  } else {
    sigma_tilde[, , i] <- tcrossprod(sqrt(w[i])*unlist(xi - mu))
  }
}

#######M step

mu    <- colSums(w/sum(w)*X_tilde)
  
#Sigma  <- apply(w * sigma_tilde, 1:2, sum)/ n
Sigma  <- apply(sigma_tilde, 1:2, sum)/ n



# if (!isSymmetric.matrix(Sigma) | max(abs(Sigma - t(Sigma))) > .Machine$double.eps) {
#   matr <- Sigma
#   matr[lower.tri(matr)] <- t(matr)[lower.tri(t(matr))]
#   Sigma <- matr
# }

theta <- n/sum(w-1)
dens <- dmsen(x = Xcom, mu = mu, Sigma = Sigma, theta = theta, formula = formula)
densI=NULL
md=0
if(any(mv)){
for(i in 1:nrow(Xinc)){
  xi <- Xinc[i, ]
  o <- !is.na(Xinc[i, ])
 
  
  mu_o <- mu[ o]
  sigma_oo <- as.matrix(Sigma[o, o])
  
  densI[i] <- dmsen(t(xi[o]), mu = mu_o, Sigma = sigma_oo, theta = theta, formula = formula)
  
}
  md=sum(log(densI))
         }




# ------------------------------------- #
# Global - Observed-data log-likelihood #
# ------------------------------------- #

llvalues <- sum(log(dens))+md

loglik[iter] <- llvalues


}

plot(loglik)
print(all(loglik==sort(loglik)))
return(list(mu=mu,Sigma=Sigma, theta=theta, X_tilde=X_tilde,loglik=loglik))

}

######USED FUNCTIONS

dmsen <- function(x, mu = rep(0, d), Sigma, theta = Inf, formula = "direct") {
  if (missing(Sigma)) {
    stop("Sigma is missing")
  }
  if (theta < 0) {
    stop("theta must be greater than, or equal to, 0")
  }
  if (is.matrix(Sigma)) {
    d <- ncol(Sigma)
  }
  if (!is.matrix(Sigma)) {
    d <- 1
  }
  if (is.vector(x)) {
    x <- matrix(x, length(x), 1)
    Sigma <- matrix(Sigma, nrow = d, ncol = d)
  }
  if (formula == "direct") {
    delta <- sapply(1:nrow(x), function(i) t(as.vector(t(x[i, ]) - mu)) %*% solve(Sigma) %*% as.vector(t(x[i, ]) - mu))
    delta <- replace(delta, delta == 0, 1 / (theta * (2 *
                                                        pi)^(d / 2) * (d / 2 + 1)) * (1 - (1 - theta)^(d / 2 +
                                                                                                         1)))
    pdfgamma <- expint::gammainc(a = (d / 2 + 1), x = 1 / 2 *
                                   delta + theta) * (1 / 2 * delta + theta)^(-(d / 2 +
                                                                                 1))
    pdfconst <- (2 * pi)^(-d / 2) * theta * exp(theta) * det(Sigma)^(-1 / 2)
    PDF <- pdfconst * pdfgamma
  }
  if (formula == "indirect") {
    delta <- sapply(1:nrow(x), function(i) t(as.vector(t(x[i, ]) - mu)) %*% solve(Sigma) %*% as.vector(t(x[i, ]) - mu))
    intf <- function(w, gamm) {
      w^(d / 2) * exp(-w * gamm)
    }
    pdfinteg <- sapply(1:nrow(x), function(i) {
      stats::integrate(intf,
                       lower = 1, upper = Inf, gamm = delta[i] / 2 + theta
      )$value
    })
    pdfconst <- (2 * pi)^(-d / 2) * theta * exp(theta) * det(Sigma)^(-1 / 2)
    PDF <- pdfconst * pdfinteg
  }
  if (formula == "series") {
    delta <- sapply(1:nrow(x), function(i) t(as.vector(t(x[i, ]) - mu)) %*% solve(Sigma) %*% as.vector(t(x[i, ]) - mu))
    delta <- replace(delta, delta == 0, 1 / (theta * (2 *
                                                        pi)^(d / 2) * (d / 2 + 1)) * (1 - (1 - theta)^(d / 2 +
                                                                                                         1)))
    n <- d / 2
    term <- sapply(1:length(delta), function(j) {
      exp(-delta[j] / 2 -
            theta) * (delta[j] / 2 + theta)^(-1) * (1 + sum(sapply(
              1:floor(n),
              function(i) {
                prod(seq(from = n, to = n - i + 1, by = -1)) *
                  (delta[j] / 2 + theta)^(-i)
              }
            )))
    })
    if (d %% 2 == 1) {
      term <- term + sapply(1:length(delta), function(j) {
        prod(seq(
          from = n,
          to = 0.5, by = -1
        )) * sqrt(pi) * 2 * 1 / (delta[j] / 2 +
                                   theta)^(floor(n) + 1 + 1 / 2) * (1 - stats::pnorm(sqrt(2) *
                                                                                       sqrt(delta[j] / 2 + theta)))
      })
    }
    PDF <- (2 * pi)^(-d / 2) * det(Sigma)^(-1 / 2) * theta *
      exp(theta) * term
  }
  return(PDF)
}

msen.WMM <- function(X, w = NULL, Rfun = "optimize"){
  
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
  
  v <- function(theta){
    theta*exp(theta)*expint::gammainc(a = 0, x = theta)
  }
  
  # kurtosis factor
  
  k <- function(theta){
    (1 - theta*exp(theta)*expint::gammainc(a = 0, x = theta))/(theta*exp(2*theta)*expint::gammainc(a = 0, x = theta)^2)
  }
  
  # objective function
  
  ktheta <- Kurt/(d*(d+2))
  if(ktheta == 1){
    
    warning("The best MSEN distribution is the normal one; for approximation, theta is fixed to 100.")
    theta <- 100
    
  }
  if(ktheta > 1){
    
    if(Rfun == "optimize"){
      
      fobj <- function(theta, ktheta){
        
        L2dist <- abs(k(theta = theta) - ktheta)^2
        
        return(L2dist)
        
      }
      
      res   <- stats::optimize(f = fobj, ktheta = ktheta, lower = 0, upper = 100) # , control=list(maxit = maxit, reltol = reltol, trace = 1)
      theta <- res$minimum
      
    }
    if(Rfun == "optim"){
      
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
##################################### Convergence ##################################

getall <- function(loglik) {
  if (length(loglik) <3){ val=1
  }else{
  n = length(loglik)
  lm1 = loglik[n]
  lm  = loglik[(n-1)]
  lm_1  = loglik[(n-2)]
  am = (lm1 - lm)/(lm - lm_1)
  lm1.Inf = lm + (lm1 - lm)/(1-am)
  val = lm1.Inf - lm
  if (is.nan(val)) val=0
  if (val < 0) val= 1}
  return( val )
}
