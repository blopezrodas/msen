
##################################### TEST ##################################
library(MixtureMissing)
library(mice)

############### Test 1 ###############
set.seed(10)
X1 = rmsen(200, mu = c(10,10,10), Sigma = diag(3), theta = 5)$X
patterns_matr <- generate_patterns(3)
X_m1 <- ampute(X1, prop = 0.5, patterns = patterns_matr)$amp
############### Test 2 ###############
X2 = rmsen(200, mu = c(10,10,10), Sigma = diag(3), theta = 0.5)$X
patterns_matr <- generate_patterns(3)
X_m2 <- ampute(X2, prop = 0.5, patterns = patterns_matr)$amp
############### Test 3 ###############
sig = matrix(0.6, 4, 4)
diag(sig) = 1
X3 = rmsen(200, mu = c(10, 10, 10, 3), Sigma = sig, theta = 12)$X
patterns_matr <- generate_patterns(4)
X_m3 <- ampute(X3, prop = 0.5, patterns = patterns_matr)$amp
############### Test 4 ###############
n <- c(500, 350)
mu <- matrix(c(10, 5), ncol = 2)
mu <- rbind(mu, c(15, 10))
Sigma <- array(, dim = c(2, 2, 2))
Sigma[, , 1] <- matrix(c(10, 5, 5, 9), ncol = 2)
Sigma[, , 2] <- matrix(c(4, -2, -2, 4), ncol = 2)
theta <- c(6, 2)
data <- c()
for (i in 1:2) {
  data <- rbind(data, rmsen(n = n[i], mu = mu[i,], Sigma = Sigma[,,i], theta = theta[i])$X)
}
row.names(data) <- NULL
X_m4 <- MixtureMissing::hide_values(data, 0.2)

#############################################
res1 = SEN_m(X_m1) # 8 w/ e = 0.0001
res2 = SEN_m(X_m2) # 200 w/ e = 0.0001
res22 = SEN_m(X_m2, epsilon = 0.1) # 3 w/ e = 0.1, 164 W/ e = 0.01
res3 = SEN_m(X_m3) # 200 w/ e = 0.0001
res4 = SEN_m(X_m4) # 3 w/ e = 0.0001

#############################################
res1 = SEN_m(X_m1, init_method = "kmeans") # 200 w/ e = 0.0001
res2 = SEN_m(X_m2, init_method = "kmeans") # 12 w/ e = 0.0001 RUN THIS ONE !!!!
res22 = SEN_m(X_m2, epsilon = 0.1, init_method = "kmeans") # 7 w/ e = 0.1, 9 W/ E = 0.01
res3 = SEN_m(X_m3, init_method = "kmeans") # 200 w/ e = 0.0001
res4 = SEN_m(X_m4, init_method = "kmeans") # 200 w/ e = 0.0001
