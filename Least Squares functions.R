x1 <- rnorm(100, mean = 20, sd = 5)
x2 <- rpois(100, lambda = 10)
x3 <- rgamma(100, shape = 10, scale = 10)

beta <- c(0.2, 0.3, 0.5)
X <- as.matrix(cbind(x1, x2,x3), ncol = 3)
Y <- X%*%beta+rnorm(100)



# constrained least squares
cls <- function(Y, X){
    one_vec <- rep(1,ncol(X))
    lambda <- (t(one_vec)%*%solve(t(X)%*%X)%*%t(X)%*%Y-1)/(t(one_vec)%*%solve(t(X)%*%X)%*%one_vec)
    betahat <- solve(t(X)%*%X)%*%(t(X)%*%Y-c(lambda)*one_vec)
    return(betahat)
}
cls(Y, X)
    