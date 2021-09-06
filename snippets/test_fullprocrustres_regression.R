
# example dat
dat <- fdasrvf::beta[,,4,]
dat <- aperm(dat, c(2,1,3))

pre <- preshape(dat)
pre <- as_complex(pre)

C <- tcrossprod(Conj(pre), pre)
R1 <- diag(nrow = 13)
R2 <- crossprod(diff(R1)) + R1 # first + zero order difference mat
R <- array( c(R1, R2), dim = c(13,13,2))

C <- Re(C)

# explicit atempt
M <- solve(R1) %*% C + solve(R2) %*% C
m <- eigen(M)$vectors[,1]
plot(m, t = "l")

loss <- function(Theta) {
  enum <- crossprod(Theta, C) %*% Theta
  enum / crossprod(Theta, R1) %*% Theta + 
    enum / crossprod(Theta, R2) %*% Theta +
    1000 * abs(sum(Theta^2) - 1) 
}

# numeric attempt
m2 <- optim( par = m, fn = loss)$par
lines(m2, col = "red", t = "l", lty = "dashed")

