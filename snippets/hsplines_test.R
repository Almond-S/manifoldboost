source("R/hsplines.R")

x <- seq(0,1, len = 1000)


# no constraint -----------------------------------------------------------

B <- hspline(x)
matplot(x, B, t = "l")


# periodicity constraint --------------------------------------------------

B <- hspline(x, periodic = TRUE)
matplot(x, B, t = "l") # note that colors match
set.seed(1808)
# three random periodic B-splines
matplot(B %*% matrix(rnorm(ncol(B)*3), ncol = 3), t = "l")


# no constant function ----------------------------------------------------

B <- hspline(x, remove.constant = TRUE)
matplot(x, B, t = "l")
abline(h = 0, col = "grey")
# check basis is orthogonal to intercept
summary( crossprod(B, rep(1, length(x))) / length(x) )
# note that this does not have to mean that on an instance of data
# the sum of each the basis functions has be precisely zero.
# other than done elsewhere, we choose a transform such that
# the function basis is actually orhtogonal to the constant function.


# exclude constant in periodic basis --------------------------------------

B <- hspline(x, remove.constant = TRUE, periodic = TRUE)
matplot(x, B, t = "l")
abline(h = 0, col = "grey")


# check out derivatives ---------------------------------------------------

dB <- hspline(x, remove.constant = TRUE, periodic = TRUE, derivs = 1)
matplot(x, dB, t = "l")
dB2 <- hspline(x, remove.constant = FALSE, periodic = TRUE, derivs = 1)
# derivatives of the Bsplines with and without constant have the same rank,
# as the derivative of the constant is zero.
stopifnot( qr(dB2)$rank == qr(dB)$rank )


# check sparse design matrices --------------------------------------------

B <- hspline(x, remove.constant = TRUE, periodic = TRUE, sparse = TRUE)
# works but doesn't make so much sense as the constraint basis won't be
# sparse anymore.


# check length(x) == 1 special case ---------------------------------------

B <- hspline(x, remove.constant = TRUE, periodic = TRUE)
b <- hspline(.5, boundary.knots = c(0,1), remove.constant = TRUE, periodic = TRUE)
matplot(x, as.matrix(B), t = "l")
points(rep(.5, ncol(b)), c(b), cex = 2, pch = 4)


# try with not equally spaces knots ---------------------------------------

set.seed(8734)
# draw inner knots
knots <- cumsum(rgamma(10, shape = 1))
knots <- knots[1:9] / tail(knots, 1)
B <- hspline(x, knots = knots, periodic = TRUE)
matplot(x, as.matrix(B), t = "l")
plot(B %*% rep(1, ncol(B)), t = "l")
# constant coef <=> constant function is still valid
# => also removing the constant still works
B <- hspline(x, knots = knots, periodic = TRUE, remove.constant = TRUE)
matplot(x, as.matrix(B), t = "l")
# check rank of first derivative
dB <- hspline(x, knots = knots, derivs = 1, 
              periodic = TRUE, remove.constant = TRUE)
stopifnot(qr(dB)$rank == qr(B)$rank )
