library(mgcv)
library(dplyr)
library(rgl)
library(manipulate)
source("R/smooth.construct.sps.smooth.spec.R")

# generate toy data
dat <- data.frame( x = 1:100 )
ps_obj <- with(dat, s(x, bs = "ps"))
B <- Predict.matrix(smooth.construct(ps_obj, dat, NULL), dat)               
set.seed(3904)
dat$y <- B %*% rnorm(ncol(B))
plot(dat, t = "l")

mod0 <- gam( y ~ s(x, bs = "sps", xt = list(skew = TRUE)), dat = dat )
lines(dat$x, predict(mod0), col = "cornflowerblue", lty = "dashed")

mod1 <- gam( y ~ s(x, bs = "sps"), dat = dat[1:50, ])
lines(dat[1:50, ]$x, predict(mod1), col = "darkred", lty = "dashed")

mod2 <- gam( y ~ s(x, bs = "sps"), dat = dat, knots = list(x = c(-50, 100)))
lines(dat$x, predict(mod2), col = "chocolate", lty = "dashed")


# make tensor product spline which is skewsym in arg1 ---------------------


Ba <- Predict.matrix(smooth.construct(s(arg1, bs = "sps", xt = list(skew = TRUE)), dat1, NULL), dat1)
B1 <- Predict.matrix(smooth.construct(s(arg1, bs = "ps"), dat1, NULL), dat1)
set.seed(9348885)
dat2$y <- c(Ba %*% tcrossprod( matrix(rnorm(ncol(Ba)*ncol(B1)), ncol = ncol(B1)), B1))

# fit data
moda <- gam( y ~ -1 + te(arg1, arg2, bs = c("sps", "ps"), xt = list(skew = TRUE), fx = TRUE), data = dat2)
dat2$pred <- predict(moda)

persp3d(dat1$arg1, dat1$arg1, dat2$y, col = "cornflowerblue", alpha = .5)
persp3d(dat1$arg1, dat1$arg1, dat2$pred, col = "darkred", add = TRUE)

# -> seems to over-smooth ??? somehow it doesn't work. 
# Maybe the basis coef don't fit?


# check out point symmetric smooths ---------------------------------------

dat1 <- data.frame(arg1 = 1:50)
dat2 <- expand.grid(arg1 = 1:50, arg2 = 1:50)

Bskew <- Predict.matrix(
  smooth.construct( 
    s(arg1, arg2, bs = "sps", xt = list(skew = TRUE)),
    data = dat2, knots = NULL ),
  data = dat2 )
Bsymm <- Predict.matrix(
  smooth.construct( 
    s(arg1, arg2, bs = "sps", xt = list(skew = FALSE)),
    data = dat2, knots = NULL ),
  data = dat2 )

set.seed(934811)
dat2$yskew <- c(Bskew %*% rnorm(ncol(Bskew)))
dat2$ysymm <- c(Bsymm %*% rnorm(ncol(Bsymm)))

modpa <- gam( I(yskew + ysymm) ~ s(arg1, arg2, bs = "sps", xt = list(skew = TRUE)) + 
                s(arg1, arg2, bs = "sps", xt = list(skew = FALSE)), data = dat2)
preds <- predict(modpa, type = "terms")
dat2$predskew <- preds[,1]
dat2$predsymm <- preds[,2]

# skew-symm part
persp3d(dat1$arg1, dat1$arg1, dat2$yskew, col = "cornflowerblue", alpha = .5)
persp3d(dat1$arg1, dat1$arg1, dat2$predskew, col = "darkred", add = TRUE)

# symm part (intercept missing)
persp3d(dat1$arg1, dat1$arg1, dat2$ysymm, col = "cornflowerblue", alpha = .5)
persp3d(dat1$arg1, dat1$arg1, dat2$predsymm, col = "darkred", add = TRUE)

# => perfect !

# cyclic (skew-)symmetric splines ---------------------------------------

B <- Predict.matrix(
  smooth.construct( 
    s(arg1, arg2, bs = "sps", xt = list(cyclic = TRUE)),
    data = dat2, NULL ),
  dat2)

# example surface
set.seed(8495)
dat2$y <- B %*% rnorm(ncol(B))
y <- matrix(dat2$y, ncol = sqrt(nrow(dat2)))
persp3d(dat1$arg1,dat1$arg1,y, col = "chocolate")
#check margins
matplot(dat1$arg1, y[, 10], t = "l")
abline( h = y[1, 10], col = "grey")

Ba <- Predict.matrix(
  smooth.construct( 
    s(arg1, arg2, bs = "sps", xt = list(cyclic = TRUE, skew = TRUE)),
    data = dat2, NULL ),
  dat2)

C <- make_summation_matrix(10, skew = TRUE, cyclic.degree = 1)
C <- array(C, dim = c(rep(sqrt(nrow(C)), 2), ncol(C)))
image(Matrix(C[,,1]))

# skew-symmetric example surface
set.seed(8333)
y <- matrix(Ba %*% rnorm(ncol(Ba)), ncol = sqrt(nrow(dat2)))
persp3d(dat1$arg1,dat1$arg1,y, col = "chocolate")
# check skew-symmetry
image(dat1$arg1,dat1$arg1,y)
#check margins
matplot(dat1$arg1, y[, 10], t = "l")
abline( h = y[1, 10], col = "grey")

# NICE !!!!


