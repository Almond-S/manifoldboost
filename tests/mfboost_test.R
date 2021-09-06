
library(manifoldboost)
library(dplyr)
library(cubelyr)

dat <- shapes::digit3.dat

# generate size as covariates
x_size <- apply(dat, 3, function(x) var(c(x)))

# prepare regular dataset
digitnames <- paste0("digit_nr", seq_len(dim(dat)[3]))
data <- list(
  y = tbl_cube(
    dimensions = list(arg = 1:13, dim = c("x", "y"), id = digitnames),
    measures = list(value = dat)
  ),
  x = x_size
)

m0 <- mfboost(y ~ 1, 
              obj.formula = value^dim ~ bbs(arg) | id, 
              data = data, family = PlanarShapeL2(pole.type = "Gaussian"))
plot(m0, t = "b", ids = 1)

# predict on dataset

pred0 <- predict(m0, type = "response")

data_grid <- expand.grid(
  value = NA,
  arg = seq(0,1, len = 100), 
  dim = c("x", "y"), 
  id = c("pole"))

# e <- environment(tangent)
# nv <- as.data.frame(e$normal_vecs)
# data <- c(data, nv)

pred0 <- predict(m0, type = "response", newdata = data)
pred0_ <- m0$family@mf$structure(pred0)
plot(pred0_[[1]], t = "b", asp = 1)

#################################

m <- mfboost(formula = y ~ bbs(x, knots = 3), 
             obj.formula = value^dim ~ bbs(arg) | id, 
             data = data, 
             family = PlanarShapeL2())

X <- extract(m$baselearner[[1]], "design")
sapply(X, dim)
B <- as.matrix(X[[2]])
pole <- m$family@mf$pole_
pole <- manifoldboost:::complex2realmat(pole[[1]])
crossprod(B,pole)

image(B)

pred <- predict(m, type = "response", newdata = data)
pred_ <- m$family@mf$structure(pred)
plot(pred_[[1]], t = "b", asp = 1)
for(i in 2:length(pred_))
  lines(pred_[[i]], t = "l", col = "grey")

##################################

# check out factorization on 'new' data
# fac <- factorize(m, newdata = data)
