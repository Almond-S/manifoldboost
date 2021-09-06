library(manifoldboost)
library(ggplot2)
library(Matrix)

# generate toy shapes
arg <- 1:20 / 20
id <- 1:8
n_arg <- length(arg)
n_id <- length(id)
x <- (rep(id,each=n_arg)/5 + 1) * cos(arg * 2*pi) + 2*cos(arg * 4*pi)
y <- (rep(id,each=n_arg)/5 + 1) * sin(arg * 2*pi) + 2*sin(arg * 4*pi)

# set up data in very long format
dat_FDboost <- data.frame(
  wert = c(x,y), 
  d = rep(c("x","y"), each = n_arg*n_id),
  time = rep(arg, 2*length(id)), 
  ID = rep(rep(id, each = n_arg), 2)
)
form <- wert^d ~ time | ID


# initialize planar shape family ------------------------------------------

fam <- PlanarShapeL2()

set.seed(39408)
dat_FDboost <- dat_FDboost[sample(1L:nrow(dat_FDboost)), ]

fam@mf$initialize(data = dat_FDboost, formula = form)
# is regular, so just take first for simplicity
fam@mf$pole_ <- fam@mf$y_[rep(1, length(fam@mf$y_))]

n <- fam@mf$get_normal(fam@mf$pole_)

# build base-learner ------------------------------------------------------

dat_FDboost$d <- factor(dat_FDboost$d)
bl <- with(dat_FDboost, brandom(d, df = Inf) %X% bbs(time))

# build and extract tangent function
tform <- fam@update_formula(form, fam@mf$pole_)
tangent <- environment(tform)$tangent

blt <- tangent(bl)
# get design matrix
X <- extract(blt, "design")
A <- crossprod(n, X)
image(A)
summary(A@x)
# => seems to work in principle

# check whether n is really normal to the unstructured pole
crossprod(n[,1], fam@mf$unstructure(fam@mf$pole_))

