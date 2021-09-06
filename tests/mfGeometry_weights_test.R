library(manifoldboost)
library(ggplot2)
library(cubelyr)
library(formula.tools)
library(dplyr)

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

plot(x,y, t = "b", col = factor(rep(id, each = n_arg)))


# check weights for one shape only ----------------------------------------

dat1 <- dat_FDboost[dat_FDboost$ID == 1, ]
s <- mfGeomPlanarShape$new(dat1, form, weight_fun = trapez_weights)
s$weights_
s$unstructure_weights(s$weights_)
stopifnot(abs(s$innerprod(s$y_) - 1) < 1e-15)

# default does not over-write weight_fun
s$initialize(dat1, form) # default 
is.null(s$weights_)

# check delayed data initialization
s <- mfGeomPlanarShape$new(weight_fun = trapez_weights)
s$initialize(dat1, form, weight_fun = NULL)
is.null(s$weights_)

# and when only passing NULL
s <- mfGeomPlanarShape$new(dat1, form, weight_fun = NULL)
# => works as well


# now, check multiple shapes -----------------------------------------------

g <- mfGeomProduct$new(mfGeom_default = s)
g$initialize(data = dat_FDboost, formula = form)
plot(g$y_[[1]], t = "b")
for(i in seq_along(g$y_)) lines(g$y_[[i]], t = "b", col = i)

unlist(g$innerprod(g$y_))
# compare with un-weighted inner product
sapply(g$y_, function(x) crossprod(Conj(x), x))

# try irregular version
set.seed(39082)
dat_FDboost <- dat_FDboost %>% group_by(ID, time) %>% 
  mutate(select = rbinom(1,1,.9)) %>% ungroup()
dat_FDboost_sub <- dat_FDboost %>% filter(select == 1)

g$initialize(data = dat_FDboost_sub, formula = form)
plot(g$y_[[1]], t = "b")
for(i in seq_along(g$y_)) lines(g$y_[[i]], t = "b", col = i)
stopifnot(all(Mod(unlist(g$innerprod(g$y_)) - 1) < 1e-15))

sapply( g$weights_, function(x) c(len = length(x), sum = sum(x)) )

# now with arg_range specified 
s <- mfGeomPlanarShape$new(weight_fun = trapez_weights, arg_range = c(0,1))
g <- mfGeomProduct$new(mfGeom_default = s, data = dat_FDboost_sub, formula = form)

sapply( g$weights_, function(x) c(len = length(x), sum = sum(x)) )
unlist(g$innerprod(g$y_))

