library(ggplot2)
library(cubelyr)
library(formula.tools)
library(manifoldboost)

# generate toy shapes
arg <- 1:20 / 20
id <- letters[1:8]
n_arg <- length(arg)
n_id <- length(id)
x <- (rep(seq_along(id),each=n_arg)/5 + 1) * cos(arg * 2*pi) + 2*cos(arg * 4*pi)
y <- (rep(seq_along(id),each=n_arg)/5 + 1) * sin(arg * 2*pi) + 2*sin(arg * 4*pi)

# set up data in very long format
dat_FDboost <- data.frame(
  wert = c(x,y), 
  d = rep(c("x","y"), each = n_arg*n_id),
  time = rep(arg, 2*length(id)), 
  ID = rep(rep(id, each = n_arg), 2)
)
form <- wert^d ~ time | ID

plot(x,y, t = "b", col = factor(rep(id, each = n_arg)))


# check one shape only ----------------------------------------------------

i <- "a"
# for(i in 1:n_id) 
{
# consider one single shape
s <- mfGeomPlanarShape$new()
dat1 <- dat_FDboost[dat_FDboost$ID == i, ]

s$initialize(dat1, form)
# circumvent direct comparison with y_ because of $register
s$pole_ <- s$structure(dat1$wert)
all(dat1$wert == s$unstructure(s$pole_))

ggplot(dat_FDboost, aes(x = time, y = wert)) + 
  geom_line() + 
  geom_point() +
  facet_wrap(~d)

# now shuffle data.frame order and try again
set.seed(9243)
shuffle <- sample(1:nrow(dat1))

ggplot(dat1[shuffle,], aes(x = time, y = wert)) + 
  geom_line() + 
  geom_point() + 
  facet_wrap(~d)

s$initialize(dat1[shuffle,], form)
s$pole_ <- s$structure(dat1[shuffle,]$wert)
plot(s$pole_)

dat1$wert2[shuffle] <- s$unstructure(s$pole_)

ggplot(dat1, aes(x = time, y = wert)) + 
  geom_line() +
  geom_line(aes(y = wert2, col = "retransform")) +
  geom_point() +
  facet_wrap(~d)
all(dat1$wert == dat1$wert2)
# nice!

# check normal vector
nvec <- s$get_normal(s$pole_)
stopifnot(
  crossprod(nvec[,1], s$unstructure(s$pole_))^2 == 
    crossprod(nvec[,1]) * crossprod(s$unstructure(s$pole_))
  )
}

# check full irregular data set -------------------------------------------

g <- mfGeomProduct$new(mfGeom_default = s)

g$initialize(data = dat_FDboost, formula = wert^d ~ time | ID)
# again circumvent registration of y_
g$pole_ <- g$structure(dat_FDboost$wert)
all(dat_FDboost$wert == g$unstructure(g$pole_))

dat_FDboost$wert2 <- g$unstructure(g$pole_)

ggplot(dat_FDboost, aes(x = time, y = wert)) + 
  geom_line() +
  geom_line(aes(y = wert2, col = "retransform")) +
  geom_point() +
  facet_wrap(~d)

plot(g$y_[[1]], t="l")
for(i in 2:length(g$y_)) lines(g$y_[[i]], col = i)

# check distance fun argument
all.equal(
  g$distance(g$y_, rev(g$y_)),
  sqrt(g$distance(g$y_, rev(g$y_), squared = TRUE))
)

# check normal vectors (at least the first)
nvec <- g$get_normal(g$pole_)
stopifnot(
  crossprod(nvec[,1], g$unstructure(g$pole_))^2 == 
    crossprod(nvec[,1]) * crossprod(g$unstructure(g$pole_))
)

# try shuffling
shuffle <- sample(1:nrow(dat_FDboost))
g$initialize(data = dat_FDboost[shuffle, ], formula = wert^d ~ time | ID)
g$pole_ <- g$structure(dat_FDboost[shuffle, "wert"])
stopifnot(
  all(dat_FDboost$wert[shuffle] == g$unstructure(g$pole_)))
plot(g$y_[[1]], t="l")
for(i in 2:length(g$y_)) lines(g$y_[[i]], col = i)
# nice!

# register manually
a <- g$register(g$pole_)
stopifnot(identical(a, g$y_))

# check normal vectors (at least the first)
nvec <- g$get_normal(g$pole_)
stopifnot(
  crossprod(nvec[,1], g$unstructure(g$pole_))^2 == 
    crossprod(nvec[,1]) * crossprod(g$unstructure(g$pole_))
)

# check cloning -----------------------------------------------------------

g2 <- g$clone(deep = TRUE)
.y_ <- g2$.__enclos_env__$private$.y_
stopifnot(!identical(g$.__enclos_env__$private$.y_, .y_))
g2$y_ <- g$pole_
g$clone(deep = TRUE)$slice(1)$plot()
g2$clone(deep = TRUE)$slice(1)$plot(y0_par = list(type = "l"))

# check slicing -----------------------------------------------------------

ids <- letters[c(7, 1, 4)]
dat_FDboost2 <- dat_FDboost[shuffle, ]
dat_FDboost2 <- dat_FDboost2[dat_FDboost2$ID %in% ids, ]
g$slice(ids)
stopifnot(length(g$y_) == length(ids))
g$plot(g$structure(dat_FDboost2$wert))
stopifnot(all.equal(g$structure(dat_FDboost2$wert), g$pole_))
stopifnot(
  all(dat_FDboost2$wert == g$unstructure(g$pole_)))
plot(g$y_[[1]], t="l")
for(i in 2:length(g$y_)) lines(g$y_[[i]], col = i)
# nice!

# try regular format ------------------------------------------------------

cube <- as.tbl_cube(dat_FDboost, dim_names = c("d", "time", "ID"))
cube_FDboost <- as_FD(cube, obj.formula = wert^d ~ time|ID)

g <- mfGeomProduct$new(mfGeom_default = s)
g$initialize(data = cube_FDboost, 
             formula = wert^d ~ time | ID)
plot(g$y_[[1]], t="l")
for(i in 2:length(g$y_)) lines(g$y_[[i]], col = i)
# => initialization seems to work
g$pole_ <- g$structure(c(cube_FDboost$wert))
stopifnot(all.equal(as.vector(g$unstructure(g$pole_)), c(cube_FDboost$wert)))
# nice!

# register manually
a <- g$register(g$pole_)
stopifnot(identical(a, g$y_))

nvec <- g$get_normal(g$pole_)
nvec <- g$get_normal(g$pole_[rep(1, length(g$pole_))])


# check $y_ <-  -----------------------------------------------------------

g$y_ <- g$register(g$pole_)










