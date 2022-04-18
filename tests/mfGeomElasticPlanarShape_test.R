library(manifoldboost)
library(Momocs)
data(bot, package = "Momocs")

b <- lapply(bot[["coo"]], as.data.frame)
b <- lapply(names(b), function(x) {
  b[[x]]$t <- (1:nrow(b[[x]]))/nrow(b[[x]])
  b[[x]]$id <- x
  b[[x]]
  })
names(b) <- names(bot[["coo"]])
# regularize
t0 <- seq(0,1, len = 101)[-101]
b <- lapply(b, function(x) {
  data.frame(
    V1 = approx(x$t, x$V1, t0)$y[c(length(t0), 2:length(t0))],
    V2 = approx(x$t, x$V2, t0)$y[c(length(t0), 2:length(t0))],
    t = t0,
    id = x$id[1]
  )
})
# format
b <- lapply(b, function(x) {
  d <- x[rep(1:nrow(x), 2), 3:4]
  d$dim <- factor(rep(c("x", "y"), each = nrow(x)))
  d$value <- c(x$V1, x$V2)
  d
})


# check single geometry ---------------------------------------------------

w <- mfGeomElasticPlanarShape_closed$new(b$franziskaner, value^dim ~ t|id #,warp_update = function(warp, ...) !warp
                                        )
w$plot()
w2 <- mfGeomElasticPlanarShape_closed$new()
w2$initialize(b$ballantines, value^dim ~ t|id)

w$pole_ <- w2$y_
arg0 <- attr(w$y_, "arg")

s <- mfGeomPlanarShape$new(b$franziskaner, value^dim ~ t|id)
s2 <- mfGeomPlanarShape$new(b$ballantines, value^dim ~ t|id)
s$pole_ <- s2$y_

# check align
.y_ <- w$y_
y_aligned <- w$align(y0_ = w$pole_)
str(list(before_alignment = .y_, after_alignment = w$y_))

plot(c(0,1), c(0,1), t = "l", asp = 1, col = "cornflowerblue")
lines(arg0, attr(w$y_, "arg_optim"), t = "l", asp = 1)
# check log
franz <- s$y_
attr(franz, "arg") <- arg0
identical(attr(w2$y_, "arg")[1:length(arg0)], arg0)
y_v <- w2$log(y0_ = w2$srv_trafo(franz))
identical(attr(w2$y_, "arg_optim")[1:length(arg0)], arg0)
lines(arg0, attr(w2$y_, "arg_optim"), col = "darkred")

par(mfrow = c(1,3))
s$plot(t = "l", main = "without warping alignment")
w$plot(t = "l", main = "with warping alignment")
w2$plot(y0_ = w2$srv_trafo(franz), t = "l", main = "with warping alignment")


# check product geometry --------------------------------------------------

mf <- mfGeomProduct$new(
  mfGeom_default = mfGeomElasticPlanarShape_closed$new(), 
  data = dplyr::bind_rows(b), formula = value^dim ~ t|id)

par(mfrow = c(5,5), mar = c(0,0,2,0))
mf$slice(which = 1:25)
mf$pole_ <- mf$y_[rep("ballantines", length(mf$pole_))]
system.time(
  mf$plot(t = "l")
)

# check model fit ---------------------------------------------------------

library(cubelyr)
bcube <- tbl_cube(
  dimensions = list(arg = t0, dim = c("x", "y"), id = names(b)),
  measures = list(value = array(sapply(b, `[[`, "value"), dim = c(length(t0), 2, length(b))))
)

bdat <- list(shape = bcube, type = factor(bot[["fac"]]$type, labels = c("whisky", "beer")))

fam <- ElasticPlanarShapeL2( )

system.time(
  m <- mfboost(shape ~ bols(type, df = Inf),
               obj.formula = value^dim ~ bbs(arg, df = Inf, knots = 50, 
                                             boundary.knots = c(0,1), cyclic = TRUE) | id,
               data = bdat,
               family = fam,
               control = boost_control(mstop = 5, nu = .2))
)

# pdf("First_PlanarShapeBoost_withWarping.pdf")
par(mfrow = c(1,2), mar = c(0,2,2,0))
plot(m, ids = which(names(b) %in% c("franziskaner", "ballantines")), 
     t = "l", main = c("beer", "whisky"))
par(mfrow = c(1,1))
panel(bot, names = TRUE, fac = "type")
# dev.off()

p <- predict(m, type = "response")
p <- fam@mf$structure(p)
for(i in seq_along(p))
  attr(p[[i]], "arg") <- c(t0, 1) 
p <- lapply(p, fam@mf$mfGeom_default$srv_trafo, inverse = T, center = T)

plot(p$franziskaner, t = "l", asp = 1)
plot(p$jackdaniels, t = "l", asp = 1)
