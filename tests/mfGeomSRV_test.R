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

w <- mfGeomSRV_closed$new(b$franziskaner, value^dim ~ t|id)
w2 <- mfGeomSRV_closed$new()
w2$initialize(b$ballantines, value^dim ~ t|id)
w$pole_ <- w2$y_
w$plot()
w$plot( y_ = structure(w2$y_, arg = NULL), y0_ = structure(w$y_, arg = NULL) )
w$plot(level = "SRV")

# check product geometry --------------------------------------------------

mf <- mfGeomProduct$new(
  mfGeom_default = mfGeomSRV_closed$new(), 
  data = dplyr::bind_rows(b), formula = value^dim ~ t|id)

par(mfrow = c(5,5), mar = c(0,0,2,0))
mf$slice(which = 1:25)
mf$pole_ <- mf$y_[rep("ballantines", length(mf$pole_))]
system.time(
  mf$plot()
)
mf$plot(level = "SRV")

mf2 <- mfGeomProduct$new(mfGeom_default = mfGeomSRV_closed$new())
mf2$initialize(data = dplyr::bind_rows(b), formula = value^dim ~ t|id)
mf2$plot()
mf2$plot(level = "SRV")

# check model fit ---------------------------------------------------------

library(cubelyr)
bcube <- tbl_cube(
  dimensions = list(arg = t0, dim = c("x", "y"), id = names(b)),
  measures = list(value = array(sapply(b, `[[`, "value"), dim = c(length(t0), 2, length(b))))
)

bdat <- list(shape = bcube, type = factor(bot[["fac"]]$type, labels = c("whisky", "beer")))

FDdat <- as_FD(bdat, formula = shape ~ bols(type, df = Inf),
               obj.formula = value^dim ~ bbs(arg, df = Inf, knots = 50, cyclic = TRUE) | id)
mf3 <- mfGeomProduct$new(mfGeom_default = mfGeomSRV_closed$new())
mf3$initialize(data = FDdat, formula = value^dim ~ arg|id)
mf3$plot()
mf3$plot(level = "SRV")


fam <- SquareRootVelocityL2()

system.time(
  m <- mfboost(shape ~ bols(type, df = Inf),
               obj.formula = value^dim ~ bbs(arg, df = Inf, knots = 50, cyclic = TRUE) | id,
               data = bdat,
               family = fam)
)

par(mfrow = c(1,2), mar = c(0,2,2,0))
plot(m, ids = which(names(b) %in% c("franziskaner", "ballantines")), t = "l", main = c("beer", "whisky"), seg_par = NA)

p <- m$family@mf$structure(predict(m))
plot(w$.__enclos_env__$private$close(w$srv_trafo(p$franziskaner, inverse = TRUE)), t = "l", asp = 1)
lines(w$.__enclos_env__$private$close(w$srv_trafo(p$jackdaniels, inverse = TRUE)), t = "l", asp = 1, col = "darkseagreen")
