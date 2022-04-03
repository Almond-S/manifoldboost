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
t0 <- seq(0,1, len = 100)
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

w <- mfGeomWarpPlanarShape$new(b$franziskaner, value^dim ~ t|id)
w$plot(t = "l")
w2 <- mfGeomWarpPlanarShape$new(b$ballantines, value^dim ~ t|id)
w$pole_ <- w2$y_
arg0 <- attr(w$y_, "arg")

s <- mfGeomPlanarShape$new(b$franziskaner, value^dim ~ t|id)
s2 <- mfGeomPlanarShape$new(b$ballantines, value^dim ~ t|id)
s$pole_ <- s2$y_

# check align 
y_aligned <- w$align(y0_ = w$pole_)
plot(c(0,1), c(0,1), t = "l", asp = 1, col = "cornflowerblue")
lines(arg0, attr(w$y_, "arg"), t = "l", asp = 1)
# check log
franz <- s$y_
attr(franz, "arg") <- arg0
identical(attr(w2$y_, "arg"), arg0)
y_v <- w2$log(y0_ = franz)
identical(attr(w2$y_, "arg"), arg0)
lines(arg0, attr(w2$y_, "arg"), col = "darkred")

par(mfrow = c(1,2))
s$plot(t = "l", main = "without warping alignment")
w$plot(t = "l", main = "with warping alignment")


# check model fit ---------------------------------------------------------

library(cubelyr)
bcube <- tbl_cube(
  dimensions = list(arg = t0, dim = c("x", "y"), id = names(b)),
  measures = list(value = array(sapply(b, `[[`, "value"), dim = c(length(t0), 2, length(b))))
)

bdat <- list(shape = bcube, type = factor(bot[["fac"]]$type, labels = c("whisky", "beer")))

fam <- manifoldboost:::WarpPlanarShapeL2()

system.time(
  m <- mfboost(shape ~ bols(type, df = Inf), 
             obj.formula = value^dim ~ bbs(arg, df = Inf, knots = 50, cyclic = TRUE) | id, 
             data = bdat, 
             family = fam, 
             control = boost_control(mstop = 3, nu = .25))
  )
pdf("First_PlanarShapeBoost_withWarping.pdf")
par(mfrow = c(1,2), mar = c(0,2,2,0))
plot(m, ids = which(names(b) %in% c("franziskaner", "ballantines")), t = "l", main = c("beer", "whisky"))
warnpar(mfrow = c(1,1))
panel(bot, names = TRUE, fac = "type")
dev.off()

plot(t0, attr(fam@mf$y_$ballantines, "arg"), t = "l")
