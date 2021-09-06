library(manifoldboost)
library(formula.tools)

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


# check geometry -----------------------------------------------------

dat1 <- dat_FDboost[dat_FDboost$ID == 1, ]
shp <- mfGeomPlanarSizeShape$new(data = dat1, form = form)
shp$plot()

shps <- mfGeomProduct$new(shp, dat_FDboost, form = form)
# re-scale and rotate shapes
set.seed(8823)
sr <- complex(re = rnorm(length(id)), im = rnorm(length(id)))
shps$y_ <- Map(`*`, shps$y_, sr)
par(mfrow=c(3,3), mar = c(0,2,2,0))
shps$plot(t="b", ylim = range(shps$unstructure(shps$y_)))
shps$pole_ <- list(shp$y_)[rep(1, length(id))]

# check functions
par(mfrow=c(3,3), mar = c(0,2,2,0))
shps$plot(t="b", seg_par = list(lty="dotted"))
logy_ <- shps$log(shps$y_, shps$pole_)
explogy_ <- shps$exp(logy_, shps$pole_)
stopifnot(all.equal(explogy_, shps$align(shps$y_, shps$pole_)))
stopifnot(
  all.equal(as.vector(sapply(shps$innerprod(logy_, shps$pole_), Im)), rep(0,8))
)
stopifnot(
  all.equal(logy_, shps$register_v(logy_, shps$pole_))
)
stopifnot(
  all(Mod(sapply(shps$y_, mean))< 1e-15)
)
stopifnot(
  all(abs(crossprod(shps$unstructure(logy_), shps$get_normal())) < 1e-13)
)
# register random tangent vector
set.seed(389)
v_ <- complex(re = rnorm(length(arg)), im = rnorm(length(arg)))
v_r <- shp$register_v(v_, y0_ = shp$y_)
shp$innerprod(v_r, shp$y_)

