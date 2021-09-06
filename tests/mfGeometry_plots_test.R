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
shp <- mfGeomPlanarShape$new(data = dat1, form = form)
shp$plot()

# add pole
shp$pole_ <- shp$register(complex(re = tail(x, 20), im = tail(y, 20))) * 
  exp(complex(i = 1) * pi / 2)
shp$plot(t = "l")
shp$plot(pch = 19, y0_par = list(type = "b", pch = 5), 
         seg_par = list(lty = "dotted"))

# check product geometry --------------------------------------------------

shps <- mfGeomProduct$new(shp, dat_FDboost, form = form)
set.seed(9403)
shps$pole_ <- Map(`*`, sample(shps$y_, length(shps$y_)),  
  exp(complex(i = 1) * runif(length(shps$y_))) )
shps$plot(y_ = shps$y_, shps$pole_)
shps$plot(type = "l", y0_par = list(type = "b")) # -> nice!

