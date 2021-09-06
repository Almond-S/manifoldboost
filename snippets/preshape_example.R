library(manifoldboost)
# load data
dat <- shapes::digit3.dat


# as landmark shape -------------------------------------------------------

plot(dat[,,1], t = "l")

shp <- as_shape(dat)

# single preshape
plot(preshape(dat[,,1]), t = "l")

shp_range <- range(c(shp[,1:2]))
opar <- par(mfrow = c(2, 2))
lapply(split(shp, shp$.id)[1:4], function(x) 
  plot(x[, 1:2], t = "l", ylim = shp_range, xlim = shp_range))
par(opar)

# plot preshapes
preshp <- preshape(shp)
preshp_range <- range(c(preshp[,1:2]))
opar <- par(mfrow = c(2,2))
lapply(split(preshp, shp$.id)[1:4], function(x) 
  plot(x[, 1:2], t = "l", ylim = preshp_range, xlim = preshp_range))
par(opar)

# successively rotate shape
prex <- preshape(dat[,,1])
xylim <- max(abs(c(prex)))
plot(prex, t = "l", 
     xlim = c(-1,1)*xylim, ylim = c(-1,1)*xylim)
theta_grid <- seq(0, .5*pi, len = 5)
lapply(lapply(theta_grid, rotate, x = prex), lines)

# rotation align rotated shapes again to the first
opar <- par(mfrow = c(2,2))
apply(rotate(
  structure( sapply(theta_grid, rotate, x = prex), 
             dim = c(dim(prex), length(theta_grid)) ), x0 = prex)[,,1:4], 
      3, function(x) 
  plot(x, t = "l", ylim = preshp_range, xlim = preshp_range))
par(opar)


# as outline shape parametrized with respect to arc-length ----------------

shp <- as_shape_long(shp)

library(dplyr)
library(ggplot2)
shp <- shp %>% as_shape_complex %>% group_by(.id) %>% 
  mutate(.arg = c(0, cumsum(Mod(diff(.value))))) %>% 
  mutate(.arg = .arg/.arg[n()]) %>% ungroup 
class(shp) <- c("shape_complex", "data.frame")
shp <- shp %>% as_shape_long

preshp <- shp %>% preshape(int.weights = "trapezoidal")
preshp %>% as_shape_wide %>% 
  ggplot(aes(.value1, .value2, group = .id)) +
  geom_path()

preshp <- shp %>% as_shape_wide %>% preshape(int.weights = "trapezoidal")
preshp %>% 
  ggplot(aes(.value1, .value2)) +
  geom_path() + facet_wrap(~.id)

preshp <- shp %>% as_shape_complex %>% preshape(int.weights = "trapezoidal")
preshp %>% as_shape_wide %>%
  ggplot(aes(.value1, .value2)) +
  geom_path() + facet_wrap(~.id)
