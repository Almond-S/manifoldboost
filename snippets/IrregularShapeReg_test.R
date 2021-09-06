library(mboost)
library(FDboost)
library(Matrix)
library(dplyr)
library(ggplot2)
library(manipulate)
library(mgcv)
library(stringr)
source("R/smooth.construct.sps.smooth.spec.R")
source("R/basic_geometry.R")
source("R/shape_long-class.R")
source("R/preshape.R")
source("R/meanshape.R")

load("data/sheep_data.RData")

# only use landmark data
sheep_bones <- as_shape_long(sheep_bones[ str_detect(rownames(sheep_bones), "landmark") , , ])

# use numeric .arg 
sheep_bones <- sheep_bones %>% group_by(.dim, .id) %>% 
  mutate( arg_id = .arg, .arg = 1:n() ) %>% ungroup
class(sheep_bones) <- c("shape_long", "data.frame")

set.seed(9034)
sheep_bones %>% as_shape_wide %>% filter( .id %in% sample(levels(.id), 12) ) %>% 
  ggplot(aes(.value1, .value2, group = .id, col = .arg)) +
  geom_path() +
  facet_wrap( ~ .id) +
  coord_fixed() +
  theme_minimal()

# clean sheep info --------------------------------------------------------

# na class in sex
sheep_info[is.na(sheep_info$Sex), "Sex"] <- "na"

## join with additional informations
sheep_lk <- sheep_bones %>% left_join(sheep_info, by = c(.id = "ID"))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# B-spline mean shapes for regular / irregular / sparse shapes -----------------------------------------------------

# estimate meanshape and its tangent space basis
mean_sheep <- bshape(sheep_bones, smoothed.cov = TRUE, arg.grid = seq(1,max(sheep_bones$.arg), len = 200), 
               return.grid = 200)
mean_sheep$meanshape %>% as_shape_wide %>% ggplot(aes(.value1, .value2, group = .id)) + coord_fixed() + geom_path()

 dat <- dat %>% mutate(.dim = factor(.dim))
class(dat) <- c("shape_long", "data.frame")
# must be orthogonalized with respect to rotation also...
# irrshp <- IrrShapeReg(pole = pole$meanshape)
.dim <- dat$.dim
bbl <-  brandom(.dim, lambda = 0) %X% bbsc(dat$.arg, lambda = 1, knots = 4)  
p <- pole$meanshape %>% as_shape_complex %>% .$.value
cmat2vmat <- function(p) cbind(c(Re(p), Im(p)), c(-Im(p), Re(p)))
p <- do.call(rbind, lapply(split(p, as_shape_complex(dat)$.id), cmat2vmat))

mean_model <- mboost_fit(list(pole$btangent,
                              buser(p)), response = dat$.value,
                      # family = irrshp,
                      control = boost_control(mstop = 500, nu = .5))
dat_mean <- dat
dat_mean$.value <- predict(mean_model, off2int = TRUE) #extract(pole$btangent, "design") %*% coef(mean_model)[[1]]
dat_mean <- dat_mean %>% as_shape_wide #%>% ggplot(aes(.value1, .value2, group = .id)) + geom_path()
dat_mean[, c(".mean1", ".mean2")] <- as_shape_wide(pole$meanshape)[, c(".value1", ".value2")]

ggplot(data = NULL, aes(.value1, .value2, group = .id)) + 
  geom_path(data = as_shape_wide(dat))
dat_mean %>% ggplot(aes(group = .id)) + #geom_path(aes(.mean1, .mean2)) +
  geom_path(aes(.value1, .value2))

