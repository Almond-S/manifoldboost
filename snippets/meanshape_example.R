source("R/basic_geometry.R")
source("R/shape_long-class.R")
source("R/preshape.R")
source("R/meanshape.R")

# # check out shapes
# library(manipulate)
# manipulate(plot(t(fdasrvf::beta[,,i,j]), t = "l"), 
#            i = slider(1, dim(fdasrvf::beta)[3]),
#            j = slider(1, dim(fdasrvf::beta)[4]))

# dat <- fdasrvf::beta[,,4,] # heart shapes
# dat <- aperm(dat, c(2,1,3))

dat <- shapes::digit3.dat

# dat <- fdasrvf::beta[,,17,] # cow shapes
# dat <- aperm(dat, c(2,1,3))
# dat[,1,c(8:9, 15)] <- -dat[,1,c(8:9, 15)]

pre <- preshape(dat)
xylim <- c(-1,1) * max(abs(pre))
plot(x = 0, y = 0, xlim = xylim, ylim = xylim)
apply(pre[,,1:7], 3, lines, t = "l", col = "grey")

lines(meanshape(pre))

# B-spline mean shapes for regular / irregular / sparse shapes -----------------------------------------------------

library(mboost)
library(FDboost)
library(Matrix)
library(dplyr)
library(ggplot2)
library(manipulate)
source("R/smooth.construct.sps.smooth.spec.R")

# parametrize with respect to (relative) arc-length
x <- pre %>% as_shape_complex %>% group_by(.id) %>% 
  mutate(arc = c(0, cumsum(Mod(diff(.value))))) %>% 
  mutate(arc = arc/arc[n()]) %>% ungroup 
class(x) <- c("shape_complex", "data.frame")
arc <- x$arc
x <- x %>% as_shape_long
x$arc <- arc


manipulate( {
  lambdas <- .1^c(Inf, 6:0)
  x_ <- x
  if(parametrize.arclen) x_$.arg <- arc
  if(sparsify.ratio) {
    set.seed(8473)
    argid <- paste(x_$.arg, x_$.id)
    subsam <- sample(argid, ceiling((1 - sparsify.ratio) * length(argid)))
    x_ <- x_[argid %in% subsam, ]
  }
  shapemean_dense <- bshape(x_, smoothed.cov = smoothed.cov,
                            lambda0 = c(lambdas[lambda0], lambdas[lambda2]), #c(0, 1e-2) #1e-2, 
                            int.weights = if(int.weights != "null") int.weights, #"simpson", #"trapezoidal"
                            cv.seed = 90348, cv.k = if(CV != "no CV") 10,
                            arg.grid = if(do.arg.grid) seq(min(x_$.arg), max(x_$.arg), len = 200), # NULL
                            return.grid = 200)$meanshape
  ggplot(NULL) + 
    pre %>% as_shape_wide %>%
    geom_path(mapping = aes(.value1, .value2, group = .id), col = "grey") +
    shapemean_dense %>% as_shape_wide %>%
    geom_path(mapping = aes(.value1, .value2)) +
    coord_fixed() +
    theme_minimal()
  }, 
  lambda0 = slider(1,8, 2), 
  lambda2 = slider(1,8, 2),
  int.weights = picker("empirical", "simpson", "trapezoidal", "null"),
  do.arg.grid = checkbox(FALSE),
  parametrize.arclen = checkbox(FALSE),
  sparsify.ratio = slider(0, 1, step = .05),
  CV = picker("no CV", "10-fold"),
  smoothed.cov = checkbox(FALSE)
)

# illustrate sparsening
manipulate( {
  x_ <- x
  if(parametrize.arclen) x_$.arg <- arc
  if(sparsify.ratio > 0) {
    set.seed(8473)
    argid <- paste(x_$.arg, x_$.id)
    subsam <- sample(argid, ceiling((1 - sparsify.ratio) * length(argid)))
    x_ <- x_[argid %in% subsam, ]
    x_ <- preshape(x_, int.weights = if(int.weights != "null") int.weights)
  }
  ggplot(NULL) + 
    pre %>% as_shape_wide %>%
    geom_path(mapping = aes(.value1, .value2, group = .id), col = "cornflowerblue") +
    pre %>% as_shape_wide %>%
    geom_point(mapping = aes(.value1, .value2, group = .id), col = "cornflowerblue") +
    x_ %>% as_shape_wide %>%
    geom_path(mapping = aes(.value1, .value2, group = .id), col = "black") +
    x_ %>% as_shape_wide %>%
    geom_point(mapping = aes(.value1, .value2, group = .id), col = "black") +
    coord_fixed() +
    facet_wrap(~.id) +
    theme_minimal()
},
  sparsify.ratio = slider(0, 1, step = .05), 
  int.weights = picker("empirical", "simpson", "trapezoidal", "null"),
  parametrize.arclen = checkbox(FALSE)
)

