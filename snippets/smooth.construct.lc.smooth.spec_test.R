library(dplyr)
library(Matrix)
library(plotly)
library(plot3D)
library(rgl)
# generate random surface with zero line constraint


t_grid <- seq(0,1, len = 50)

dat0 <- data.frame( t = t_grid, y = 1)

library(mgcv)

mod0 <- gam(y ~ s(t, bs = "ps"), data = dat0, fit = FALSE)

B <- mod0$X
P <- mod0$S[[1]]


# generate symmetric surface and test symm smooth -------------------------

set.seed(4590)
Coefs <- matrix( rnorm(ncol(B)^2), ncol = ncol(B) )
D <- chol(solve(100*P + diag(nrow = ncol(B)-1)))
Coefs[-1,-1] <- D %*% Coefs[-1,-1] %*% D
Coefs <- forceSymmetric(Coefs) %>% as.matrix

persp3d( x = t_grid, y = t_grid, z= tcrossprod( B %*% Coefs, B), 
         col = "cornflowerblue", theta = 30, phi = 50 )

# now fit the symmetric smooth
library(sparseFLMM)

data <- expand.grid(t1 = t_grid, t2 = t_grid)
data$y <- c(tcrossprod( B %*% Coefs, B))

smooth.construct.symm.smooth.spec <- sparseFLMM:::smooth.construct.symm.smooth.spec

mod <- gam(y ~ s(t1, t2, bs = "symm"), data = data)
Pred <- predict(mod) %>% matrix(nrow = nrow(B))

persp3d( x = t_grid, y = t_grid, z= Pred, col = "orange", 
         theta = 30, phi = 50, add = TRUE )

mod2 <- gam(y ~ -1 + s(t1, t2, bs = "symm"), data = data)
Pred2 <- predict(mod2) %>% matrix(nrow = nrow(B))

persp3d( x = t_grid, y = t_grid, z= Pred2, col = "darkred", 
         theta = 30, phi = 50, add = TRUE )

# => when intercept is excluded, there isn't another intercept inlcuded into the p-spline


# check out standard point constraint -------------------------------------

mod_pc <- gam(y ~ -1+ te(t1, t2, bs = "ps", pc = c(1,1)), 
              data = data)
Pred_pc <- predict(mod_pc) %>% matrix(nrow = nrow(B))

persp3d( x = t_grid, y = t_grid, z= Pred, col = "cornflowerblue", alpha = .5,
         theta = 30, phi = 50 )
persp3d( x = t_grid, y = t_grid, z= Pred_pc, col = "orange",
         theta = 30, phi = 50, add = TRUE )



# generate zero line constraint surface and test lc smooth ----------------

lmod <- gam(y ~ -1 + s(t, bs = "ps", pc = 0), data = dat0, fit = FALSE)
lB <- lmod$X
matplot(lB, t = "l")

set.seed(4590)
lCoefs <- matrix( rnorm(ncol(lB) * ncol(B)), nrow = ncol(lB) )

lY <- tcrossprod( lB %*% lCoefs, B)
persp3d( x = t_grid, y = t_grid, z= lY, 
         col = "cornflowerblue", theta = 30, phi = 50 )
data$ly <- c(lY)

# try to fit it
source("R/smooth.construct.lc.smooth.spec.R")
# use pc to specify a line constraint rather than a point constraint
# for each of the dimentsions t1 and t2 setting the repsective spline value to 
# zero for all the whole axis, such that, e.g., c(0,NA) forces
# s(0, t2) = 0 for all t2. The NA for t2 specifies no corresponding constraint for t2.
lmod1 <- gam(ly ~ -1 + s(t1, t2, bs = "lc", pc = c(0,NA)), data = data)
lPred <- predict(lmod1)
plot(lmod1)

matplot(lmod1$smooth[[1]]$margin[[1]]$X[1:50, ], t = "l")
matplot(lmod1$smooth[[1]]$margin[[2]]$X[(0:49) + 50*(0:49), ], t = "l")

persp3d( x = t_grid, y = t_grid, z= lY, col = "cornflowerblue", alpha = .5,
         theta = 30, phi = 50, xlab = "t1", ylab = "t2" )
persp3d( x = t_grid, y = t_grid, z= lPred, col = "orange",
         theta = 30, phi = 50, xlab = "t1", ylab = "t2", add = TRUE )

# => pefect fit!


# now, check out on incomplete data only living on triangle ---------------------------------------

data %>% ggplot(aes(t1, t2, fill = t1 <= t2 & t1 <= 1-t2)) + geom_raster()

data_tri <- filter(data, t1 <= t2 & t1 <= 1-t2)

lmod2 <- gam(ly ~ -1 + s(t1, t2, bs = "lc", pc = c(0,NA)), data = data_tri)
data_tri$lp_tri <- c(predict(lmod2))
data <- data %>% left_join(select(data_tri, t1, t2, lp_tri))

plot(lmod2)

persp3d( x = t_grid, y = t_grid, z= lY, col = "cornflowerblue", alpha = .5,
         theta = 30, phi = 50, xlab = "t1", ylab = "t2" )
persp3d( x = t_grid, y = t_grid, z= data$lp_tri, col = "orange",
         theta = 30, phi = 50, xlab = "t1", ylab = "t2", add = TRUE )

# => NICE !!!
