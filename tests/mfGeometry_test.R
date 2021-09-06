library(manifoldboost)

# check first example:

o <- c(0,0,0)
p <- c(0,0,1)
v <- c(pi/3,0,0) 
y <- cos(pi/4) * p + sin(pi/4) * c(0,1,0)

s <- mfGeomUnitSphere$new()

# define pseudo-residual and negative gradient ----------------------------

ngradient <- function(p, v, y) {
  mu <- s$exp(v, p)
  epsilon <- s$log(y, mu)
  norm_v <- sqrt(s$innerprod(v))
  norm_epsilon <- sqrt(s$innerprod(epsilon))
  dvLog <- - sin(norm_v) * tcrossprod(p,v) / norm_v + 
    (norm_v*cos(norm_v) - sin(norm_v)) / norm_v^3 * tcrossprod(v) +
    sin(norm_v) / norm_v * diag(length(p))
  return( norm_epsilon / sin(norm_epsilon) * crossprod(y, dvLog) )
}

residual <- function(p, v, y) {
  mu <- s$exp(v, p)
  epsilon <- s$log(y, mu)
  return( s$transport(epsilon, mu, p) )
}


# test computation --------------------------------------------------------


mu <- s$exp(v,p)
epsilon <- s$log(y,mu)
res <- residual(p, v, y)
ngrad <- c(ngradient(p, v, y))

Exp_v <- lapply(seq(0,1,len = 20), function(t) s$exp(t*v, p))
Exp_eps <- lapply(seq(0,1,len = 20), function(t) s$exp(t*epsilon, mu))

library(rgl)
rgl::spheres3d(0, col = "lightblue", alpha = .7)
rgl::arrow3d(o, p, col = "black")
rgl::arrow3d(o, mu, col = "darkgreen")
for(pt in Exp_v) rgl::points3d(pt[1], pt[2], pt[3], col = "darkgreen")
rgl::arrow3d(o, y, col = "darkred")
for(pt in Exp_eps) rgl::points3d(pt[1], pt[2], pt[3], col = "darkred")
rgl::arrow3d(mu, mu + epsilon, col = "darkred")
rgl::arrow3d(p, p + res, col = "pink")
rgl::arrow3d(p, p + ngrad, col = "grey")

rgl::arrow3d(p, p + ngrad - s$innerprod(p, ngrad) * p, col = "grey")

# do ngrad and res point in the same direction?
acos( s$innerprod( res, ngrad) / sqrt(s$innerprod(res)) / sqrt(s$innerprod(ngrad)) )

# and the tangent part of ngrad only?
tngrad <- ngrad - s$innerprod(p,ngrad) * p
acos( s$innerprod( res, tngrad) / sqrt(s$innerprod(res)) / sqrt(s$innerprod(tngrad)) )
# -> almost


# test rotation alignment -------------------------------------------------


