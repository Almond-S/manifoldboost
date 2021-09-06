# test gradient hypothesis empirically

ngradient <- function(p, v, y) {
  mu <- Exp(v, p)
  epsilon <- Log(y, mu)
  norm_v <- sqrt(innerprod(v))
  norm_epsilon <- sqrt(innerprod(epsilon))
  dvExp <- - sin(norm_v) * tcrossprod(p,v) / norm_v + 
    (norm_v*cos(norm_v) - sin(norm_v)) / norm_v^3 * tcrossprod(v) +
    sin(norm_v) / norm_v * diag(length(p))
  return( norm_epsilon / sin(norm_epsilon) * crossprod(y, dvExp) )
}


ngradient2 <- function(p, v, y) {
  mu <- Exp(v, p)
  epsilon <- Log(y, mu)
  norm_v <- sqrt(innerprod(v))
  norm_epsilon <- sqrt(innerprod(epsilon))
  dvExp <- - sin(norm_v) * tcrossprod(p,v) / norm_v + 
    (norm_v*cos(norm_v) - sin(norm_v)) / norm_v^3 * tcrossprod(v) +
    sin(norm_v) / norm_v * diag(length(p))
  y_ <- y - innerprod(mu, y) * mu
  norm_y_ <- sqrt(innerprod(y_))
  dmuPart <- tcrossprod(y_, y) / norm_y_ -
    innerprod(mu, y) / norm_y_ *
    (diag(length(y)) - tcrossprod(y_) / norm_y_^2 + tcrossprod(mu) / norm_y_^2) *
    norm_epsilon
  return( t(epsilon) %*% dmuPart %*% dvExp )
}

residual <- function(p, v, y) {
  mu <- Exp(v, p)
  epsilon <- Log(y, mu)
  return( Transport(epsilon, mu, p) )
}

# check first example:

o <- c(0,0,0)
p <- c(0,0,1)
v <- c(pi/3,0,0) 
y <- cos(pi/4) * p + sin(pi/4) * c(0,1,0)

mu <- Exp(v,p)
epsilon <- Log(y,mu)
res <- residual(p, v, y)
ngrad <- c(ngradient(p, v, y))
ngrad2 <- c(ngradient2(p, v, y))

Expv <- lapply(seq(0,1,len = 20), function(t) Exp(t*v, p))
Expeps <- lapply(seq(0,1,len = 20), function(t) Exp(t*epsilon, mu))

library(rgl)
rgl::spheres3d(0, col = "lightblue", alpha = .7)
rgl::arrow3d(o, p, col = "black")
rgl::arrow3d(o, mu, col = "darkgreen")
for(pt in Expv) rgl::points3d(pt[1], pt[2], pt[3], col = "darkgreen")
rgl::arrow3d(o, y, col = "darkred")
for(pt in Expeps) rgl::points3d(pt[1], pt[2], pt[3], col = "darkred")
rgl::arrow3d(mu, mu + epsilon, col = "darkred")
rgl::arrow3d(p, p + res, col = "pink")
rgl::arrow3d(p, p + ngrad, col = "grey")
rgl::arrow3d(p, p + ngrad2, col = "cornflowerblue")

rgl::arrow3d(p, p + ngrad - innerprod(p, ngrad) * p, col = "grey")
rgl::arrow3d(p, p + ngrad2 - innerprod(p, ngrad2) * p, col = "cornflowerblue")

# do ngrad and res point in the same direction?
acos( innerprod( res, ngrad) / sqrt(innerprod(res)) / sqrt(innerprod(ngrad)) )

# and the tangent part of ngrad only?
tngrad <- ngrad - innerprod(p,ngrad) * p
acos( innerprod( res, tngrad) / sqrt(innerprod(res)) / sqrt(innerprod(tngrad)) )
# -> almost

# and the negative gradient calculated with the alternative procedure
cbind(ngrad, ngrad2)
# not precisely the same length
acos( innerprod( ngrad, ngrad2) / sqrt(innerprod(ngrad)) / sqrt(innerprod(ngrad2)) )
# but precisely the same direction
