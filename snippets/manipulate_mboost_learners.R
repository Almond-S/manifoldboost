
library(mboost)

dat <- USArrests[order(USArrests$UrbanPop), ]
dat$Rape_discrete <- findInterval(dat$Rape, quantile(dat$Rape, seq(0,1,len = 5)))
dat$Rape_discrete <- factor(dat$Rape_discrete)

mod <- mboost(Murder ~ bols(Assault) + bols(Rape_discrete) + bbs(UrbanPop), data = dat)
mod2 <- eval(mod$call)

plot(mod)

# add dublicate of UrbanPop learner ---------------------------------------

# extract central environment
e <- environment(mod$coef)
w <- e$thiswhich("UrbanPop")

# dublicate basemodels = bl
mod$basemodel$`bbs(UrbanPop)-2` <- mod$basemodel$`bbs(UrbanPop)`
e$bl$`bbs(UrbanPop)-2` <- e$bl$`bbs(UrbanPop)`
# dublicate baselearners = blg
mod$baselearner$`bbs(UrbanPop)-2` <- mod$baselearner$`bbs(UrbanPop)`
e$blg$`bbs(UrbanPop)-2` <- e$blg$`bbs(UrbanPop)`
# add fit containers
nselect <- sum(e$xselect==w)
e$xselect <- c(e$xselect, rep(4, nselect))
e$ens <- c(e$ens, e$ens[e$xselect == w])
# change design of 
# extend bnames list
e$bnames <- names(e$bl)

plot(mod[length(e$ens)])

# implant other design matrix and coefs -----------------------------------

# orhtogonalize UrbanPop Design Matrix
X <- extract(mod, "design", which = "UrbanPop")[[1]]
cf <- coef(mod, which = "UrbanPop")[[1]]

QR <- qr(X)
Q <- qr.Q(QR)
R <- qr.R(QR)[, order(QR$pivot)]

# implant design
d <- environment(mod$baselearner[[w]]$dpp)
d$X <- Q

# Tadaaa! Both broken :-)
opar <- layout(t(1:2))  
for(i in 0:1) {
  p <- predict(mod, which = w+i)
  plot(dat$UrbanPop, p, t = "b") 
}
par(opar)

# implant coefs in double
for(i in which(e$xselect == 4)) {
  e$ens[[i]]$model <- R %*% e$ens[[i]]$model
}


plot(mod2, which = w)
lines(dat$UrbanPop, 
       predict(mod, which = w), pch = 3, col = "darkred")
lines(dat$UrbanPop,
      predict(mod, which = 4), pch = 3, col = "cornflowerblue")
plot(mod)
# nice, everything under control 3-D
