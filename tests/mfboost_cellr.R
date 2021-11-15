
# load irregular cell data
data("cellr", package = "manifoldboost")

# subsample (one for each covariate combination)
cellsub <- as.data.frame(cellr[-which(names(cellr)=="response")])
cellsub$myd <- factor(with(cellsub,
                           paste0("a=", a, " r=", r, " b=", b, " m=", m)))
subids <- match(unique(cellsub$myd), cellsub$myd)
cellsub <- as.list(cellsub[subids, ])
cellsub$response <- cellr$response
cellsub$response$id <- cellsub$response$id[subids]

# fit model
cell_model <- mfboost(
  formula = response ~ bbsc(a, df = 3, knots = 5) + 
    bbsc(r, df = 3, knots = 5) + 
    bbsc(b, df = 3, knots = 5) + 
    bbsc(m, df = 3, knots = 5),
  obj.formula = value^dim ~ 
    bbs(arg, df = 1, differences = 0, knots = 5, 
        boundary.knots = c(0,1), cyclic = TRUE) | id, 
  data = cellr,
  family = PlanarShapeL2(),
  control = boost_control(mstop = 300)
  )

cross-validation
set.seed(8768)
sheep_cv <- cvrisk(sheep_model, 
                   folds = cvMa(ydim = sheep_model$ydim,
                                type = "kfold"), 
                   grid = 0:mstop(sheep_model))
cell_model[mstop(cell_cv)]

# plot first four predictions
par(mfrow = c(2,2), mar = rep(2, 4) )
plot(cell_model, ids = 1:4, t = "l", 
     main = cells$myd[1:4], 
     seg_par = list(lty = "dashed"))
legend(x = "bottomright", lty = c(1,1, 2),
       legend = c("intercept", "prediction", "point correspondence"), 
       col = c("grey", "black", "grey"))

# compare with data
plot(cell_model, ids = 1:4, t = "l", y0_ = cell_model$family@mf$y_[1:4], 
     main = cellsub$myd[1:4], 
     seg_par = list(lty = "dashed"))
legend(x = "bottomright", lty = c(1,1, 2),
       legend = c("observation", "prediction", "point correspondence"), 
       col = c("grey", "black", "grey"))

# predict dense cells on grids
cellgrid <- cellsub
cellgrid$response <- with(cellgrid$response, expand.grid(
  id = unique(id), 
  arg = seq(0,1, len = 100),
  dim = unique(dim),
  value = NA))
cellgrid$response$value <- predict(cell_model, 
                                   newdata = cellgrid, type = "response")


# factorize effects
cell_fac <- factorize(cell_model)

vimp <- varimp(cell_fac$cov)
plot(vimp, auto.key = FALSE)

# plot two most important effect directions
this <- cell_fac$cov$which(head(names(vimp)[order(vimp, decreasing = TRUE)], 2))
par(mfcol = c(2,2))
plot(cell_fac$resp, which = this, y0_par = list(type="l"))
plot(cell_fac$cov, which = this)


