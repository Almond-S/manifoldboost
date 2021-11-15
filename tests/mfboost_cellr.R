

# modeling the SHAPE of REGULAR curves -------------------

# load regular cell data
data("cellr", package = "manifoldboost")

# subsample (one for each covariate combination)
cellsub <- as.data.frame(cellr[-which(names(cellr)=="response")])
cellsub$myd <- factor(with(cellsub,
                           paste0("a=", a, " r=", r, " b=", b, " m=", m)))
subids <- match(unique(cellsub$myd), cellsub$myd)
cellsub <- as.list(cellsub[subids, ])
cellsub$response <- cellr$response
cellsub$response$dims$id <- ordered(cellsub$response$dims$id[subids], 
                                    levels = unique(cellsub$response$dims$id[subids]))
cellsub$response$mets$value <- cellsub$response$mets$value[,,subids]
class(cellsub$response) <- "tbl_cube" 

# fit SHAPE model
cell_model <- mfboost(
  formula = response ~ bbsc(a, df = 3, knots = 5) + 
    bbsc(r, df = 3, knots = 5) + 
    bbsc(b, df = 3, knots = 5) + 
    bbsc(m, df = 3, knots = 5),
  obj.formula = value^dim ~ 
    bbs(arg, df = 1, differences = 0, knots = 5, 
        boundary.knots = c(0,70), cyclic = TRUE) | id, 
  data = cellsub,
  family = PlanarShapeL2(),
  control = boost_control(mstop = 100)
  )

# # cross-validation
# set.seed(8768)
# cell_cv <- cvrisk(cell_model,
#                    folds = cvMa(ydim = cell_model$ydim,
#                                 type = "kfold"),
#                    grid = 0:mstop(cell_model))
# cell_model[mstop(cell_cv)]

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
cellgrid$response <- cubelyr::tbl_cube(
  dimensions = list(
    id = cellsub$response$dims$id, 
    arg = seq(0,70, len = 100),
    dim = unique(cellsub$response$dims$dim)
  ),
  measures = list(
    value = array(NA, dim = c(29, 100, 2))))

cellgrid$response$mets$value <- array(predict(cell_model, 
                                   newdata = cellgrid, type = "response"), 
                                   dim = c(29,100,2))

for(i in 1:4) 
  plot(cellgrid$response$mets$value[i,,], t = "l")

# factorize effects
cell_fac <- factorize(cell_model)

vimp <- varimp(cell_fac$cov)
plot(vimp, auto.key = FALSE)

# plot two most important effect directions
this <- cell_fac$cov$which(head(names(vimp)[order(vimp, decreasing = TRUE)], 2))
par(mfcol = c(2,2))
plot(cell_fac$resp, which = this, y0_par = list(type="l"))
plot(cell_fac$cov, which = this)


