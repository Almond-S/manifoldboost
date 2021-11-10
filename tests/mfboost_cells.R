
data("cells")

# subsample
cdat <- as.data.frame(cells[-which(names(cells)=="response")])
set.seed(9320)
subids <- sample(1:nrow(cdat), size = 100)
cdat <- as.list(cdat[subids, ])
cdat$response <- cells$response[as.numeric(cells$response$id) %in% subids, ]

# fit model
cell_model <- mfboost(
  formula = response ~ bbsc(a, df = 1, differences = 0, knots = 5) + 
    bbsc(r, df = 1, differences = 0, knots = 5) + 
    bbsc(b, df = 1, differences = 0, knots = 5) + 
    bbsc(m, df = 1, differences = 0, knots = 5),
  obj.formula = value^dim ~ 
    bbs(arg, df = 2, differences = 0, knots = 5, 
        boundary.knots = c(0,1), cyclic = TRUE) | id, 
  data = cdat,
  family = PlanarShapeL2()
  )

plot(cell_model)


