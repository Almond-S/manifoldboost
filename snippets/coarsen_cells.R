library(dplyr)
library(ggplot2)
library(magrittr)

source("data_coarsening.R")

load("cell.RData")

# format cell data as needed for coarsen()
cell <- as.data.table(cell)
cell[, time := arg]

cell <- split(cell, cell$contour)

blacklists <- lapply(cell, coarsen, do.blacklist = TRUE, 
                    tol = .1, do.plot = FALSE, 
                    verbose = FALSE,  method = "polygon")

blacklists %<>% bind_rows(.id = "contour")
blacklists %<>% group_by(contour) %>% 
  mutate( coarsen_SSE = squ_error,  
          coarsen_relSSE = coarsen_SSE/total_var )

# merging blacklist and cell data
cell <- bind_rows(cell) 
cell <- blacklists %>% 
  select(contour, time, coarsen_SSE, coarsen_relSSE) %>% 
  right_join(cell, by = c("contour", "time"))

cell[is.na(cell$coarsen_SSE), ] %<>% 
  mutate( coarsen_SSE = Inf, coarsen_relSSE = Inf )

cell %<>% ungroup %>% select(-time, -score)

save(cell, file = "cell.RData")
