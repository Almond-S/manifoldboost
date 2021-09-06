
library(R6)

Geom <- R6Class("Geom", lock_objects = FALSE,
  public = list(
    initialize = function(x) {
      if(missing(x)) return(invisible(self))
      private$point <- x
      invisible(self)
    },
    translate = function(x) {
      private$point <- private$point + x
    },
    datum_point = NA,
    plot = function() {
      if(is.null(dev.list())) plot.new()
      segments(0, 0, private$point[1], private$point[2])
    }
  ),
  private = list(
    point = NULL,
    hamster = "Nager"
  )
)

Geom2 <- R6Class("Geom2", inherit = Geom,
                 public = list(
                   initialize = function(x, datum_point = NULL) {
                     super$initialize(x)
                     if(!is.null(datum_point)) 
                       self$datum_point <- datum_point
                   }
                 ))

a <- Geom2$new(datum_point = function(x) 5)
a$datum_point(1)
a$initialize(datum_point = c(3,2))
a$.__enclos_env__$private$hamster
