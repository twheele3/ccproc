#
#
#
#
#
#
#
#
#
#
# Class: Crypt
#
# Description: A class to hold vector information for a crypt respective to cells from a CellCounter object.
#
#
library(R6)
library(ggplot2)

Crypt <-  R6Class("CryptObj",
  public = list(

    origin = NA,
    # Index of cells respective to their row in cc$cells of CellCounter object from Image
    cells = NA,
    # Image that cells originated from
    Image = NA,
    # Position of a cell respective to origin (scalar integer)
    position = NA,
    # Index of edges between cells in self$cells, such that edge x is edge formed between cells[x:x+1]
    E.index = NA,

    initialize = function(origin,cells,Image,E.index){
      self$origin <- origin
      self$cells <- cells
      self$Image <- Image
      self$position <- (1:length(self$cells)) - which(self$cells == self$origin)
      self$E.index <- E.index
    },

    pos = function(x){
      # Returns position of cell x with respect to origin cell.
      return(sapply(x,function(y){which(self$cells==y)-which(self$cells==self$origin)}))
    },

    pos.abs = function(x){
      # Returns position of cell x upward from origin cell.
      return(abs(self$pos(x)))
    },

    in.crypt = function(x){
      # Returns boolean for whether index/indices are in crypt object.
      return(x %in% self$cells)
    },

    neighbors = function(x,dist=1){
      # Returns neighbor cell indices at a distance d to cell x.
      xIndex <- which(self$cells==x)
      n <- c(-dist,dist) + xIndex
      return(n[n >= 1 & n < length(self$cells)])
    }
  )
)

