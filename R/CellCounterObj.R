#
#
#
#
# CellCounter Object
#
# For storing and processing of individual CellCounter XML files.
#
#





CellCounter <-  R6Class("CellCounterObj",
  public = list(

    data = NA,
    cells = NA,
    metadata = NA ,

    initialize = function(ccfile){
      # '
      # Initialize CellCounterObj by loading CellCounter XML file.
      # '

      #Checks if valid CellCounter Marker File. Fails if true
      ccCheck <- all(XML::xmlToDataFrame(nodes = XML::getNodeSet(XML::xmlParse(ccfile), "//CellCounter_Marker_File")) == 0)
      if(ccCheck){
        return(message(paste0("Invalid CellCounter file:  ", ccfile)))
      }

      self$metadata <- as.list(XML::xmlToDataFrame(nodes = XML::getNodeSet(XML::xmlParse(ccfile), "//Image_Properties")))
      self$data <- XML::xmlToDataFrame(nodes = XML::getNodeSet(XML::xmlParse(ccfile), c("//Name","//Marker")),
                             colClasses = c("character","numeric","numeric","numeric"))
      # Fills NA values down in text(Marker) column to adapt from XML formatting, factors to clear empty levels.
      self$data$text <- factor(self$data$text[
        which(!is.na(self$data$text))[c(1,1:sum(!is.na(self$data$text)))][cumsum(!is.na(self$data$text))+1]])
      # Further data formatting cleanup.
      colnames(self$data) <- c("Marker","X","Y","Z")
      self$data <- self$data[!is.na(self$data$X),]
      rownames(self$data) <- NULL


      # Attempts comprehension of various formats cells may be marked in, with error check.
      cellMarkerNames <- c("nucleus","nuclei","cell","cells","dapi","hoechst")
      cellMarkerIndex <- which(gsub("\\s|_", "", tolower(levels(self$data$Marker))) %in% cellMarkerNames)
      if(length(cellMarkerIndex) != 1){
        return(message(paste0("Unable to parse:  ", as.character(ccfile),
                       "\n   ", as.character(length(cellMarkerIndex)), " base cell types detected!")))
      }
      self$cells <- self$data[which(as.numeric(self$data$Marker)==cellMarkerIndex),-1]

      # Attempts comprehension of various formats that origin/pos0 may be marked in. Adds origin cell indices to metadata
      originMarkerNames <- c("pos0","position0","base","origin","0","origins")
      originMarkerIndex <- which(gsub("\\s|_", "", tolower(levels(self$data$Marker))) %in% originMarkerNames)
      marked <- which(self$data$Marker==levels(self$data$Marker)[originMarkerIndex])
      toMark <- apply(self$data[marked,c("X","Y","Z")],1, function(x){
        order(apply(self$cells[,c("X","Y","Z")], 1, function(y) sum((y - x)**2)))[1] })
      self$metadata$OriginCells <- c(unlist(toMark), use.names=FALSE)

      # Other dummy markers (Type 2,Type3 etc) are being carried over into marker values, so find and store for elimination in next step.
      dummyMarkerIndex <- which(grepl("(Type [0-9]{1})",levels(self$data$Marker)))

      # Generate cells dataframe, fill marker columns with remaining markers. Adds Markers to metadata for later easy reference.
      self$metadata$Markers <- levels(self$data$Marker)[-c(cellMarkerIndex,originMarkerIndex,dummyMarkerIndex)]
      for(marker in self$metadata$Markers){
        self$cells[,marker] <- FALSE
        marked <- which(self$data$Marker==marker)
        # Identifies nearest cell per mark by radius squared.
        toMark <- apply(self$data[marked,c("X","Y","Z")],1, function(x){
          order(apply(self$cells[,c("X","Y","Z")], 1, function(y) sum((y - x)**2)))[1] })
        self$cells[toMark,marker] <- TRUE
      }

      # Process with calibration data if calibration is present.
      if("X_Calibration" %in% names(self$metadata)){
        self$cells$X <- self$cells$X * as.numeric(levels(self$metadata$X_Calibration))
        self$cells$Y <- self$cells$Y * as.numeric(levels(self$metadata$Y_Calibration))
        self$cells$Z <- self$cells$Z * as.numeric(levels(self$metadata$Z_Calibration))
      }

    },

    dist = function(){
      # "Return distance matrix for cells"
      stats::dist(self$cells[,c("X","Y","Z")])
    }
  ),
  active = list(
    Markers = function(v){
      if(missing(v)){ return(self$metadata$Markers) }
      else if(length(v)==length(self$metadata$Markers)){
        colnames(self$cells)[which(colnames(self$cells) %in% self$metadata$Markers)] <- v
        self$metadata$Markers <- v
      }
    }
  )
)
