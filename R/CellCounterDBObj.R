#
#
#
#
#
#
# CellCounterDatabase
#
# An object for storing and processing CellCounter files
#
#
#
library(R6)

CellCounterDB <-  R6Class('CellCounterDatabaseObj',
  public = list(
    cells = NA,
    CCFiles = NA,
    crypts = NA,
    metadata = NA,
    Markers = NA,
    # MarkerColors = NA,

    initialize = function(filelist){
      self$CCFiles <- list()
      self$crypts <- list()
      self$metadata <- list()
      self$Markers <- c()
      for(file in filelist){
        self$addCC(file)
      }
    },

    addCC = function(ccfile,replace=FALSE){
      # '
      # Adds new CellCounterObj to list, then adds cells to database data.frame.
      # Need method to run a list of files through, potential integration with other methods for post-processing.
      # '
      cc <- ccproc::CellCounter$new(ccfile)
      # private$create_triangulation(cc)
      image_name <- as.character(cc$metadata$Image_Filename)
      # Remove if replacing
      if( replace & (sum(image_name %in% names(self$CCFiles)) > 0) ){
        self$removeCC(image_name)
      }
      # Check that image is not already loaded.
      if(sum(image_name %in% names(self$CCFiles)) == 0){
        tryCatch({
          self$CCFiles[[image_name]] <- cc
          newcells <- cbind(cc$cells,Image=rep(image_name,nrow(cc$cells)),Index=1:nrow(cc$cells))
          if("Attributes" %in% names(self$metadata)){
            for(attribute in names(self$metadata)){
              tryCatch({
              newcells[,attribute] <- factor(
                private$string_extract(
                  as.character(newcells[,self$metadata$Attributes[[attribute]]$source_col]),
                  self$metadata$Attributes[[attribute]]$Regex))},
              error = function(err){message(paste("Error in tagging attributes for ",image_name,": ",err))},
              warning = function(w){message(paste("Warning in tagging attributes for ",image_name,": ",w))})
            }
            newcells$Sample <-factor(private$string_extract(as.character(newcells$Image),self$metadata$SampleRegex)[1])
          }
          if(all(is.na(self$cells))){
            self$cells <- newcells
          }else{
            self$cells <- private$rbind_fill(self$cells, newcells)
          }
          # Add new markers to Markers
          self$Markers <- sort(unique(c(self$Markers,cc$Markers)))
          # Reorder cells columns with Markers last.
          self$cells <- self$cells[,c(colnames(self$cells)[!(colnames(self$cells) %in% self$Markers)],self$Markers)]
        },
        error = function(err){message(paste("Error in loading ",image_name,": ",err))},
        warning = function(w){message(paste("Warning in loading ",image_name,": ",w))})
      }else{message(paste0(image_name, " not loaded, already found in database."))}
    },

    removeCC = function(Image){
      # '
      # Remove CellCounterObject from data set. Requires image name(s).
      # '
      to_remove <- names(self$CCFiles) %in% Image
      if(any(to_remove)){
        self$CCFiles[[which(to_remove)]] <- NULL
        self$cells <- self$cells[-which(self$cells$Image %in% Image),]
      }
    },

    save = function(file=""){
      # '
      # Save CellCounterDB to file in R external object format.
      # '
      base::saveRDS(self, file = file)
    },

    processFeatures = function(ccfiles="all",reprocess=FALSE,feature="Crypt"){
      # '
      # Process loaded cells into features. Expand in future to non-crypt stuff.
      # '

      # Add feature column if not yet present.
      if(!(feature %in% colnames(self$cells))){
        self$cells[,feature] <- 0
        self$cells[,paste0("Position.",feature)] <- NA
      }
      self$cells[is.na(self$cells[,feature]),feature] <- 0
      # TODO: add name comprehension and reprocessing checks.
      if(ccfiles=="all"){
        ccfiles = names(self$CCFiles)
      }
      # Check that all image names are valid.
      not_in_images <- ccfiles[!(ccfiles %in% names(self$CCFiles))]
      if(length(not_in_images) > 0 ){
        message(paste("Warning: Images not found in database:",not_in_images))
        ccfiles = ccfiles[(ccfiles %in% names(self$CCFiles))]
      }
      # TODO: Put this in loop below properly
      if(reprocess){
        to_null <- unique(self$cells[which(self$cells$Image %in% ccfiles),feature])
        to_null <- to_null[to_null > 0]
        for(i in to_null){
          self[[feature]][[i]]$Image <- NULL
        }
      }
      # If not reprocessing, removes any names from the list that are processed. Checks this based on if all cells are assigned to crypts.
      if(!reprocess){
        exclude <- c()
        for(img in ccfiles){
          if(sum(img==self$cells$Image) == sum(self$cells[img==self$cells$Image,feature]>0)){
            exclude <- c(exclude,img)
          }
        }
        ccfiles <- ccfiles[!(ccfiles %in% exclude)]
      }

      for(img in ccfiles){
        tryCatch({
          # Create triangulation for CCFiles that don't have it yet.
          if(!any("triangulation" %in% names(self$CCFiles[[img]]$metadata))){
            private$create_triangulation(self$CCFiles[[img]])
          }
          # Recreate triangulation for CCFiles that have outdated triangulations, based on rowcounts.
          if(nrow(self$CCFiles[[img]]$metadata$triangulation$P) != nrow(self$CCFiles[[img]]$cells)){
            private$create_triangulation(self$CCFiles[[img]])
          }
          cc <- self$CCFiles[[img]]

          tri <- cc$metadata$triangulation
          # Creates a rolling count of edges for each vector, sorted by distance increasing.
          tri$E.dist.roll <- private$rollcount(tri$E.dist,tri$E)
          # Creates a rolling count of edges for each vector, sorted by beta increasing.
          tri$E.beta.roll <- private$rollcount(tri$E.beta,tri$E)
          # First define a subset of edges (E.subset) based on paring edges where both have more than 3 connections as sorted by distance.
          E.subset <- private$pare_edges_by_rolling(tri$E.beta.roll,connections.max=2)
          # Further refine E.subset by paring edges as sorted by beta.
          E.subset <- E.subset[private$pare_edges_by_rolling(tri$E.dist.roll[E.subset,],connections.max=2)]
          # Removes additional cyclical vectors by determining which edges have vertices both with too many connections.
          E.subset <- E.subset[private$pare_edges_by_count(tri$E[E.subset,],connections.max=2)]
          # Create an aggregate statistic from normalized distance (shifted to 1 for multiplier effect) and beta.
          dist.mean = mean(tri$E.dist[E.subset])
          dist.sd = sd(tri$E.dist[E.subset])
          dist.normalized <- (tri$E.dist[E.subset]-dist.mean)/dist.sd + 1
          # These stats can theoretically be weighted in the future if needed.
          hybrid_stat <- dist.normalized + tri$E.beta[E.subset]

          E.subset <- E.subset[private$pare_edges_by_stat(stat=hybrid_stat, edgeset=tri$E[E.subset,], decreasing=TRUE, connections.max=2)]

          # Sort remaining edges into ordered vectors.
          ordered_sets <- private$pared_edgeset_to_ordered_vectors(tri$E[E.subset,])

          # Break ordered sets into separate crypts by origins.
          crypts <- list()
          unassigned_sets <- list()
          # Cleanup of ordered sets into crypts
          for(i in 1:length(ordered_sets)){
            set <- ordered_sets[[i]]
            origins <- cc$metadata$OriginCells[cc$metadata$OriginCells %in% set]
            looping <- private$check_looping(set) & (length(origins) > 0)
            # Detect for looping, reorient to center on first origin.
            if(looping){
              set <- private$cycle_vector(set,origins[1])
            }
            cleavepoints <- c()
            ind_or <- sort(which(set %in% origins))
            if(length(ind_or)>2){
              # Looks for cleave points between origins.
              for(j in 1:(length(ind_or)-1)){
                subset <- set[ind_or[j]:ind_or[j+1]]
                cut_point <- private$find_cut_point(subset,img,dist.mean,dist.sd)
                cleavepoints <- c(cleavepoints,
                                  which(set==subset[cut_point])[1],
                                  which(set==subset[cut_point])[1]+1)
              }
              # Break loop around first cleavepoint and adjust set vector
              set <- set[-length(set)]
              set <- private$cycle_vector(set,set[cleavepoints[2]])
              cleavepoints <- cleavepoints[-(1:2)]-cleavepoints[1]
              # Set up paired indices of start and endpoints for each set in a 2-column matrix
              cleavepoints <- matrix(c(1,cleavepoints,length(set)),ncol=2,byrow=TRUE)
              for(j in 1:nrow(cleavepoints)){
                subset <- set[cleavepoints[j,1]:cleavepoints[j,2]]
                if(any(origins %in% subset)){
                  crypts <- private$list_append(crypts,subset)
                }else{
                  unassigned_sets <- private$list_append(unassigned_sets,subset)
                }
              }
              # If just a single crypt, processing at this step should be done.
            }else if(length(ind_or)==2){
              crypts <- private$list_append(crypts,set)
              # Put left-overs into unassigned.
            }else{
              unassigned_sets <- private$list_append(unassigned_sets,set)
            }
          }

          # Make crypts loop properly before assigning unassigned.
          for(i in 1:length(crypts)){
            crypts[[i]] <- private$cut_open_crypt(crypts[[i]],img,dist.mean,dist.sd)
          }

          # Assign unassigned sets by relinking with separated crypts.
          # Still experimental, needs testing on real occurrences.
          if(length(unassigned_sets)>0){
            message(paste0("Experimental error-correction algorithm used for ",img,", please review results with 'plot_debug'!"))
            print("Crypts:")
            for(i in crypts){print(i)}
            print("Unassigned:")
            for(i in unassigned_sets){print(i)}
            while(!all(unlist(lapply(unassigned_sets,function(x){is.null(x)})))){

              # First, pull all the endpoints from established crypts.
              endpoints <- c()
              for(i in crypts){
                endpoints <- c(endpoints,i[1],i[length(i)])
              }

              # Grab set. Looping by nullifying as going through.
              to_grab <- which(!unlist(sapply(unassigned_sets,function(x){is.null(x)})))[1]
              set <- unassigned_sets[[to_grab]]
              if(set[length(set)] == set[1]){
                set <- set[1:(length(set)-1)]
              }
              unassigned_sets[[to_grab]] <- NULL
              # Should extract candidate vertices from set
              testedges <- matrix(c(rep(set,each=length(endpoints)),rep(endpoints,times=length(set))),ncol=2)
              edges <- private$areEdges(testedges,tri$E)
              edges <- edges[edges>0]
              # Picks a start edge based on minimum distance * theta
              start.edge <- tri$E[edges[order(tri$E.beta[edges] * tri$E.dist[edges])[1]],]
              # Pull index from unassigned for start edge
              set.start.index <- which(set %in% start.edge[start.edge %in% set])
              # End point should be adjacent in ordered vector, with looping assumed.
              set.end.candidates <- set[c((length(set) + set.start.index + 1) %% length(set),(length(set) + set.start.index - 1) %% length(set))]
              # From candidates, determine which edge is the most likely fit.
              edge.candidates <- tri$E[edges,][private$which_vert(set.end.candidates,tri$E[edges,]),]
              edge.candidate.index <- private$areEdges(edge.candidates,tri$E)
              end.edge <- tri$E[edge.candidate.index[order(tri$E.beta[edge.candidate.index] * tri$E.dist[edge.candidate.index])[1]],]
              set.end.index <- which(set %in% end.edge[end.edge %in% set])
              # Cycle and orient vector to merge with others.
              set.start <- set[set.start.index]
              set.end <- set[set.end.index]
              set <- private$cycle_vector(set,set.start)
              if(set.end != set[length(set)]){
                set <- private$cycle_vector(set[length(set):1],set.start)
              }

              # Figure out which crypt sets to pull.
              crypt.start <- ceiling((1:length(endpoints))/2)[which(endpoints %in% start.edge[start.edge %in% endpoints])]
              crypt.end <- ceiling((1:length(endpoints))/2)[which(endpoints %in% end.edge[end.edge %in% endpoints])]

              if(crypt.start == crypt.end){
                cryptset <- crypts[[crypt.start]]
                # Reorient cryptset vector
                if(which(cryptset %in% start.edge) != length(cryptset)){
                  cryptset <- cryptset[length(cryptset):1]
                }
                crypts[[crypt.start]] <- private$cut_open_crypt(c(cryptset,set),img,dist.mean,dist.sd)
              }else{
                # Process as per prior splitting method. Expecting two crypts so fixed output.
                crypt1 <- crypts[[crypt.start]]
                crypt2 <- crypts[[crypt.end]]
                # Reorient crypt vectors
                if(which(crypt1 %in% start.edge) != length(crypt1)){
                  crypt1 <- crypt1[length(crypt1):1]
                }
                if(which(crypt2 %in% end.edge) != 1){
                  crypt2 <- crypt2[length(crypt2):1]
                }
                set <- c(crypt1,set,crypt2,crypt1[1])

                # Copy paste from previous cleaving algorithm. TODO: Try to get into a separate function to modularize.
                origins <- cc$metadata$OriginCells[cc$metadata$OriginCells %in% set]
                set <- private$cycle_vector(set,origins[1])
                cleavepoints <- c()
                ind_or <- sort(which(set %in% origins))
                # Looks for cleave points between origins.
                for(j in 1:(length(ind_or)-1)){
                  subset <- set[ind_or[j]:ind_or[j+1]]
                  cut_point <- private$find_cut_point(subset,img,dist.mean,dist.sd)
                  cleavepoints <- c(cleavepoints,
                                    which(set==subset[cut_point])[1],
                                    which(set==subset[cut_point])[1]+1)
                }
                # Break loop around first cleavepoint and adjust set vector
                set <- set[-length(set)]
                set <- private$cycle_vector(set,set[cleavepoints[2]])
                cleavepoints <- cleavepoints[-(1:2)]-cleavepoints[1]
                # Set up paired indices of start and endpoints for each set in a 2-column matrix
                cleavepoints <- matrix(c(1,cleavepoints,length(set)),ncol=2,byrow=TRUE)
                set[cleavepoints[j,1]:cleavepoints[j,2]]
                crypts[[crypt.start]] <- private$list_append(crypts,set[cleavepoints[1,1]:cleavepoints[1,2]])
                crypts[[crypt.end]] <- private$list_append(crypts,set[cleavepoints[2,1]:cleavepoints[2,2]])
              }
            }
          }

          # Final processing of crypt vectors and sorting of crypts into CryptObjs and addition to cell database.
          for(v in crypts){
            i <- length(self$crypts) + 1
            i0 <- cc$metadata$OriginCells[cc$metadata$OriginCells %in% v]
            c <- ccproc::Crypt$new(origin=i0,cells=v,Image=cc$metadata$Image_Filename,E.index=private$areEdges(private$vector_to_edges(v),tri$E))
            self$crypts[[i]] <- c
            indices <- self$which_cells(list("Image" = c$Image, "Index" = v))
            self$cells[indices,feature] <- i
            self$cells[indices,paste0("Position.",feature)] <- c$pos(self$cells$Index[indices])
          }
        },
      error = function(err){message(paste("Error: Unable to process ",img,"\n",err,"\n"))},
      warning = function(w){message(paste("Warning in",img,":",w))})
      }
    },

    plot_debug = function(images,savedir=NA){
      for(image in images){
        # Load cells of interest
        cc <- self$cells[self$cells$Image==image,]
        # Reorder based on file index.
        cc <- cc[order(cc$Index),]
        # Rescale to pixel values.
        cc[,"X"] <- cc[,"X"] / as.numeric(levels(self$CCFiles[[image]]$metadata$X_Calibration)[1])
        cc[,"Y"] <- cc[,"Y"] / as.numeric(levels(self$CCFiles[[image]]$metadata$Y_Calibration)[1])
        # Flip Y axis
        cc[,"Y"] <- abs(cc[,"Y"] - max(cc[,"Y"]))
        # Create edges for all cells in crypts.
        edges <- c()
        for(i in unique(cc$Crypt[!is.na(cc$Crypt) & (cc$Crypt > 0)])){
          edges <- c(edges,t(private$vector_to_edges(self$crypts[[i]]$cells)))
        }
        if(length(edges) > 0){
          edges <- matrix(edges,ncol=2,byrow=TRUE)
        }
        # Save to file if save directory specified.
        if(is.character(savedir)){

          tryCatch({
            dev.new()
            par(mar=c(0,0,1,0))
            plot( x = c( min(cc[,"X"]) - 50, max(cc[,"X"]) + 50), y = c( min(cc[,"Y"]) - 50, max(cc[,"Y"]) + 50), col="white",
                  axes = FALSE, main = image, bg = "white") +
              segments( x0 = cc[edges[,1],"X"], y0 = cc[edges[,1],"Y"], x1 = cc[edges[,2],"X"], y1 = cc[edges[,2],"Y"] ) +
              points( x = cc[,"X"], y = cc[,"Y"], col=factor(cc$Crypt), cex=3, lwd=1, pch=21, bg="white") +
              text( x = cc[,"X"], y = cc[,"Y"], labels = as.character(abs(cc$Position.Crypt)), pos = 3, offset = 0.06, cex = 0.75) +
              text( x = cc[,"X"], y = cc[,"Y"], labels = as.character(abs(cc$Index)), pos = 1, offset = 0.09, cex = 0.75)

            dev.copy(png, paste0(savedir, image," debug plot.png"),
                     width = (max(cc[,"X"]) + 100 - min(cc[,"X"])), height = (max(cc[,"Y"]) + 220 - min(cc[,"Y"])), res=120, type="cairo")
            dev.off()
            graphics.off()
          },
          error = function(err){message(paste0("Could not plot ",image,"/n",err,"/n"))})
        }else{
          tryCatch({

            plot( x = c( min(cc[,"X"]) - 50, max(cc[,"X"]) + 50), y = c( min(cc[,"Y"]) - 50, max(cc[,"Y"]) + 50), col="white",
                  axes = FALSE, main = image, bg = "white") +
              segments( x0 = cc[edges[,1],"X"], y0 = cc[edges[,1],"Y"], x1 = cc[edges[,2],"X"], y1 = cc[edges[,2],"Y"] ) +
              points( x = cc[,"X"], y = cc[,"Y"], col=factor(cc$Crypt), cex=3, lwd=1, pch=21, bg="white") +
              text( x = cc[,"X"], y = cc[,"Y"], labels = as.character(abs(cc$Position.Crypt)), pos = 3, offset = 0.06, cex = 0.75) +
              text( x = cc[,"X"], y = cc[,"Y"], labels = as.character(abs(cc$Index)), pos = 1, offset = 0.09, cex = 0.75)
          },
          error = function(err){message(paste0("Could not plot ",image,"/n",err,"/n"))})
        }
      }
    },

    remove_duplicates = function(r=5){
      # TODO: Test this function.
      # Removes duplicate cells within CCFiles on the assumption that any clicks within the test radius are duplicates.
      for(img in names(self$CCFiles)){
        # Radius given in pixels and converted to real distance based on calibration.
        radius = r * as.numeric(levels(self$CCFiles[[img]]$metadata$X_Calibration)[1])

        # Calculate pairs that are closer than stated radius.
        indices <- which(as.matrix(stats::dist(self$CCFiles[[img]]$cells[,c("X","Y","Z")])) <= radius, arr.ind=TRUE)
        if(length(indices) > 0){
          indices = indices[indices[,1] > indices[,2],]
        }

        i <- length(indices)/2
        # Proceed if any are close enough.
        # Iterative process to remove clusters then check whether there are any still too close.
        while(length(indices)>0){
          merged <- t(data.frame(lapply(1:(length(indices)/2), function(x){
            c(colMeans(self$CCFiles[[img]]$cells[c(indices[c(x,x+(length(indices)/2))]),!(colnames(self$CCFiles[[img]]$cells) %in% self$CCFiles[[img]]$metadata$Markers)]),
              colSums( self$CCFiles[[img]]$cells[c(indices[c(x,x+(length(indices)/2))]),self$CCFiles[[img]]$metadata$Markers]) > 0)
          })))

          self$CCFiles[[img]]$cells = rbind(self$CCFiles[[img]]$cells[-unique(c(indices)),],merged)
          rownames(cc) <- NULL
          indices = which(as.matrix(stats::dist(self$CCFiles[[img]]$cells[,c("X","Y","Z")])) <= radius, arr.ind=TRUE)
          if(length(indices)>0){
            indices = indices[indices[,1] > indices[,2],]
          }
          i = i + length(indices)/2
        }
        if(i > 0){
          # Set markers back to boolean
          self$CCFiles[[img]]$cells[,self$CCFiles[[img]]$metadata$Markers] <- self$CCFiles[[img]]$cells[,self$CCFiles[[img]]$metadata$Markers] > 0
          # Remove old cells from database.
          self$cells = self$cells[!(self$cells$Image==img),]
          # Insert new cells.
          self$cells <- private$rbind_fill(self$cells, self$CCFiles[[img]]$cells)
          # Alert that some duplicates detected.
          message(paste(i,"duplicate cells merged in",img))
        }
      }
    },

    which_cells = function(listconditions){
      # Returns indices in self$cells that match conditions specified in list.
      return(which(apply(sapply(names(listconditions),function(x){self$cells[,x] %in% listconditions[[x]] }),1,prod)>0))
    },

    tag_attribute = function(Regex="",attribute="", source_col="Image"){
      # Inserts a column into dataframe as attribute, labeling cells by regex pattern as applied to source_col string.
      if(Regex==""){
        return(message("Error: Attribute tagging requires a RegEx pattern to apply to image names."))
      }
      if(attribute=="" | any(attribute %in% colnames(self$cells))){
        return(message("Error: Attribute must be named and not already in use in cells table."))
      }
      if(!any(colnames(self$cells) %in% source_col)){
        return(message("Error: Source column must be in $cells dataframe."))
      }
      # self$metadata[[paste0(attribute,"Regex")]] <- Regex
      if(!("Attributes" %in% names(self$metadata))){
        self$metadata[["Attributes"]] <- list()
      }
      self$metadata[["Attributes"]][[attribute]] <- list("Regex" = Regex,
                                            "source_col" = source_col)
      self$cells[,attribute] <- factor(
        private$string_extract(
          as.character(self$cells[,self$metadata$Attributes[[attribute]]$source_col]),
          self$metadata$Attributes[[attribute]]$Regex))
    },

    traceDistAlongFeature = function(feature="Crypt",origin=0){
      # Calculate total distance along cells in a feature from origin position by tracing along edges.
      tracecol <- paste0("Dist.From.",origin,".",feature)
      self$cells[,tracecol] <- NA
      for(i in na.omit( unique(self$cells$Crypt))){
        crypt <- self$crypts[[i]]
        if(!is.null(crypt$Image)){
          img <- as.character(self$crypts[[i]]$Image)
          leftward <- ave(
            self$CCFiles[[img]]$metadata$triangulation$E.dist[
              private$areEdges(
                private$vector_to_edges(
                  crypt$cells[
                    which(crypt$cells == crypt$origin): 1]),
                self$CCFiles[[img]]$metadata$triangulation$E)],
            FUN = cumsum)
          rightward <- ave(
            self$CCFiles[[img]]$metadata$triangulation$E.dist[
              private$areEdges(
                private$vector_to_edges(
                  crypt$cells[
                    which(crypt$cells == crypt$origin):length(crypt$cells)]),
                self$CCFiles[[img]]$metadata$triangulation$E)],
            FUN = cumsum)
          dists <- c(leftward[length(leftward):1],0,rightward)
          indices <- self$which_cells(list("Image" = img, "Index" = crypt$cells))
          self$cells[indices,tracecol][order(self$cells[indices,"Index"])] <- dists[order(crypt$cells)]
        }
      }
    },

    traceAngleAlongFeature = function(feature="Crypt"){
      # Calculate total distance along cells in a feature from origin position by tracing along edges.
      tracecol <- paste0("Angle.Along.",feature)
      self$cells[,tracecol] <- NA
      for(i in na.omit(unique(self$cells$Crypt))){
        crypt <- self$crypts[[i]]
        if(!is.null(crypt$Image)){
          img <- as.character(self$crypts[[i]]$Image)
          angles <- pi - private$rolling_angle(crypt$cells,self$CCFiles[[img]])
          indices <- self$which_cells(list("Image" = img, "Index" = crypt$cells))
          self$cells[indices,tracecol][order(self$cells[indices,"Index"])] <- angles[order(crypt$cells)]
        }
      }
    },

    add_attributes_with_key = function(key_attribute, keydb, key_column = "key", whichcells = NA, is.new=TRUE){
      # A function to add attributes in bulk with a key table.
      # key_attribute : a string specifying a source column from self$cells table
      # keydb : a table containing attribute info, with column titles. May be a string pointing to a csv file.
      # key_column : the column name within keydb that contains key information.
      # whichcells : a list of parameters to constrain attribute assignment to, as per self$which_cells function.
      # is.new : If TRUE, will add a metadata entry to parse future entries automatically.

      # Input error checking
      if(is.character(keydb)){
        tryCatch(
          keydb <- read.csv(keydb,header = T),
          error = function(err){message(paste0("Key database could not be parsed as data frame: /n",err,"/n"))})
      }

      if(!is.data.frame(keydb)){
        return(message("Key database could not be parsed as data frame."))
      }

      if(!any(colnames(keydb) == key_column)){
        return(message(paste0("No columns in key dataframe specified as '", key_column,"'.")))
      }

      # Experiment_regex = NA by default because may fall in key_regex
      if(is.na(whichcells)){
        cell_indices <-  1:nrow(self$cells)
      }else{
        cell_indices <-  self$which_cells(whichcells)
      }

      # For future use, extract all other column names aside from key.
      valnames <- colnames(keydb)[!(colnames(keydb) %in% key_column)]

      # Coerce all key columns into factors.
      for(i in names(keydb)) { keydb[,i] <- factor(keydb[,i]) }

      # Add new factor levels to new columns
      for(i in valnames) {
        # Add column for new attribute types
        if(!(i %in% colnames(self$cells))) { self$cells[,i] <- NA }
        # Coerce target column to factor if not already
        if(class(self$cells[,i]) != "factor") { self$cells[,i] <- factor(self$cells[,i]) }
        #Add new factor levels to column.
        levels(self$cells[,i]) <- c( levels(self$cells[,i]), levels(keydb[,i]) )
      }

      for(i in 1:nrow(keydb)){
        key <- keydb[,key_column][i]
        key_indices <- which(self$cells[,key_attribute] %in% key)
        indices <- key_indices[which(key_indices %in% cell_indices)]
        self$cells[indices,valnames] = keydb[i,valnames]
      }

      if(!("Unblinding" %in% names(self$metadata)) ){
        self$metadata$Unblinding = list()
      }

      if(is.new){
        i <- length(self$metadata$Unblinding) + 1
        self$metadata$Unblinding[[i]] <- list("whichcells" = whichcells,
                                              "keydb" = keydb,
                                              "key_attribute" = key_attribute,
                                              "key_column" = key_column)
      }
    },

    aggregate_by_feature = function(group.by,measure.by,melt.data=TRUE){
      # Uses table functions to aggregate data into a dataframe by grouping attributes

      ### Fcn Start
      all.by <- c( group.by, measure.by )

      ### Error checking
      if(any(!(all.by %in% colnames(self$cells)))){
        return(message(paste("Error in aggregate_by_feature: Specified column names not found in cells dataframe: ",
                             c(all.by)[!(c(all.by) %in% colnames(self$cells))])))
      }

      val.table <- table(self$cells[,all.by])
      val.key <- which(rep(
        data.frame(margin.table( val.table, 1:length(group.by) ))$Freq > 0,
        times=(2**length(measure.by))
      ))
      val.df <- data.frame(val.table)[val.key,]
      ratio.df <- data.frame(prop.table( val.table, which(group.by %in% all.by) ))[val.key,]
      val.df[, paste(c("ratio.vs",group.by), collapse = "." ) ] <- ratio.df$Freq
      colnames(val.df)[which(colnames(val.df)=="Freq")] <- "Count"
      colnames(val.df)[1:length(all.by)] <- all.by

      if(melt.data){
        # val.df <- val.df[which(base::rowSums( val.df[,measure.by]
        #                                       == matrix( rep(melt.by,times=nrow(val.df)),ncol=length(melt.by),byrow = TRUE )
        # )==length(melt.by)),]
        # Creating boolean permutation array for tagging responses.
        perms <- !(expand.grid(rep(list(c(TRUE,FALSE)),length(measure.by))))
        tags <- rep(0,nrow(perms))
        for(i in 1:nrow(perms)){
          tags[i] <- paste0(c("0.","1.")[perms[i,]+1],1:length(measure.by),collapse="|")
        }
        # Combining to factors based on boolean response.
        val.df$MarkerCombo <- 0
        for(i in 1:length(measure.by)){
          val.df$MarkerCombo <- val.df$MarkerCombo + (as.logical(val.df[,measure.by[i]]) * 2**(i-1))
        }
        val.df$MarkerCombo <- factor(val.df$MarkerCombo, levels = 0:(2**length(measure.by) - 1),labels=tags)
        # Adding labels for markers for melting.
        for(i in 1:length(measure.by)){
          val.df[,paste0("Marker",i)] <- measure.by[i]
        }
        # Tidying df order
        val.df <- val.df[,c(colnames(val.df)[1:length(group.by)],
                            "MarkerCombo",
                            paste0("Marker",1:length(measure.by)),
                            colnames(val.df)[(length(all.by)+1):(length(all.by)+2)])]
      }
      rownames(val.df) <- NULL
      return(val.df)
    },

    aggregate_by_features = function(group.by, measure.by, melt.data=TRUE){
      # Aggregates data based on grouping variables, using a list of features or markers to measure by.
      # Wrapper for self$aggregate_by_feature
      if(any(!(c(group.by,unlist(measure.by)) %in% colnames(self$cells)))){
        return(message(paste("Error in aggregate_by_features: Specified column names not found in cells dataframe: ",
                             c(group.by,unlist(measure.by))[!(c(group.by,unlist(measure.by)) %in% colnames(self$cells))])))
      }
      if(class(measure.by) != "list" ){
        return(message("Error in aggregate_by_features: measure.by parameter must be a list of markers/attributes."))
      }

      to_return <- NULL
      for(i in 1:length(measure.by)){
        to_add <- NULL
        tryCatch({
          to_add <- self$aggregate_by_feature(group.by, measure.by[[i]], melt.data)
          if(is.null(to_return)){
            to_return = to_add
          }else{
            to_return = private$rbind_fill(to_return, to_add)
          }
        },
        error = function(e){message(paste("Error in aggregate_by_features: ",e))},
        warning = function(w){message(paste("Warning in aggregate_by_features: ",w))})
      }
      return(to_return)
    }
    #
    # setMarkerColors = function(list=list()){
    #   # Provide a list of Marker names color names or hexcode RGB values, to format to color.
    #   # Will fill any unassigned markers by Brewer paired palette.
    #   # Form of list is Marker="MarkerColor", ex list=list(Reg4="red",Lgr5="#00EE00")
    #   brewer.paired <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928")
    #   if(!any("MarkerColors" %in% names(self$metadata))){
    #     self$metadata[["MarkerColors"]] <- brewer.paired
    #   }
    #
    # }
  ),
  private = list(
    angle = function(x,y){
      # Calcs angle between two vectors in radians, nondirectional max value pi. Assumes origin of 0^n
      dot.prod <- x%*%y
      norm.x <- base::norm(x,type="2")
      norm.y <- base::norm(y,type="2")
      theta <- acos(round(dot.prod / (norm.x * norm.y),digits = 15))
      return(as.numeric(theta))
    },

    twolines_angle = function(a,b,c,d){
      # Finds the angle between two lines AB and CD
      # Assumes that all connect A-B-C-D, with B and C serving as origins for A and D.
      x <- NaN
      try(x <- private$angle(a-b,d-c),silent=TRUE)
      if(is.nan(x)){return(pi)}
      else{return(x)}
    },

    rolling_angle = function(v,cc){
      tryCatch({
        return(mapply(function(x1,y1,x2,y2,x3,y3,x4,y4){private$twolines_angle(c(x1,y1),c(x2,y2),c(x3,y3),c(x4,y4))},
                      x1 = cc$cells[v[c(1,1:(length(v)-1))],c("X")],
                      y1 = cc$cells[v[c(1,1:(length(v)-1))],c("Y")],
                      x2 = cc$cells[v[c(1:length(v))],c("X")],
                      y2 = cc$cells[v[c(1:length(v))],c("Y")],
                      x3 = cc$cells[v[c(2:length(v),length(v)+1)],c("X")],
                      y3 = cc$cells[v[c(2:length(v),length(v)+1)],c("Y")],
                      x4 = cc$cells[v[c(3:length(v),length(v)+1,length(v)+2)],c("X")],
                      y4 = cc$cells[v[c(3:length(v),length(v)+1,length(v)+2)],c("Y")]))},
        # Including error/warning functions to send back vector of 1s on assumption that supplied v is too short to work.
        error = function(e){return(rep(1,times=length(v)))},
        warning = function(w){return(rep(1,times=length(v)))})
    },

    rolling_dist = function(v,cc){
      tryCatch({
        return(mapply(function(x1,y1,x2,y2){sqrt((x1-x2)**2 + (y1-y2)**2)},
                      x1 = cc$cells[v[c(1:length(v))],c("X")],
                      y1 = cc$cells[v[c(1:length(v))],c("Y")],
                      x2 = cc$cells[v[c(2:length(v),length(v)+1)],c("X")],
                      y2 = cc$cells[v[c(2:length(v),length(v)+1)],c("Y")]))},
        # Including error/warning functions to send back vector of 1s on assumption that supplied v is too short to work.
        error = function(e){return(rep(1,times=length(v)))},
        warning = function(w){return(rep(1,times=length(v)))})
    },

    cycle_vector = function(v,i){
      # Cycles vector so that value i is at position 1 (and pos(last) if looped)
      if(v[1]==i){
        return(v)
      }else if(v[1]==v[length(v)]){
        return(c(v[which(v==i):length(v)],v[2:(which(v==i)-1)],i))
      }else{
        return(c(v[which(v==i):length(v)],v[1:(which(v==i)-1)]))
      }
    },

    connections = function(i,x){
      # Returns all vertices in set x connected to vertices i
      return(unique(unlist(lapply(i,function(y){lapply(private$which_vert(y,x),function(z){sum(x[z,])-y})}))))
    },

    which_vert = function(i,x){
      # Returns which vertices in set x contain indices i. Assumes vertices are rows.
      return(which(rowSums(matrix(x %in% i,ncol=ncol(x)))>0))
    },

    check_looping = function(v){
      # Determines whether vector loops or not based on if v[1]==v
      return(v[1]==v[length(v)])
    },

    list_append = function(list,v){
      # Robustly appends vector to list as single element, at last position.
      list[[length(list)+1]] <- v
      return(list)
    },

    create_triangulation = function(cc){
      # Create Delaunay triangulation object from CellCounter object, and modify with additional useful information regarding likely contours.
      tri <- RTriangle::triangulate(RTriangle::pslg(cc$cells[,c("X","Y")]))
      # Generate all possible combos of triangles from permutations of Tri$T, which describes all directional connections.
      tri$T.perm <- rbind(tri$T[,c(1,2,3)],tri$T[,c(1,3,2)],tri$T[,c(2,1,3)],tri$T[,c(2,3,1)],tri$T[,c(3,1,2)],tri$T[,c(3,2,1)])
      # Generate angle made between points as represented by columns 1-3-2.
      tri$T.theta <- mapply(function(x1,y1,x2,y2)
        {as.numeric(acos(round(c(x1,y1)%*%c(x2,y2) / (base::norm(c(x1,y1),type="2") * base::norm(c(x2,y2),type="2")),digits=15)))},
        cc$cells[tri$T.perm[,1],"X"]-cc$cells[tri$T.perm[,3],"X"],
        cc$cells[tri$T.perm[,1],"Y"]-cc$cells[tri$T.perm[,3],"Y"],
        cc$cells[tri$T.perm[,2],"X"]-cc$cells[tri$T.perm[,3],"X"],
        cc$cells[tri$T.perm[,2],"Y"]-cc$cells[tri$T.perm[,3],"Y"],
        USE.NAMES = FALSE)
      # Generate max(theta) for all possible connections to 1,2 above (respect to tri$E).
      # This represents a beta threshold respective to a beta skeleton.
      tri$E.beta <- mapply(function(x,y)
        {max(tri$T.theta[which((tri$T.perm[,1]==x) & (tri$T.perm[,2]==y))])},
        tri$E[,1],
        tri$E[,2],
        USE.NAMES=FALSE)
      # Generate distance of each edge.
      tri$E.dist <- mapply(function(x,y){base::sqrt(x**2 + y**2)},
                           cc$cells[tri$E[,1],"X"]-cc$cells[tri$E[,2],"X"],
                           cc$cells[tri$E[,1],"Y"]-cc$cells[tri$E[,2],"Y"],
                           USE.NAMES = FALSE)
      cc$metadata[["triangulation"]] <- tri
    },

    string_extract = function(string,regex){
      # Returns instance of regex match from string.
      return(base::regmatches(string, base::regexpr(regex, string)))
    },

    rollcount = function(v.stat,vector_set,decreasing=FALSE){
      # Gives a rolling count of instances of vertices per vector, sorted by vector statistic v.stat.
      # Expects v.stat to be a numerical vector of length = rows vector_set. vector_set must be a data.frame or matrix with vectors expressed row-wise.
      array_flat <- as.vector(t(vector_set[order(v.stat, decreasing=decreasing),]))
      return(matrix(ave(
                      rep(TRUE,length(array_flat)),
                      array_flat,FUN=cumsum),
                    nrow=nrow(vector_set),
                    ncol=ncol(vector_set),
                    byrow=TRUE)[order(order(v.stat, decreasing=decreasing)),])
    },

    pare_edges_by_rolling = function(rollset, connections.max, threshold=0){
      # Returns indices of rollset after pruning by a rolling set of connections.
      # Rollset is a casting of a vector set with connections per vertex ordered by a test stat as per private$rollcount.
      # Threshold defines the minimum number of connections per vector needed to satisfy test. Must be less than ncol(rollset)
      return(which(base::rowSums(rollset <= connections.max) > threshold))
    },

    pare_edges_by_stat = function(stat, edgeset, decreasing=TRUE, connections.max){
      # Returns indices of edgeset after pruning based on stat and connections.max.
      # By default retains lowest stat (decreasing=TRUE)
      indices <- rep(TRUE,nrow(edgeset))
      too_many <- which(table(as.vector(edgeset))>connections.max)
      for(i in too_many){
        to_test <- which(base::rowSums(edgeset==i)>0)
        to_remove <- order(stat[to_test],decreasing=decreasing)
        if(length(to_remove)>connections.max){
          indices[to_test[to_remove[1:(length(to_test)-connections.max)]]] <- FALSE
        }
      }
      return(which(indices))
    },

    pare_edges_by_count = function(edgeset,connections.max){
      # Returns indices of edgeset after pruning based on both vertices having more than connections.max.
      indices <- rep(TRUE,nrow(edgeset))
      too_many <- which(table(as.vector(edgeset))>connections.max)
      for(i in too_many){
        to_test <- private$which_vert(i,edgeset)
        # Determines indices of vertices that both have too many connections, ergo cutting connection safely.
        to_remove <- to_test[which(apply(matrix(edgeset[to_test,] %in% too_many,ncol=2,byrow=FALSE),1,FUN=prod)==1)]
        if(length(to_remove)>0){
          indices[to_remove] <- FALSE
        }
      }
      return(which(indices))
    },

    pared_edgeset_to_ordered_vectors = function(edgeset, make_looping=TRUE){
      # Converts an unordered set of edges to a list of ordered vectors. Assumes that it has been pared down to a linearizable set.
      # Output is given as looped ordered vectors.
      ordered_sets <- list()
      # Start with single connections
      for(i in which(table(as.vector(edgeset))==1)){
        # Skip if already part of a set.
        if(i %in% unlist(ordered_sets)){next}
        ordered <- c(private$connections(i,edgeset),i)
        # Works because it should be strictly linear at this point
        while(!all(private$connections(ordered[1],edgeset) %in% ordered)){
          to_try <- private$connections(ordered[1],edgeset)
          if(!all(to_try %in% ordered)){
            ordered <- append(ordered, to_try[-which(to_try %in% ordered)], after=0)
          }
        }
        # Append first to last to indicate looping. Easier cleanup later.
        if((ordered[1]!=ordered[length(ordered)]) & make_looping ){
          ordered <- c(ordered,ordered[1])
        }
        ordered_sets[[length(ordered_sets)+1]] <- ordered
      }
      # Paring sorted values from set.
      if(length(ordered_sets) > 0){
        subset <- (1:nrow(edgeset))[-private$which_vert(unlist(ordered_sets),edgeset)]
      }else{
        subset <- 1:nrow(edgeset)
      }
      # Grabbing random vectors and creating ordered sets from there.
      while(length(subset)>0){
        ordered <- c(edgeset[subset,][1,2],edgeset[subset,][1,1])
        # Works because it's a directional assessment
        while(!all(private$connections(ordered[1],edgeset[subset,]) %in% ordered)){
          to_try <- private$connections(ordered[1],edgeset[subset,])
          if(!all(to_try %in% ordered)){
            ordered <- append(ordered, to_try[-which(to_try %in% ordered)], after=0)
          }
        }
        subset <- subset[-private$which_vert(ordered,edgeset[subset,])]
        #Adding last index to first position as well to indicate looping vector.
        ordered <- append(ordered,ordered[length(ordered)],after=0)
        ordered_sets[[length(ordered_sets)+1]] <- ordered
      }
      return(ordered_sets)
    },

    isEdge = function(vector,edgeset){
      # Returns index of edge formed by vector components in edgeset, or 0 if not in edgeset.
      tryCatch(
        {
          index = (base::rowSums((edgeset==vector[1])+(edgeset==vector[2]))==2)
          if(any(index)){
            return(which(index))
          }
          else{
            return(0)
          }
        },
        error = function(e){message("Error in isEdge: vector must be length 2.")},
        warning = function(w){}
      )
    },

    areEdges = function(m,edgeset){
      # Returns index of edge formed by vector components in each row of matrix m in edgeset, or 0 if not in edgeset.
      # m is a matrix with 2 columns.
      # A mapply wrapper for isEdge.
      tryCatch(
        {
          return(mapply(function(v1,v2){private$isEdge(c(v1,v2),edgeset)},
                        v1=m[,1],
                        v2=m[,2]))
        },
        error = function(e){message("Error in areEdges: Matrix not compatible.")},
        warning = function(w){}
      )
    },

    vector_to_edges = function(v){
      # Returns a matrix of ordered edges, with edges moving between column 1 to column 2 down through rows.
      # v is an ordered vector of vertex indices.
      tryCatch(
        {
          return(matrix(c(v[1:(length(v)-1)],v[2:length(v)]),ncol=2))
        },
        error = function(e){message("Error in vector_to_edges: Vector not compatible.")},
        warning = function(w){}
      )
    },

    rbind_fill = function(d1,d2){
      # Returns merged dataframes with empty columns in either filled with NA. Sets column class with priority for d1
      d1.names <- names(d1)
      d2.names <- names(d2)
      d1.names.u <- d1.names[!(d1.names %in% d2.names)]
      d2.names.u <- d2.names[!(d2.names %in% d1.names)]
      names.all <- c(d1.names,d2.names.u)
      d2.class <- sapply(d2, class)
      for(name in d2.names.u){
        d1[,name] <- NA
        if(d2.class[name]=="factor"){
          d1[,name] <- factor(d1[,name])
        }else{
          class(d1[,name]) <- d2.class[name]
        }
      }
      d2[,d1.names.u] <- NA
      d1.class <- sapply(d1,class)
      for(name in names.all){
        if(d1.class[name]=="factor"){
          d2[,name] <- factor(d2[,name])
        }else{
          class(d2[,name]) <- d1.class[name]
        }
      }
      return(rbind(d1,d2))
    },

    stats_per_feature = function(feature="Crypt"){
      # TODO: create aggregate stats around specified feature, further functionality, with specified table
    },

    find_cut_point = function(v,img,dist.mean,dist.sd){
      # Returns the optimal cut point in an ordered vector
      # Inputs:
      # v = ordered vector describing points.
      # cc = cellcounter object
      # edges = edges E from a triangulation object
      # dist.mean = mean distance between points
      # dist.sd = sd of distance between points
      cc <- self$CCFiles[[img]]
      testfwd <- log(1:length(v) + 0.1)
      testrev <- log(length(v):1 + 0.1)
      angular_strain <- pi - private$rolling_angle(v,cc)
      dist_strain <- exp((private$rolling_dist(v,cc)-dist.mean)/dist.sd)
      # Checks that edges are present in Delaunay triangulation, applies strong penalty if not.
      edgecheck <- (c(private$areEdges(private$vector_to_edges(v),cc$metadata$triangulation$E),1)==0)*100000
      test_stat <- testfwd*testrev*angular_strain*angular_strain*dist_strain + edgecheck
      return(order(test_stat,decreasing=TRUE)[1])
    },

    cut_open_crypt = function(v,img,dist.mean,dist.sd){
      i0 <- self$CCFiles[[img]]$metadata$OriginCells[self$CCFiles[[img]]$metadata$OriginCells %in% v]
      # Cycle vector to origin and make loop.
      v <- private$cycle_vector(v,i0)
      if(v[1]!=v[length(v)]){
        v <- c(v,v[1])
      }
      # Test to verify optimal opening point for crypt.
      cut_point <- private$find_cut_point(v,img,dist.mean,dist.sd)
      v <- c(v[(cut_point+1):(length(v)-1)],v[1:cut_point])
      return(v)
    }

  )#,
  # active = list(
  #
  #   Markers = function(v){
  # TODO: add dynamic marker element.
  #     # Returns either all unique markers in CCfiles, or which CCFiles have all specified markers
  #     if(missing(v)){ return(unique(unlist(lapply(self$crypts, function(x){x$Markers})))) }
  #     else{
  #       return(which(lapply(self$CCFiles, function(x){ v %in% x$Markers})))
  #     }
  #   }
  # )
)
