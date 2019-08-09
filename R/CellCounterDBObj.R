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
    # Markers = NA,
    # MarkerColors = NA,

    initialize = function(filelist){
      self$CCFiles <- list()
      self$crypts <- list()
      self$metadata <- list()
      for(file in filelist){
        self$addCC(file)
      }
    },

    addCC = function(ccfile){
      # '
      # Adds new CellCounterObj to list, then adds cells to database data.frame.
      # Need method to run a list of files through, potential integration with other methods for post-processing.
      # '
      cc <- ccproc::CellCounter$new(ccfile)
      private$create_triangulation(cc)
      image_name <- as.character(cc$metadata$Image_Filename)
      if(sum(image_name %in% names(self$CCFiles)) == 0){
        self$CCFiles[[image_name]] <- cc
        newcells <- cbind(cc$cells,Image=rep(image_name,nrow(cc$cells)),Index=1:nrow(cc$cells))
        if(any("SampleRegex" %in% names(self$metadata))){
          newcells$Sample <-factor(private$string_extract(as.character(newcells$Image),self$metadata$SampleRegex)[1])
        }
        if(all(is.na(self$cells))){
          self$cells <- newcells
        }else{
          self$cells <- plyr::rbind.fill(self$cells, newcells)
        }
      }
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

    processFeatures = function(feature="Crypt"){
      # '
      # Process loaded cells into features. Expand in future to non-crypt stuff.
      # '

      if(!(feature %in% colnames(self$cells))){
        self$cells[,feature] <- 0
        self$cells[,paste0("Pos.",feature)] <- NA}
      for(cc in self$CCFiles){
        tryCatch({
          if(!any("triangulation" %in% names(cc$metadata))){
            private$create_triangulation(cc)
          }
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
                testfwd <- log(1:length(subset))
                testrev <- log(length(subset):1)
                angular_strain <- pi - private$rolling_angle(subset,cc)
                dist_strain <- exp((private$rolling_dist(subset,cc)-dist.mean)/dist.sd)
                # Checks that edges are present in Delaunay triangulation, applies strong penalty if not.
                edgecheck <- (c(private$areEdges(private$vector_to_edges(subset),tri$E),1)==0)*100000
                test_stat <- testfwd*testrev*angular_strain*angular_strain*dist_strain + edgecheck
                cleavepoints <- c(cleavepoints,
                                  which(set==subset[order(test_stat,decreasing=T)[1]])[1],
                                  which(set==subset[order(test_stat,decreasing=T)[1]])[1]+1)
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
              unassigned_sets <- private$list_append(unassigned_sets,subset)
            }
          }


          # Add column for feature and position if not present.
          if(!all(c(feature,paste0("Pos.",feature)) %in% colnames(self$cells))){
            self$cells[,feature] <- NA
            self$cells[,paste0("Pos.",feature)] <- NA
          }

          # Final processing of crypt vectors and sorting of crypts into CryptObjs and addition to cell database.
          for(v in crypts){
            i <- length(self$crypts) + 1
            i0 <- cc$metadata$OriginCells[cc$metadata$OriginCells %in% v]
            # Cycle vector and make loop.
            v <- private$cycle_vector(v,i0)
            if(v[1]!=v[length(v)]){
              v <- c(v,v[1])
            }
            # Re-testing to verify optimal opening point for crypt.
            testfwd <- log(1:length(v))
            testrev <- log(length(v):1)
            angular_strain <- pi - private$rolling_angle(v,cc)
            dist_strain <- exp((private$rolling_dist(v,cc)-dist.mean)/dist.sd)
            # Checks that edges are present in Delaunay triangulation, applies strong penalty if not.
            edgecheck <- (c(private$areEdges(private$vector_to_edges(v),tri$E),1)==0)*100000
            test_stat <- testfwd*testrev*angular_strain*dist_strain + edgecheck
            cut_point <- order(test_stat,decreasing=TRUE)[1]
            v <- c(v[(cut_point+1):(length(v)-1)],v[1:cut_point])
            c <- ccproc::Crypt$new(origin=i0,cells=v,Image=cc$metadata$Image_Filename,E.index=private$areEdges(private$vector_to_edges(v),tri$E))
            self$crypts[[i]] <- c
            indices <- self$which_cells(list("Image" = c$Image, "Index" = v))
            self$cells[indices,feature] <- i
            self$cells[indices,paste0("Pos.",feature)] <- c$pos(self$cells$Index[indices])
          }
        },
      error = function(err){message(paste("Error: Unable to process ",cc$metadata$Image_Filename))},
      warning = function(w){})
      }
    },

    plot_debug = function(images,savedir="debug_images/"){
      for(image in images){
        tryCatch({

          subset <- self$cells[self$which_cells(list("Image"=image)),]
          p <- ggplot(subset[order(subset$Pos.Crypt),], aes(x=X,y=Y,group=Crypt,label=abs(Pos.Crypt),fill=as.factor(Crypt))) +
            geom_path() +
            geom_point(shape=21,color="black",size=6,stroke=1) +
            geom_text() +
            scale_x_reverse() + scale_y_reverse() +
            theme(axis.line=element_blank(),
                  axis.text.x=element_blank(),
                  axis.text.y=element_blank(),
                  axis.ticks=element_blank(),
                  axis.title.x=element_blank(),
                  axis.title.y=element_blank(),
                  legend.position="none",
                  panel.background=element_blank(),
                  panel.border=element_blank(),
                  panel.grid.major=element_blank(),
                  panel.grid.minor=element_blank(),
                  plot.background=element_blank())

          png(paste0(savedir, image," debug plot.png"),
              width=(max(subset$X)/as.numeric(levels(self$CCFiles[[image]]$metadata$X_Calibration)[1]) -
                       min(subset$X)/as.numeric(levels(self$CCFiles[[image]]$metadata$X_Calibration)[1]) + 50),
              height=(max(subset$Y)/as.numeric(levels(self$CCFiles[[image]]$metadata$Y_Calibration)[1]) -
                        min(subset$Y)/as.numeric(levels(self$CCFiles[[image]]$metadata$Y_Calibration)[1]) + 50),
              res=120, type="cairo")
          plot(p)
          dev.off()
        },
        error = function(err){message(paste0("Could not plot ",image))})
      }
    },

    which_cells = function(listconditions){
      # Returns indices in self$cells that match conditions specified in list.
      return(which(apply(sapply(names(listconditions),function(x){self$cells[,x] %in% listconditions[[x]] }),1,prod)>0))
    },

    tag_attribute = function(Regex="",attribute=""){
      # Inserts a column into dataframe as attribute, labeling cells by regex pattern as applied to Image string.
      if(Regex==""){
        return(message("Error: Attribute tagging requires a RegEx pattern to apply to image names!"))
      }
      if(attribute=="" | any(attribute %in% colnames(self$cells))){
        return(message("Error: Attribute must be named and not already in use in cells table."))
      }
      self$metadata[[paste0(attribute,"Regex")]] <- Regex
      self$cells[,attribute] <- factor(private$string_extract(as.character(self$cells$Image),self$metadata[[paste0(attribute,"Regex")]])[1])
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
      theta <- acos(dot.prod / (norm.x * norm.y))
      return(as.numeric(theta))
    },

    twolines_angle = function(a,b,c,d){
      # Finds the angle between two lines AB and CD
      # Assumes that all connect A-B-C-D, with B and C serving as origins for A and D.
      x <- NaN
      try(x <- private$angle(a-b,d-c),silent=TRUE)
      if(is.nan(x)){return(pi)}
      else(return(x))
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
        {as.numeric(acos(c(x1,y1)%*%c(x2,y2) / (base::norm(c(x1,y1),type="2") * base::norm(c(x2,y2),type="2"))))},
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
    }
  )#,
  # active = list(
  #
  #   Markers = function(v){
  #     # Returns either all unique markers in CCfiles, or which CCFiles have all specified markers
  #     if(missing(v)){ return(unique(unlist(lapply(self$crypts, function(x){x$Markers})))) }
  #     else{
  #       return(which(lapply(self$CCFiles, function(x){ v %in% x$Markers})))
  #     }
  #   }
  # )
)
