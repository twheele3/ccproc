#
#
#
#
#
#
# Miscellaneous functions for algorithm development.
#
#
#
#

angle <- function(x,y){
  # Calcs angle between two vectors in radians, nondirectional max value pi. Assumes origin of 0^n
  dot.prod <- x%*%y
  norm.x <- base::norm(x,type="2")
  norm.y <- base::norm(y,type="2")
  theta <- acos(round(dot.prod / (norm.x * norm.y),digits = 15))
  return(as.numeric(theta))
}

twolines_angle <- function(a,b,c,d){
  # Finds the angle between two lines AB and CD
  # Assumes that all connect A-B-C-D, with B and C serving as origins for A and D.
  x <- NaN
  try(x <- ccproc::angle(a-b,d-c),silent=TRUE)
  if(is.nan(x)){return(pi)}
  else{return(x)}
}

rolling_angle <- function(v,cc){
  tryCatch({
    return(mapply(function(x1,y1,x2,y2,x3,y3,x4,y4){ccproc::twolines_angle(c(x1,y1),c(x2,y2),c(x3,y3),c(x4,y4))},
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
}

rolling_dist <- function(v,cc){
  tryCatch({
    return(mapply(function(x1,y1,x2,y2){sqrt((x1-x2)**2 + (y1-y2)**2)},
                  x1 = cc$cells[v[c(1:length(v))],c("X")],
                  y1 = cc$cells[v[c(1:length(v))],c("Y")],
                  x2 = cc$cells[v[c(2:length(v),length(v)+1)],c("X")],
                  y2 = cc$cells[v[c(2:length(v),length(v)+1)],c("Y")]))},
    # Including error/warning functions to send back vector of 1s on assumption that supplied v is too short to work.
    error = function(e){return(rep(1,times=length(v)))},
    warning = function(w){return(rep(1,times=length(v)))})
}

cycle_vector = function(v,i){
  # Cycles vector so that value i is at position 1 (and pos(last) if looped)
  if(v[1]==i){
    return(v)
  }else if(v[1]==v[length(v)]){
    return(c(v[which(v==i):(length(v)-1)],v[1:which(v==i)]))
  }else{
    return(c(v[which(v==i):length(v)],v[1:(which(v==i)-1)]))
  }
}

connections <- function(i,x){
  # Returns all vertices in set x connected to vertices i
  return(unique(unlist(lapply(i,function(y){lapply(ccproc::which_vert(y,x),function(z){sum(x[z,])-y})}))))
}

which_vert <- function(i,x,both=FALSE){
  # Returns which vertices in set x contain indices i. Assumes vertices are rows.
  # Optional parameter both checks whether both verts in edges are in set.
  return(which(rowSums(matrix(x %in% i,ncol=ncol(x)))>both))
}

check_looping <- function(v){
  # Determines whether vector loops or not based on if v[1]==v
  return(v[1]==v[length(v)])
}

list_append <- function(list,v){
  # Robustly appends vector to list as single element, at last position.
  list[[length(list)+1]] <- v
  return(list)
}

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
  return(cc)
}

string_extract <- function(string,regex){
  # Returns instance of regex match from string.
  return(base::regmatches(string, base::regexpr(regex, string)))
}

rollcount <- function(v.stat,vector_set,decreasing=FALSE){
  # Gives a rolling count of instances of vertices per vector, sorted by vector statistic v.stat.
  # Expects v.stat to be a numerical vector of length = rows vector_set. vector_set must be a data.frame or matrix with vectors expressed row-wise.
  array_flat <- as.vector(t(vector_set[order(v.stat, decreasing=decreasing),]))
  return(matrix(ave(
    rep(TRUE,length(array_flat)),
    array_flat,FUN=cumsum),
    nrow=nrow(vector_set),
    ncol=ncol(vector_set),
    byrow=TRUE)[order(order(v.stat, decreasing=decreasing)),])
}

pare_edges_by_rolling <- function(rollset, connections.max, threshold=0){
  # Returns indices of rollset after pruning by a rolling set of connections.
  # Rollset is a casting of a vector set with connections per vertex ordered by a test stat as per ccproc::rollcount.
  # Threshold defines the minimum number of connections per vector needed to satisfy test. Must be less than ncol(rollset)
  return(which(base::rowSums(rollset <= connections.max) > threshold))
}

pare_edges_by_stat <- function(stat, edgeset, decreasing=TRUE, connections.max){
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
}

pare_edges_by_count <- function(edgeset,connections.max){
  # Returns indices of edgeset after pruning based on both vertices having more than connections.max.
  indices <- rep(TRUE,nrow(edgeset))
  too_many <- which(table(as.vector(edgeset))>connections.max)
  for(i in too_many){
    to_test <- ccproc::which_vert(i,edgeset)
    # Determines indices of vertices that both have too many connections, ergo cutting connection safely.
    to_remove <- to_test[which(apply(matrix(edgeset[to_test,] %in% too_many,ncol=2,byrow=FALSE),1,FUN=prod)==1)]
    if(length(to_remove)>0){
      indices[to_remove] <- FALSE
    }
  }
  return(which(indices))
}

order_edges <- function(edgeset){
  # Returns an edgeset ordered with smaller,larger indices
  return(t(mapply(function(x,y){c(min(x,y),max(x,y))},
                  x=edgeset[,1], y=edgeset[,2])))
}

next_edge <- function(ordered_edges){
  # Returns a list of ordered connections using a sorted edgeset
  return(mapply(function(c2){which(ordered_edges[,1] == c2)}, c2=ordered_edges[,2]))
}

pared_edgeset_to_ordered_vectors <- function(edgeset, make_looping=TRUE){
  # Converts an unordered set of edges to a list of ordered vectors. Assumes that it has been pared down to a linearizable set.
  # Output is given as looped ordered vectors.
  oedges <- ccproc::order_edges(edgeset)
  if(any(table(oedges)>2) | any(as.integer(names(which(table(oedges[,1])>1))) %in% oedges[,2])){
    return(message("Edge set could not be coerced into ordered vectors."))
  }
  # two tiered sort, by second then first
  oedges <- oedges[order(oedges[,2]),]
  oedges <- oedges[order(oedges[,1]),]
  # Origins of edgesets
  o.edges <- unique(oedges[,1])[which(!(unique(oedges[,1]) %in% unique(oedges[,2])))]

  # List of connections to next edge
  nextedge <- ccproc::next_edge(oedges)

  sets <- list()
  ends <- matrix(0,ncol=2,nrow=length(o.edges))
  for(i in 1:length(o.edges)){
    which.o <- which(oedges[,1] == o.edges[i])
    set <- oedges[which.o[1],1]
    for(j in 1:length(which.o)){
      set <- append(set,oedges[which.o[j],2], after = length(set)*(j==1))
      nextlink <- nextedge[[which.o[j]]]
      while(length(nextlink) == 1){
        to.add <- oedges[nextlink,2]
        set <- append(set,oedges[nextlink,2], after = length(set)*(j==1))
        nextlink <- nextedge[[nextlink]]
      }
    }
    sets[[i]] <- set
    ends[i,] <- c(set[1],set[length(set)])
  }
  while(length(as.integer(names(which(table(ends[(ends[,1]!=ends[,2]),])>1))))){
    i <- as.integer(names(which(table(ends[(ends[,1]!=ends[,2]),])>1)))[1]
    indices <- which(t(ends) == i)
    sides <- indices %% 2
    rownum <- ceiling(indices/2)
    sets[[rownum[1]]] <- c(sets[[rownum[1]]][abs(1:length(sets[[rownum[1]]]) -
                                                   (sides[1]==0)*(1+length(sets[[rownum[1]]])))],
                           sets[[rownum[2]]][abs(1:length(sets[[rownum[2]]]) -
                                                   sides[2]*(1+length(sets[[rownum[2]]])) )][2:length(sets[[rownum[2]]])])
    ends[rownum[1],] <- sets[[rownum[1]]][c(1,length(sets[[rownum[1]]]))]
    sets[[rownum[2]]] <- NULL
    ends <- matrix(c(ends[-rownum[2],]),ncol=2)
  }
  if(make_looping){
    for(i in which(ends[,1]!=ends[,2])){
      sets[[i]] <- append(sets[[i]],sets[[i]][1],after=length(sets[[i]]))
    }
  }
  return(sets)
}

isEdge <- function(vector,edgeset){
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
}

areEdges <- function(m,edgeset){
  # Returns index of edge formed by vector components in each row of matrix m in edgeset, or 0 if not in edgeset.
  # m is a matrix with 2 columns.
  # A mapply wrapper for isEdge.
  tryCatch(
    {
      return(unlist(mapply(function(v1,v2){ccproc::isEdge(c(v1,v2),edgeset)},
                    v1=m[,1],
                    v2=m[,2])))
    },
    error = function(e){message("Error in areEdges: Matrix not compatible.")},
    warning = function(w){}
  )
}

vector_to_edges <- function(v){
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

rbind_fill <- function(d1,d2){
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
}

stats_per_feature <- function(feature="Crypt"){
  # TODO: create aggregate stats around specified feature, further functionality, with specified table
}

find_cut_point <- function(v,cc,dist.mean,dist.sd){
  # Returns the optimal cut point in an ordered vector
  # Inputs:
  # v = ordered vector describing points.
  # cc = cellcounter object
  # edges = edges E from a triangulation object
  # dist.mean = mean distance between points
  # dist.sd = sd of distance between points
  v <- as.integer(unlist(v))
  testfwd <- log(1:length(v) + 0.1)
  testrev <- log(length(v):1 + 0.1)
  angular_strain <- pi - ccproc::rolling_angle(v,cc)
  dist_strain <- exp((ccproc::rolling_dist(v,cc)-dist.mean)/dist.sd)
  # Checks that edges are present in Delaunay triangulation, applies strong penalty if not.
  edgecheck <- (c(ccproc::areEdges(ccproc::vector_to_edges(v),cc$metadata$triangulation$E),1)==0)*100000
  test_stat <- testfwd*testrev*angular_strain*angular_strain*dist_strain + edgecheck
  return(order(test_stat,decreasing=TRUE)[1])
}

cut_open_crypt <- function(v,cc,dist.mean,dist.sd){
  i0 <- cc$metadata$OriginCells[cc$metadata$OriginCells %in% v]
  # Cycle vector to origin and make loop.
  v <- ccproc::cycle_vector(v,i0)
  if(v[1]!=v[length(v)]){
    v <- c(v,v[1])
  }
  # Test to verify optimal opening point for crypt.
  cut_point <- ccproc::find_cut_point(v,cc,dist.mean,dist.sd)
  v <- c(v[(cut_point+1):(length(v)-1)],v[1:cut_point])
  return(v)
}

path_check <-  function(edges, origins){
  # Checks whether all vertices in matrix edges have some pathway to one or more origin vertices
  # Returns a list of vectors of vertices, connected and disconnected respectively.
  vertices <- unique(c(edges))
  vertexcheck <- rep(FALSE,length(vertices))
  vertexcheck[which(vertices %in% origins)] = TRUE
  lastsum <- 0
  while(lastsum < sum(vertexcheck)){
    lastsum = sum(vertexcheck)
    vertexcheck[which(vertices %in% unique(c(edges[ccproc::which_vert(vertices[vertexcheck],edges),])))] = TRUE
  }
  return(list("connected" = vertices[vertexcheck],
              "disconnected" = vertices[!vertexcheck]))
}

fix_pathing <-  function(edgeset,subset,origins){
  # Reconnect disconnected pieces by first getting all affected, edges, then looking for crossover with connected set.
  checked.verts <- ccproc::path_check(edgeset[subset,],origins)
  to.reconnect1 <- ccproc::which_vert(checked.verts$disconnected, edgeset)
  to.reconnect2 <- ccproc::which_vert(checked.verts$connected, edgeset)
  edges.to.reconnect <- to.reconnect1[which(to.reconnect1 %in% to.reconnect2)]
  return( sort( c( edges.to.reconnect , subset ) ) )
}
