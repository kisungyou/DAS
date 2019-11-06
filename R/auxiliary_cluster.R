############# Auxiliary R Functions #########################
# (08) rclust_index_Sw      : within-cluster
#      rculst_index_Sb      : between-cluster sum of distances given a distance matrix
# (09) rclust_index_Ns      : number of counts information
# (10) rclust_membership    : membership of clusterings
# (11) rclust_concordant    : s(+) and s(-)
# (12) rclust_Z             : membership matrix of size (n x q)
# (13) rclust_index_density : use for SDbw method


# (08) sum of distances ---------------------------------------------------
#   rclust_index_Sw : within-cluster
#   rclust_index_Sb : between-cluster sum of distances given a distance matrix
#' @keywords internal
#' @noRd
rclust_index_Sw <- function(distmat, label, mfdname, memberlist){
  ulabel = unique(label)
  q      = length(ulabel)
  output = 0
  for (i in 1:q){
    idlabel = memberlist[[which(ulabel==ulabel[i])]]
    if (length(idlabel)!=1){
      partmat = distmat[idlabel,idlabel]
      output  = output + sum(as.vector(partmat[lower.tri(partmat)]))
    }
  }
  return(output)
}
#' @keywords internal
#' @noRd
rclust_index_Sb <- function(distmat, label, mfdname, memberlist){
  ulabel = unique(label)
  q      = length(ulabel)
  score  = 0.0
  for (k in 1:(q-1)){
    label.k = memberlist[[which(ulabel==ulabel[k])]]
    for (j in (k+1):q){
      label.j = memberlist[[which(ulabel==ulabel[j])]]
      score   = score + sum(as.vector(distmat[label.k, label.j]))
    }
  }
  return(score)
}

# (09) rclust_index_Ns    : number of counts information
#' @keywords internal
#' @noRd
rclust_index_Ns <- function(label, memberlist){
  # basic parameters
  ulabel = unique(label)
  n = length(label)
  q = length(ulabel)
  # Nt
  Nt = as.integer(n*(n-1)/2)
  # Nw
  Nw = 0
  for (i in 1:q){
    Nk = length(memberlist[[which(ulabel==ulabel[i])]])
    Nw = as.integer(Nk*(Nk-1)/2) + Nw
  }
  # Nb
  Nb = Nt - Nw
  # return output
  output = list()
  output$Nt = Nt
  output$Nb = Nb
  output$Nw = Nw
  return(output)
}

# (10) rclust_membership  : membership of clusterings
#' @keywords internal
#' @noRd
rclust_membership <- function(label){
  ulabel = sort(unique(label))
  memvec = list()
  for (i in 1:length(ulabel)){
    memvec[[i]] = which(label==ulabel[i])
  }
  return(memvec)
}

# (11) rclust_concordant  : s(+) and s(-)
#' @keywords internal
#' @noRd
rclust_concordant <- function(distmat, label, membership){
  # get some parameters
  q      = length(membership)
  ulabel = unique(label)
  # 1. compute all within-distance vector
  vec.within = c()
  for (i in 1:q){
    idx = membership[[which(ulabel==ulabel[i])]]
    partmat = distmat[idx,idx]
    vec.within = c(vec.within, as.vector(partmat[upper.tri(partmat)]))
  }
  # 2. compute all between-distance vector
  vec.between = c()
  for (i in 1:(q-1)){
    id1 = membership[[which(ulabel==ulabel[i])]]
    for (j in (i+1):q){
      id2 = membership[[which(ulabel==ulabel[j])]]
      partmat = distmat[id1,id2]
      vec.between = c(vec.between, as.vector(partmat))
    }
  }
  # 3. count all numbers
  niter = length(vec.within)
  count.con = 0
  count.dis = 0
  for (i in 1:niter){
    count.con = count.con + sum((vec.between>vec.within[i]))
    count.dis = count.dis + sum((vec.between<vec.within[i]))
  }
  # return
  output = list()
  output$con = count.con
  output$dis = count.dis
  return(output)
}

# (12) rclust_Z           : membership matrix of size (n x q)
#' @keywords internal
#' @noRd
rclust_Z <- function(label, membership){
  n = length(label)
  q = length(membership)
  ulabel = unique(label)
  
  output = array(0,c(n,q))
  for (i in 1:q){
    id = membership[[which(ulabel==ulabel[i])]]
    output[id,i] = 1
  }
  return(output)
}

# (13) rclust_index_density : use for SDbw method in extrinsic manner
#' @keywords internal
#' @noRd
rclust_index_density <- function(meanvec, datamat, thr){
  centered = datamat - matrix(rep(meanvec,nrow(datamat)),nrow=nrow(datamat),byrow=TRUE)
  distvec  = apply(centered, 1, function(x){sqrt(sum(x^2))})
  return(sum(distvec<thr))
}