# List of Auxiliary Functions ---------------------------------------------
# (1) check_input     : check symmetric matrix / dist object
# (2) sample_group    : given a vector of (positive) integers, sample random groupings
# (3) label2idlist    : given a label vector, return a list containing id of each class
# (4) dist_arrange    : given distance matrix and label, return an arranged distance mat with vec.ng
# (5) arrange_2sample : rearrange the distance matrix according to the label
# (6) aux_pseudomean  : distance to pseudomean from 1st to other observations

# (1) check_input ---------------------------------------------------------
#' @keywords internal
#' @noRd
check_input <- function(x){
  if (inherits(x, "dist")){
    x = as.matrix(x)
  } else {
    if ((!is.matrix(x))||(nrow(x)!=ncol(x))||(!isSymmetric(x))){
      stop("* CHECK_INPUT : input should be a symmetric matrix.")
    }
  }
  if (any(abs(diag(x)) > .Machine$double.eps)){
    stop("* CHECK_INPUT : diagonal elements should all be zeros.")
  }
  return(x)
}


# (2) sample_group --------------------------------------------------------
#' @keywords internal
#' @noRd
sample_group <- function(intvec){
  n = sum(intvec)
  G = length(intvec)
  output = list()
  allvec = 1:n
  for (g in 1:(G-1)){
    samid = sample(allvec, intvec[g])
    output[[g]] = samid
    allvec = setdiff(allvec, samid)
  }
  output[[G]] = allvec
  return(output)
}


# (3) label2idlist --------------------------------------------------------
#' @keywords internal
#' @noRd
label2idlist <- function(label){
  ulabel = unique(label)
  nulabs = length(ulabel)
  output = list()
  for (i in 1:nulabs){
    output[[i]] = which(label==ulabel[i])
  }
  return(output)
}


# (4) dist_arrange --------------------------------------------------------
#' @keywords internal
#' @noRd
dist_arrange <- function(dmat, label){
  idx = label2idlist(label)
  idxnum = unlist(idx)
  
  out.ng = unlist(lapply(idx, length))
  out.dd = dmat[idxnum, idxnum]
  
  output = list()
  output$dmat   = out.dd
  output$vec.ng = out.ng
  return(output)
}

# (5) arrange_2sample -----------------------------------------------------
#' @keywords internal
#' @noRd
arrange_2sample <- function(x, label){
  output = list()
  ulabel = unique(label)
  id1 = which(label==ulabel[1])
  id2 = which(label==ulabel[2])
  
  hormat1 = cbind(x[id1,id1], x[id1,id2])
  hormat2 = cbind(x[id2,id1], x[id2,id2])
  outmat  = rbind(hormat1, hormat2)
  m       = length(id1)
  n       = length(id2)
  
  output = list()
  output$mat = outmat
  output$m   = m         # size of class 1
  output$n   = n         # size of class 2
  return(output)
}

# (6) aux_pseudomean ------------------------------------------------------
#' @keywords internal
#' @noRd
aux_pseudomean <- function(dmat){
  # we need embedding .. umm .. automatic dimension selection
  if (nrow(dmat)==1){
    stop("* aux_pseudomean : error..")
  } else if (nrow(dmat)==2){
    return(dmat[1,2])
  } else {
    embedded = aux_pseudomean_auto(dmat)
    n = nrow(embedded)
    p = ncol(embedded)
    
    # centering based on other points
    emcenter = as.vector(base::colMeans(embedded[2:n,]))
    embednew = embedded - matrix(rep(emcenter,n), ncol=p, byrow=TRUE)
    
    # compute scalar
    d1mat = dmat[2:n,2:n]                          # d(x,y)
    d2mat = as.matrix(stats::dist(embednew[2:n,])) # ||x-y||
    d12mat = (d1mat*d2mat)
    d22mat = (d2mat^2)
    dlower = base::lower.tri(d12mat)
    cstar =sum(d12mat[dlower])/sum(d22mat[dlower])
    
    # update embednew and compute 
    erow1 = cstar*as.vector(embednew[1,])
    return(sqrt(sum(erow1^2))) 
  }
}
#' @keywords internal
#' @noRd
aux_pseudomean_auto <- function(dmat){ # only positive eigenvalues' part
  n = nrow(dmat)
  J = diag(rep(1,n))-(1/n)*outer(rep(1,n),rep(1,n))
  B = -(J%*%(dmat^2)%*%J)/2.0
  eigB = base::eigen(B, symmetric = TRUE) # decreasing order
  
  m = max(length(which(eigB$values > 0)),2)
  X = (eigB$vectors[,1:m])%*%(base::diag(sqrt(eigB$values[1:m])))
  return(X)
}