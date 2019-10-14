#' Cluster Validity Indices based on Pairwise Distance
#' 
#' 
#' @export
cvindex <- function(x, label, index=c("Silhouette")){
  ##################################################3
  # Check Input and Transform
  x  = check_input(x)
  n  = nrow(x)
  label = round(as.integer(label))
  if ((length(label)!=n)||(!is.vector(label))){
    stop("* cvindex : input 'label' should be a vector of length = nrow(x).")
  }
  memberlist = rclust_membership(label)
  
  # Index Issue
  allindices = c("silhouette")
  index = match.arg(tolower(index), allindices)
  
  ##################################################
  # Main Computation
  if (all(index=="silhouette")){
    return(index_silhouette(x, label, memberlist))
  }
}



# copying from RiemBaseExt ------------------------------------------------
# (01) Silhouette Index / Max / Rousseeux (1987)




# (01) Silhouette Index ---------------------------------------------------
#' @keywords internal
#' @noRd
index_silhouette <- function(x, label, memberlist){
  # compute 'pdist'
  distmat = x
  # labeling care
  ulabel = unique(label)
  k      = length(ulabel)
  if (k==1){
    stop("* cvindex : for a single-label clustering, Silhouette index is not defined.")
  }
  # let's iterate, take a lot of time !
  n = length(label)
  vec.a = rep(0,n)
  vec.b = rep(0,n)
  for (i in 1:n){ # for each data object
    # label of i-th element
    idlabel = which(ulabel==label[i])
    # compute a(i)
    idx.samelabel = setdiff(memberlist[[idlabel]],i)
    vec.a[i] = base::mean(as.vector(distmat[i,idx.samelabel]))
    # compute b(j)
    otherlabels = setdiff(ulabel, label[i])
    dics = rep(0,k-1)
    for (j in 1:(k-1)){
      jlabel = which(ulabel==otherlabels[j])
      idx.jlabel = memberlist[[jlabel]]
      dics[j] = base::mean(as.vector(distmat[i,idx.jlabel]))
    }
    vec.b[i] = min(dics)
  }
  # now compute using vectorized operations
  vec.s = ((vec.b-vec.a)/(base::pmax(vec.a, vec.b)))
  # return
  return(mean(vec.s))
}
