#' K-Sets Clustering Algorithm
#' 
#' sort of autotuning is happening..
#' what about.. just.. recursive partitioning?
#' 
#' @references 
#' \insertRef{chang_mathematical_2016}{DAS}
#' 
#' @export
ksets <- function(x, k=2, init=c("Random","Medoids"), maxiter=20, label.init=NULL){
  ##################################################3
  # Check Input and Transform
  x  = check_input(x)
  n  = nrow(x)
  K  = round(k)
  if (K >= n){
    stop("* ksets : 'k' should be smaller than the number of observations.")
  }
  allinit  = c("random","medoids")
  init     = base::match.arg(tolower(init), allinit)
  iter.max = round(maxiter)
  
  # lesson from example.. if not a connected graph.. make an arbitrary large embedding
  x[is.infinite(x)] = max(x[!is.infinite(x)])*10
  
  ##################################################3
  # initialize
  if ((length(label.init)==0)&&(is.null(label.init))){
    if (all(init=="random")){
      labels = c(1:K, base::sample(1:K, n-K, replace=TRUE))
      labels.old = base::sample(labels)
    } else if (all(init=="medoids")){
      labels.old = as.vector(cluster::pam(stats::as.dist(x), k=K)$clustering)
    }
  } else {
    labels.old = as.integer(as.factor(label.init))
  }
  
  
  
  ##################################################
  # iterate
  history = list()
  for (iter in 1:iter.max){
    
    history[[iter]] = labels.old
    
    # membership list
    membership = rclust_membership(labels.old)
    # let's compute triangular distances : (n x k)
    tridist = compute_tridist(x, membership)
    # update labels
    labels.new = rep(0,n)
    for (i in 1:n){
      tgtvec = as.vector(tridist[i,])
      idmin  = which.min(tgtvec)
      if (length(idmin) > 1){
        labels.new[i] = sample(idmin, 1)
      } else {
        labels.new[i] = idmin
      }
    }
    # difference
    labels.new = as.integer(as.factor(labels.new))
    diffs = as.double(mclustcomp::mclustcomp(labels.old, labels.new, types="rand")[2])
    if (diffs > 0.99){
      break
    }
    labels.old = labels.new
    print(paste(" iteration ",iter," complete..",sep=""))
  }
  
  ##################################################
  # return
  output = list()
  output$cluster = labels.old
  output$history = history
  return(output)
}


#   -----------------------------------------------------------------------
#' @keywords internal
#' @noRd
compute_tridist <- function(x, membership){
  # paramters & variables
  N = nrow(x)
  K = length(membership)
  output = array(0,c(N,K))

  for (k in 1:K){
    idk = membership[[k]]
    S   = base::length(idk)
    ccS = sum(x[idk,idk])/(S^2)
    for (n in 1:N){
      output[n,k] = (sum(x[n,idk])*2/S) - ccS
    }
  }
  # return
  return(output)
}


# # personal test
# library(maotai)
# library(cluster)
# npt = 300
# inner.r = runif(npt, min=10, max=12); inner.theta = runif(npt, min=0, max=2*pi)
# outer.r = runif(npt, min=20, max=22); outer.theta = runif(npt, min=0, max=2*pi)
# X1 = cbind(inner.r*cos(inner.theta), inner.r*sin(inner.theta))
# X2 = cbind(outer.r*cos(outer.theta), outer.r*sin(outer.theta))
# X  = rbind(X1,X2)
# 
# cl.kmeans = kmeans(X,4)$cluster
# # exactly follow the protocol from the paper
# dX = as.matrix(dist(X))
# ddX = array(0,c(nrow(dX),nrow(dX)))
# for (i in 1:(nrow(dX)-1)){
#   for (j in (i+1):nrow(dX)){
#     if (dX[i,j] < 5){
#       ddX[i,j] <- ddX[j,i] <- 1
#     }
#   }
# }
# dpath       = maotai::shortestpath(ddX)
# cl.kmedoids = cluster::pam(as.dist(dX),k=4)$clustering
# cl.ksets    = ksets(as.dist(dpath),k=4, init="random")$cluster
# cl.kmedgr   = cluster::pam(as.dist(dpath),k=4)$clustering
# 
# 
# par(mfrow=c(2,2),pty="s")
# plot(X[,1],X[,2],col=cl.kmeans, pch=19, main="naive k-means", cex=0.25)
# plot(X[,1],X[,2],col=cl.kmedoids, pch=19, main="naive k-medoids", cex=0.25)
# plot(X[,1],X[,2],col=cl.ksets, pch=19, main="graph k-sets", cex=0.25)
# plot(X[,1],X[,2],col=cl.kmedgr, pch=19, main="graph k-medoids", cex=0.25)
# 
# kset's history
# ksethist = ksets(dpath,k=4,init="random")
# nhist = length(ksethist$history)
# imrow = floor(sqrt(nhist))
# imcol = ceiling(nhist/imrow)
# par(mfrow=c(imrow,imcol),pty="s")
# for (i in 1:nhist){
#   pm = paste("nclust=",length(unique(ksethist$history[[i]])),sep="")
#   plot(X[,1],X[,2],col=ksethist$history[[i]], main=pm, cex=0.4, pch=19)
# }
