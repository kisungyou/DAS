#' Approximate K-Means Clustering
#' 
#' 
#' @export
akmeans <- function(x, k=2, init=c("Random","Medoids"), maxiter=20, label.init=NULL){
  ##################################################3
  # Check Input and Transform
  x  = check_input(x)
  n  = nrow(x)
  K  = round(k)
  if (K >= n){
    stop("* akmeans : 'k' should be smaller than the number of observations.")
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
  nmin = 3                          # for test, I use 3-point
  for (iter in 1:iter.max){
    # nmin-point bounding polytope
    voronoilist = list()
    for (i in 1:K){
      idnow = which(labels.old==i)
      if (length(idnow) < nmin){    # for a singleton or two-point cluster
        voronoilist[[i]] = idnow
      } else {
        voronoilist[[i]] = idnow[order(base::rowSums(x[idnow,idnow]))[1:nmin]]
      }
    }
    # need to compute pseudodistance
    pseudodist = array(0,c(n,K))
    for (k in 1:K){
      idnow = voronoilist[[k]]
      for (i in 1:n){
        idnew = c(i,idnow)
        if (length(idnow)==1){
          pseudodist[i,k] = as.double(x[idnow,idnow])
        } else {
          pseudodist[i,k] = aux_pseudomean(x[idnew,idnew])  
        }
      }
    }
    
    # update labels
    labels.new = rep(0,n)
    for (i in 1:n){
      tgtvec = as.vector(pseudodist[i,])
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
# # test 1. compare naive kmeans and approximate k-means clustering
# cl.kmeans  = kmeans(X,4)$cluster
# cl.akclust = akmeans(as.matrix(dist(X)), k=4)$cluster
# par(mfrow=c(1,2), pty="s")
# plot(X[,1],X[,2],col=cl.kmeans, main="naive k-means")
# plot(X[,1],X[,2],col=cl.akclust,main="approximate k-means")
# 
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
# cl.akcl2    = akmeans(as.dist(dpath),k=3,init="random")$cluster
# cl.ksets    = ksets(as.dist(dpath),k=3, init="random")$cluster
# cl.kmedgr   = cluster::pam(as.dist(dpath),k=3)$clustering
# 
# 
# par(mfrow=c(1,3),pty="s")
# plot(X[,1],X[,2],col=cl.akcl2, pch=19, main="graph approximate k-means", cex=0.25)
# plot(X[,1],X[,2],col=cl.ksets, pch=19, main="graph k-sets", cex=0.25)
# plot(X[,1],X[,2],col=cl.kmedgr, pch=19, main="graph k-medoids", cex=0.25)
