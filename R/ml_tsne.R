#' t-SNE embedding
#' 
#' 
#' 
#' @examples 
#' ## use simple example of iris dataset 
#' data(iris) 
#' dat  = iris[,1:4]
#' dat.n = nrow(dat)
#' dat.p = ncol(dat)
#' dat  = dat + matrix(rnorm(dat.n*dat.p, sd=0.1), ncol=dat.p)
#' dmat = as.matrix(stats::dist(dat)) # distance matrix
#' 
#' ## run t-SNE and MDS for comparison
#' iris.cmds = cmds(dmat, ndim=2)
#' iris.tsne = tsne(dmat, ndim=2)
#' 
#' ## extract coordinates and class information
#' cx = iris.cmds$embed # embedded coordinates of CMDS
#' tx = iris.tsne$embed #                         t-SNE
#' icol = iris[,5]      # class information
#' 
#' 
#' ## visualize
#' mc = paste("CMDS with STRESS=",round(iris.cmds$stress,4),sep="")
#' mt = paste("tSNE with STRESS=",round(iris.tsne$stress,4),sep="")
#' 
#' opar <- par(mfrow=c(1,2), pty="s")
#' plot(cx, col=icol,pch=19,main=mc)
#' plot(tx, col=icol,pch=19,main=mt)
#' par(opar)
#' 
#' @export
tsne <- function(x, ndim=2, ...){
  ##################################################3
  # Check Input and Transform
  k  = as.integer(ndim)
  x  = check_input(x)
  if ((length(ndim)>1)||(ndim<1)||(ndim>=nrow(x))){
    stop("* DAS::cmds - 'ndim' should be an integer in [1,nrow(x)). ")
  }
  
  ##################################################
  # Pass to 'Rtsne'
  dx  = stats::as.dist(x)
  tmpout = Rtsne::Rtsne(dx, dims=k, ..., is_distance=TRUE)
  Y   = tmpout$Y
  DY  = as.matrix(stats::dist(Y))
  
  ##################################################
  # Return
  output = list()
  output$embed  = Y
  output$stress = compute_stress(x, DY)
  return(output)
  
}