#' Classical Multidimensional Scaling
#' 
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
#' ## run the algorithm
#' iris.cmds = cmds(dmat, ndim=2)
#' 
#' ## extract coordinates and class information
#' cx = iris.cmds$embed # embedded coordinates of CMDS
#' icol = iris[,5]      # class information
#' 
#' ## visualize
#' par(pty="s")
#' mc = paste("CMDS with STRESS=",round(iris.cmds$stress,4),sep="")
#' plot(cx, col=icol,pch=19,main=mc)
#' 
#' @references 
#' \insertRef{torgerson_multidimensional_1952}{DAS}
#' 
#' @author Kisung You
#' @export
cmds <- function(x, ndim=2){
  ##################################################3
  # Check Input and Transform
  ndim = round(ndim)
  k  = as.integer(ndim)
  x  = check_input(x)
  D2 = (x^2)          # now squared matrix
  n  = nrow(D2)
  if ((length(ndim)>1)||(ndim<1)||(ndim>=nrow(x))){
    stop("* DAS::cmds - 'ndim' should be an integer in [1,nrow(x)). ")
  }
  
  ##################################################3
  # Computation
  J = diag(n) - (1/n)*outer(rep(1,n),rep(1,n))
  B = -0.5*J%*%D2%*%J
  eigB = eigen(B)
  
  LL = eigB$values[1:k]
  EE = eigB$vectors[,1:k]
  
  # Y  = as.matrix(base::scale((EE%*%diag(sqrt(LL))), center=TRUE, scale=FALSE))
  Y  = EE%*%diag(sqrt(LL))
  DY = as.matrix(stats::dist(Y))
  
  output = list()
  output$embed  = Y
  output$stress = compute_stress(x, DY)
  return(output)
}