#' Hausdorff Distance
#' 
#' @examples 
#' \dontrun{
#' ## create two datasets from bivariate normal
#' ## let's try to see the evolution of Hausdorff distance
#' nmax = 1000
#' X = matrix(rnorm(nmax*2),ncol=2) # obs. for X
#' Y = matrix(rnorm(nmax*2),ncol=2) # obs. for Y
#' 
#' ## compute cross-distance between X and Y
#' dXY = array(0,c(nmax,nmax))
#' for (i in 1:nmax){
#'   vx = as.vector(X[i,])
#'   for (j in 1:nmax){
#'     vy  = as.vector(Y[j,])
#'     dXY[i,j] = sqrt(sum((vx-vy)^2))
#'   }
#' }
#' 
#' ## compute
#' ndraw = 500
#' xgrid = 2:ndraw
#' ygrid = rep(0,ndraw-1)
#' for (i in 1:(ndraw-1)){
#'   ytmps = rep(0,10)
#'   for (j in 1:10){
#'     id1 = base::sample(1:nmax, i+1)
#'     id2 = base::sample(1:nmax, i+1)
#'     pXY = dXY[id1,id2]
#'     ytmps[j] = hausdorff(pXY)$distance
#'   }
#'   ygrid[i] = base::mean(ytmps)
#'   print(paste("Iteration ",i+1,"/",ndraw," Complete..",sep=""))
#' }
#' 
#' ## visualize
#' opar <- par(pty="s")
#' plot(xgrid, ygrid, "b", lwd=1, main="Evolution of Average Hausdorff Distance from 10 Runs",
#'      xlab="number of samples", ylab="distance", pch=18)
#' par(opar)
#' }
#' 
#' @export
hausdorff <- function(dxy){
  ##################################################
  # Preprocessing
  if ((!is.matrix(dxy))||(any(dxy<0))||(any(is.na(dxy)))||(any(is.infinite(dxy)))){
    stop("* hausdorff : input 'dxy' should be a matrix of nonnegative real numbers.")
  }
  
  ##################################################
  # Main Computation
  val1 = hausdorff.single(dxy)
  val2 = hausdorff.single(t(dxy))
  
  ##################################################
  # Computation
  output = list()
  output$distance = base::max(c(val1,val2))
  return(output)
}


#   -----------------------------------------------------------------------
#  https://en.wikipedia.org/wiki/Hausdorff_distance
#' @keywords internal
#' @noRd
hausdorff.single <- function(x){
  return(base::max(base::apply(x, 2, base::min)))
}
