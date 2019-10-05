#' optimal transport; returns p-Wasserstein distance
#' 
#' c(x,y) = dxy^p; ground cost/metric
#' 
#' @examples 
#' ## create two small datasets from bivariate normal
#' X = matrix(rnorm(5*2),ncol=2) # 5 obs. for X
#' Y = matrix(rnorm(5*2),ncol=2) # 5 obs. for Y
#' 
#' ## compute cross-distance between X and Y
#' dXY = array(0,c(5,5))
#' for (i in 1:5){
#'   vx = as.vector(X[i,])
#'   for (j in 1:5){
#'     vy  = as.vector(Y[j,])
#'     dXY[i,j] = sqrt(sum((vx-vy)^2))
#'   }
#' }
#' 
#' ## compute the distance and report
#' output = coupling(dXY, p=2) # 2-Wasserstein distance
#' image(output$coupling, main=paste("distance=",round(output$distance,4),sep=""))
#' 
#' \dontrun{
#' ## create two datasets from bivariate normal
#' ## let's try to see the evolution of 2-Wasserstein distance
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
#' xgrid = 2:nmax
#' ygrid = rep(0,nmax-1)
#' for (i in 1:(nmax-1)){
#'   pXY = dXY[1:(i+1),1:(i+1)]
#'   ygrid[i] = coupling(pXY, p=2)$distance
#'   print(paste("Iteration ",i+1,"/",nmax," Complete..",sep=""))
#' }
#' 
#' ## visualize
#' plot(xgrid, ygrid, "b", lwd=1, main="Evolution of 2-Wasserstein Distances",
#'      xlab="number of samples", ylab="distance", pch=18)
#' }
#' 
#' @export
coupling <- function(dxy, p=1, wx, wy,
                        method = c("networkflow", "shortsimplex", "revsimplex", "primaldual")){
  ##################################################
  # Preprocessing
  # 1. dxy
  if ((!is.matrix(dxy))||(any(dxy<0))||(any(is.na(dxy)))||(any(is.infinite(dxy)))){
    stop("* coupling : input 'dxy' should be a matrix of nonnegative real numbers.")
  }
  n = nrow(dxy)
  m = ncol(dxy)
  # 2. p; pass
  if ((length(p)>1)||(p<=0)){
    stop("* coupling : 'p' should be a nonnegative real number of Inf.")
  }
  # 3. wx and wy
  if (missing(wx)){wx = rep(1/n,n)}; wx=wx/sum(wx)
  if (missing(wy)){wy = rep(1/m,m)}; wy=wy/sum(wy)
  if ((!is.vector(wx))||(!check_Sto1(wx))||(length(wx)!=n)||(any(wx<0))){
    stop("* coupling : 'wx' should be a vector of nonnegative numbers summing to 1.")
  }
  if ((!is.vector(wy))||(!check_Sto1(wx))||(length(wy)!=m)||(any(wy<0))){
    stop("* coupling : 'wy' should be a vector of nonnegative numbers summing to 1.")
  }
  # 4. method
  mymethod = match.arg(method)
  
  ##################################################
  # Main Computation
  if (is.infinite(p)){ # p=Inf
    output = list()
    output$distance = max(dxy)
  } else {
    cxy = dxy^p
    output = compute.coupling(cxy, transport::transport(a=wx, b=wy, costm = cxy))
    output$distance = (output$distance^(1/p))
    return(output)
  }
}



# auxiliary functions -----------------------------------------------------
#' @keywords internal
#' @noRd
compute.coupling <- function(costm, plan){
  nx = nrow(costm)
  ny = ncol(costm)
  output = array(0,c(nx,ny))
  
  nplans = nrow(plan)
  for (i in 1:nplans){
    id.from = plan[i,1]
    id.to   = plan[i,2]
    vals    = plan[i,3]
    
    output[id.from,id.to] = vals
  }
  
  result = list()
  result$distance = sum(output*costm)
  result$coupling = output
  return(result)
}

#' @keywords internal
#' @noRd
check_Sto1 <- function(w){
  return((abs(sum(w)-1) < 100*.Machine$double.eps))
}