#' Biswas and Ghosh (2014)
#' 
#' 
#' @examples 
#' ## small test for CRAN submission
#' dat1 <- matrix(rnorm(60, mean= 1), ncol=2) # group 1 : 30 obs of mean  1
#' dat2 <- matrix(rnorm(50, mean=-1), ncol=2) # group 2 : 25 obs of mean -1
#' 
#' dmat <- as.matrix(dist(rbind(dat1, dat2)))  # Euclidean distance matrix
#' lab  <- c(rep(1,30), rep(2,25))             # corresponding label
#' 
#' test2.2014bg(dmat, lab)                         # run the code !
#' 
#' \dontrun{
#' ## WARNING: computationally heavy. 
#' #  Let's compute empirical Type 1 error at alpha=0.05
#' niter = 496  
#' pvals = rep(0,niter)
#' for (i in 1:niter){
#'   dat = matrix(rnorm(200),ncol=2)
#'   lab = c(rep(1,50), rep(2,50))
#'   pvals[i] = test2.2014bg(as.matrix(dist(dat)), lab)$p.value
#'   print(paste("iteration ",i," complete..",sep=""))
#' }
#' print(paste("* Empirical Type 1 Error : ",sum(pvals<=0.05)/niter,sep=""))
#' 
#' #  Visualize the above at multiple significance levels
#' alphas = seq(from=0.001, to=0.999, length.out=100)
#' errors = rep(0,100)
#' for (i in 1:100){
#'    errors[i] = sum(pvals<=alphas[i])/niter
#' }
#' plot(alphas, errors, "b", main="Empirical Type 1 Errors", 
#'      xlab="alpha", ylab="error", lwd=2)
#' abline(v=0.05, lwd=2, col="red")
#' }
#' 
#' @export
test2.2014bg <- function(x, label, mc.iter=999){
  ##################################################
  # Check Input and Transform
  DNAME = deparse(substitute(x))
  x  = check_input(x)
  if (length(label)!=nrow(x)){
    stop("* test2.2014bg : length of 'label' vector should equal to the number of observations.")
  }
  if (length(unique(label))!=2){ # also, we can take only 2-sample data
    stop("* test2.2014bg : this method only supports the data with two classes.")
  }
  
  ##################################################
  # Rearrange data and set as Biswas & Ghosh code from SHT
  x2 = arrange_2sample(x, label)
  m   = x2$m
  n   = x2$n
  DXY = x2$mat
  
  ##################################################
  # H0 computation
  DX0 = DXY[1:m,1:m]
  DY0 = DXY[(m+1):(m+n),(m+1):(m+n)]
  DZ0 = DXY[1:m,(m+1):(m+n)]
  Tmn = R_eqdist_2014BG_statistic(DX0,DY0,DZ0)
  
  ##############################################################
  # Monte Carlo Computation
  nreps = mc.iter
  Tvec  = rep(0,nreps)
  for (i in 1:nreps){
    idx = sample(1:(m+n), m, replace=FALSE)
    idy = setdiff(1:(m+n), idx)
    
    DX1 = DXY[idx,idx]
    DY1 = DXY[idy,idy]
    DZ1 = DXY[idx,idy]
    Tvec[i] = R_eqdist_2014BG_statistic(DX1,DY1,DZ1)
  }
  pvalue = (sum(Tvec>=Tmn)+1)/(nreps+1)
 
  
  ##############################################################
  # REPORT
  thestat = Tmn
  hname   = "Test for Equality of Two Distributions by Biswas and Ghosh (2014)"
  Ha    = "two distributions are not equal"
  names(thestat) = "Tmn"
  res   = list(statistic=thestat, p.value=pvalue, alternative = Ha, method=hname, data.name = DNAME)
  class(res) = "htest"
  return(res)
}



# auxiliary functions -----------------------------------------------------
#' @keywords internal
#' @noRd
R_eqdist_2014BG_statistic <- function(DX,DY,DXY){
  m = nrow(DXY)
  n = ncol(DXY)
  
  muff = sum(DX[upper.tri(DX)])/(m*(m-1)/2)
  mufg = sum(DXY)/(m*n)
  mugg = sum(DY[upper.tri(DY)])/(n*(n-1)/2)
  
  vec1 = c(muff,mufg)
  vec2 = c(mufg,mugg)
  output = sum((vec1-vec2)^2)
  return(output)
}
