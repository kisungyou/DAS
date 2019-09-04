#' sum of squared distances
#' 
#' 
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
#' testk.ssd(dmat, lab)                        # run the code !
#' 
#' \dontrun{
#' ## WARNING: computationally heavy. 
#' #  Let's compute empirical Type 1 error at alpha=0.05
#' niter = 496  
#' pvals = rep(0,niter)
#' for (i in 1:niter){
#'   dat = matrix(rnorm(200),ncol=2)
#'   lab = c(rep(1,50), rep(2,50))
#'   pvals[i] = testk.ssd(as.matrix(dist(dat)), lab)$p.value
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
testk.ssd <- function(x, label, mc.iter=999){
  ##################################################3
  # Check Input and Transform
  DNAME = deparse(substitute(x))
  x  = check_input(x)
  D2 = (x^2)          # now squared matrix
  N  = nrow(D2)
  if (length(label)!=N){
    stop("* testk.ssd : length of 'label' vector should equal to the number of observations.")
  }
  
  ##################################################
  # Preliminary items
  idlist = label2idlist(label)            # id of observations per group
  vec.ng = unlist(lapply(idlist, length)) # number of observations per group
  G      = length(vec.ng)
  if (G < 2){
    stop("* testk.ssd : given label says we only have one group. Nothing to test.")
  }
  if (any(vec.ng < 2)){
    stop("* testk.ssd : we have a group of size 1. This case is not valid for distribution test.")
  }
  thestat = testk.ssd.statistic(D2, idlist)
  
  ##################################################
  # Permutations
  permvec  = rep(0, mc.iter)
  for (i in 1:mc.iter){
    permidlist = sample_group(vec.ng)
    permvec[i] = testk.ssd.statistic(D2, permidlist)
  }
  
  ##############################################################
  # REPORT
  pvalue  = (sum(permvec <= thestat)+1)/(mc.iter+1)
  hname   = "Test for Equality of Distributions with Sum of Squared Distances."
  Ha      = "one of equalities does not hold."
  names(thestat) = "SSD"
  res   = list(statistic=thestat, p.value=pvalue, alternative = Ha, method=hname, data.name = DNAME, permvec=permvec)
  class(res) = "htest"
  return(res)
}


# auxiliary functions -----------------------------------------------------
#' @keywords internal
#' @noRd
testk.ssd.statistic <- function(D2, idlist){
  vec.nm = unlist(lapply(idlist, length))
  statistic = 0
  for (i in 1:length(vec.nm)){
    idnow = idlist[[i]]  # current id
    nm    = vec.nm[i]    # number of current id
    if (nm > 1){
      statistic = statistic + sum(D2[idnow,idnow])/(2*nm*(nm-1))
    }
  }
  return(statistic)
}
