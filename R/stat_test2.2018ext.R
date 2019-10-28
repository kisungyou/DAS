#' extended : seems like it's working
#' 
#' @examples 
#' ## small test for CRAN submission
#' dat1 <- matrix(rnorm(60, mean= 1), ncol=2) # group 1 : 30 obs of mean  1
#' dat2 <- matrix(rnorm(50, mean=-1), ncol=2) # group 2 : 25 obs of mean -1
#' 
#' dmat <- as.matrix(dist(rbind(dat1, dat2)))  # Euclidean distance matrix
#' lab  <- c(rep(1,30), rep(2,25))             # corresponding label
#' 
#' test2.2018ext(dmat, lab)                    # run the code !
#' 
#' \dontrun{
#' ## WARNING: computationally heavy. 
#' #  Let's compute empirical Type 1 error at alpha=0.05 : should not be on 45 line
#' niter = 496  
#' pvals1 = rep(0,niter) # 2018 method
#' pvals2 = rep(0,niter) # extended
#' for (i in 1:niter){
#'   dat1 = matrix(rnorm(60*3,mean= 1), ncol=3) 
#'   dat2 = matrix(rnorm(40*3,mean=-1), ncol=3)
#'   dat  = rbind(dat1, dat2) 
#'   lab = c(rep(1,60), rep(2,40))
#'   pvals1[i] = test2.2018bbw(as.matrix(dist(dat)), lab)$p.value
#'   pvals2[i] = test2.2018ext(as.matrix(dist(dat)), lab)$p.value
#'   print(paste("iteration ",i," complete..",sep=""))
#' }
#' 
#' #  Visualize the above at multiple significance levels
#' alphas = seq(from=0.001, to=0.999, length.out=100)
#' errors1 = rep(0,100)
#' errors2 = rep(0,100)
#' for (i in 1:100){
#'    errors1[i] = sum(pvals1<=alphas[i])/niter
#'    errors2[i] = sum(pvals2<=alphas[i])/niter
#' }
#' par(mfrow=c(1,2))
#' plot(alphas, errors1, "b", main="Type 1 : 2018 method", 
#'      xlab="alpha", ylab="error", lwd=2)
#' abline(0,1, lwd=2, col="red")
#' plot(alphas, errors2, "b", main="Type 1 : extension", 
#'      xlab="alpha", ylab="error", lwd=2)
#' abline(0,1, lwd=2, col="red")
#' }
#' 
#' @export
test2.2018ext <- function(x, label){
  ##################################################
  # Check Input and Transform
  DNAME = deparse(substitute(x))
  x  = check_input(x)
  if (length(label)!=nrow(x)){
    stop("* test2.2018bbw : length of 'label' vector should equal to the number of observations.")
  }
  ulabel = unique(label)
  if (length(ulabel)!=2){ # also, we can take only 2-sample data
    stop("* test2.2018bbw : this method only supports the data with two classes.")
  }
  
  ##################################################
  # Computation : idea, from the global minimum, test two
  # 1. Separate the distance matrices
  id.x = (label==ulabel[1])
  id.y = (label==ulabel[2])
  m = length(id.x)
  n = length(id.y)
  
  idmin = which.min(rowSums(x^2))[1]
  tgtx  = sort(as.vector(x[idmin,id.x]))
  tgty  = sort(as.vector(x[idmin,id.y]))
  
  # 3. select the target 
  tgtx = tgtx[2:m]
  tgty = tgty[2:n]
  
  # 4. compute p-value for two-sample kolmogorov smirnov test
  ksout  = stats::ks.test(tgtx, tgty)
  Tmn    = as.double(ksout$statistic)
  pvalue = stats::ks.test(tgtx, tgty)$p.value
  
  ##############################################################
  # REPORT
  thestat = Tmn
  hname   = "Extension to  Blumberg, Bhaumik, and Walker (2018)"
  Ha    = "two distributions are not equal"
  names(thestat) = "Tmn"
  res   = list(statistic=thestat, p.value=pvalue, alternative = Ha, method=hname, data.name = DNAME)
  class(res) = "htest"
  return(res)
}