#' MDF test
#' 
#' 
#' Given the distance matrix of two or more samples, it tests
#' \deqn{H_0 : F_{g_1} = \cdots = F_{g_K} F_Y\quad vs\quad H_1 : \textrm{ not } H_0}
#' using the procedure by Minas et al. (2011) in a nonparametric way 
#' that depends on \code{mc.iter} number of permutations. 
#' 
#' @param x an \eqn{(N\times N)} matrix or \code{dist} object.
#' @param label length-\eqn{N} vector of class/grouping label.
#' @param mc.iter the number of permutations to be taken.
#'
#' @return a (list) object of \code{S3} class \code{htest} containing: \describe{
#' \item{statistic}{a test statistic.}
#' \item{p.value}{\eqn{p}-value under \eqn{H_0}.}
#' \item{alternative}{alternative hypothesis.}
#' \item{method}{name of the test.}
#' \item{data.name}{name(s) of provided sample data.}
#' \item{permvec}{statistics from each permutation.}
#' }
#' 
#' @examples 
#' ## small test for CRAN submission
#' dat1 <- matrix(rnorm(60, mean= 1), ncol=2) # group 1 : 30 obs of mean  1
#' dat2 <- matrix(rnorm(50, mean=-1), ncol=2) # group 2 : 25 obs of mean -1
#' 
#' dmat <- as.matrix(dist(rbind(dat1, dat2)))  # Euclidean distance matrix
#' lab  <- c(rep(1,30), rep(2,25))             # corresponding label
#' 
#' testk.mdf(dmat, lab)                        # run the code !
#' 
#' \dontrun{
#' ## WARNING: computationally heavy. 
#' #  Let's compute empirical Type 1 error at alpha=0.05
#' niter = 496  
#' pvals = rep(0,niter)
#' for (i in 1:niter){
#'   dat = matrix(rnorm(200),ncol=2)
#'   lab = c(rep(1,50), rep(2,50))
#'   pvals[i] = testk.mdf(as.matrix(dist(dat)), lab)$p.value
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
#' @references 
#' \insertRef{minas_distance-based_2011}{DAS}
#' 
#' @export
testk.mdf <- function(x, label, mc.iter=999){
  ##################################################3
  # Check Input and Transform
  DNAME = deparse(substitute(x))
  x  = check_input(x)
  D2 = (x^2)          # now squared matrix
  N  = nrow(D2)
  if (length(label)!=N){
    stop("* test.mdf : length of 'label' vector should equal to the number of observations.")
  }
  
  ##################################################
  # Preliminary items
  idlist = label2idlist(label)            # id of observations per group
  vec.ng = unlist(lapply(idlist, length)) # number of observations per group
  G      = length(vec.ng)
  if (G < 2){
    
    stop("* test.mdf : given label says we only have one group. Nothing to test.")
  }
  J      = (1/N)*outer(rep(1,N),rep(1,N)) # (1/N)*ones(N,N)
  
  Tdel = sum(D2)/(2*N)   # always fixed : total variation
  Hc   = array(0,c(N,N)) # always fixed : something like scaler
  start = 1
  for (g in 1:G){
    ng = vec.ng[g]
    Hc[start:(start+ng-1),start:(start+ng-1)] = (1/ng)*outer(rep(1,ng),rep(1,ng))
    start = start + ng
  }
  Hc   = Hc - J
  InJ  = diag(N)-J
  
  Gdel     = (InJ%*%(-0.5*D2)%*%InJ) # current statistic
  trHcGdel = sum(diag(Hc%*%Gdel))
  thestat  = trHcGdel/(Tdel - trHcGdel)
  
  ##################################################
  # Permutations
  permvec  = rep(0, mc.iter)
  for (i in 1:mc.iter){
    permid = unlist(sample_group(vec.ng))
    tmpD2  = D2[permid,permid]
    Gdel     = (InJ%*%(-0.5*tmpD2)%*%InJ) # current statistic
    trHcGdel = sum(diag(Hc%*%Gdel))
    permvec[i] = trHcGdel/(Tdel - trHcGdel)
  }
  
  ##############################################################
  # REPORT
  pvalue  = (sum(permvec >= thestat)+1)/(mc.iter+1)
  hname   = "Test for Equality of Distributions by Minas et al. (2011)"
  Ha      = "one of equalities does not hold."
  names(thestat) = "F"
  res   = list(statistic=thestat, p.value=pvalue, alternative = Ha, method=hname, data.name = DNAME, permvec=permvec)
  class(res) = "htest"
  return(res)
}
