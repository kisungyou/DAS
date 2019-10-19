#' Test for Equality of Two Distributions by Blumberg, Bhaumik, and Walker (2018)
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
#' test2.2018bbw(dmat, lab)                    # run the code !
#' 
#' \dontrun{
#' ## WARNING: computationally heavy. 
#' #  Let's compute empirical Type 1 error at alpha=0.05
#' niter = 496  
#' pvals = rep(0,niter)
#' for (i in 1:niter){
#'   dat = matrix(rnorm(300),ncol=3)
#'   lab = c(rep(1,60), rep(2,40))
#'   pvals[i] = test2.2018bbw(as.matrix(dist(dat)), lab)$p.value
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
#' abline(0,1, lwd=2, col="red")
#' }
#' 
#' @references 
#' \insertRef{blumberg_testing_2018}{DAS}
#' 
#' @export
test2.2018bbw <- function(x, label){
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
  # Computation
  # 1. Separate the distance matrices
  id.x = (label==ulabel[1])
  id.y = (label==ulabel[2])
  Mx = x[id.x,id.x]
  My = x[id.y,id.y]
  m  = nrow(Mx)
  n  = nrow(My)
  
  # 2. find the minimal index
  idrowMx = which.min(rowSums(Mx))
  idrowMy = which.min(rowSums(My))
  if (length(idrowMx) > 1){
    idrowMx = base::sample(idrowMx, 1)
  }
  if (length(idrowMy) > 1){
    idrowMy = base::sample(idrowMy, 1)
  }
  
  # 3. select the target 
  tgtx = sort(as.vector(Mx[idrowMx,])); tgtx = tgtx[2:m]
  tgty = sort(as.vector(My[idrowMy,])); tgty = tgty[2:n]
  
  # 4. compute p-value for two-sample kolmogorov smirnov test
  ksout  = stats::ks.test(tgtx, tgty)
  Tmn    = as.double(ksout$statistic)
  pvalue = stats::ks.test(tgtx, tgty)$p.value
  
  ##############################################################
  # REPORT
  thestat = Tmn
  hname   = "Test for Equality of Two Distributions by Blumberg, Bhaumik, and Walker (2018)"
  Ha    = "two distributions are not equal"
  names(thestat) = "Tmn"
  res   = list(statistic=thestat, p.value=pvalue, alternative = Ha, method=hname, data.name = DNAME)
  class(res) = "htest"
  return(res)
}