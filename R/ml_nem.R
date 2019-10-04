#' Negative Eigenvalue Magnitude
#' 
#' NEM == 0; data is euclidean
#' NEM >> 0; far
#' 
#' Gram matrix -1/2 JD^2J be positive semidefinite == EUclidean
#' 
#' @examples 
#' ## use simple example of iris dataset 
#' data(iris) 
#' dat  = iris[,1:4]
#' dmat = as.matrix(stats::dist(dat))
#' 
#' ## run the algorithm
#' nem(dmat)
#' 
#' @references 
#' \insertRef{pekalska_non-euclidean_2006}{DAS}
#' 
#' @seealso \code{\link{nef}}
#' @author Kisung You
#' @export
nem <- function(x){
  ##################################################3
  # Check Input and Transform
  x  = check_input(x)
  D2 = (x^2)
  n  = nrow(D2)
  
  ##################################################3
  # Computation
  H = diag(n) - (1/n)*base::outer(rep(1,n),rep(1,n))
  S = -0.5*(H%*%D2%*%H)
  eigS = base::eigen(S)
  evals = eigS$values
  
  ##################################################3
  # Finalize
  output = abs(min(evals))/max(evals)
  return(output)
}