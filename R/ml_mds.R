#' Classical Multidimensional Scaling
#' 
#' 
#' @export
mds <- function(x, k=2){
  ##################################################3
  # Check Input and Transform
  x  = check_input(x)
  D2 = (x^2)          # now squared matrix
  n  = nrow(D2)
  
  ##################################################3
  # Computation
  J = diag(n) - (1/n)*outer(rep(1,n),rep(1,n))
  B = -0.5*J%*%D2%*%J
  eigB = eigen(B)
  
  LL = eigB$values[1:k]
  EE = eigB$vectors[,1:k]
  
  Y  = EE%*%diag(sqrt(LL))
  return(Y)
}