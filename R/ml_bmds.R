#' Bayesian Multidimensional Scaling by Oh and Raftery (2001)
#'  
#' @examples 
#' ## use simple example of iris dataset with perturbation
#' data(iris) 
#' dmat = as.matrix(stats::dist(iris[,1:4]))
#' 
#' ## run Bayesian MDS
#' #  let's run 49 iterations (CRAN) quietly
#' iris.cmds =  mds(dmat, ndim=2)
#' iris.bmds = bmds(dmat, ndim=2, mc.iter=49, par.step=(2.38^2)) 
#' 
#' ## extract coordinates and class information
#' cx = iris.cmds$embed # embedded coordinates of CMDS
#' bx = iris.bmds$embed #                         BMDS
#' icol = iris[,5]      # class information
#' 
#' ## visualize
#' par(mfrow=c(1,2),pty="s")
#' mc = paste("CMDS with STRESS=",round(iris.cmds$stress,4),sep="")
#' mb = paste("BMDS with STRESS=",round(iris.bmds$stress,4),sep="")
#' plot(cx, col=icol,pch=19,main=mc)
#' plot(bx, col=icol,pch=19,main=mb)
#' 
#' @export
bmds <- function(x, ndim=2, par.a=5, par.alpha=0.5, par.step=1, mc.iter=8128, verbose=FALSE){
  ######################################################
  # Initialization
  x    = check_input(x)
  ndim = round(ndim)

  if ((length(ndim)>1)||(ndim<1)||(ndim>=nrow(x))){
    stop("* DAS::mds - 'ndim' should be an integer in [1,nrow(x)). ")
  }
  
  n = nrow(x)
  m = n*(n-1)/2
  
  ######################################################
  # Preliminary Computation
  # 1. apply CMDS for initialization
  y     = as.matrix(base::scale(mds(x, ndim)$embed, # (N x ndim) centered 
                                center=TRUE, scale=FALSE)) 
  Delta = as.matrix(stats::dist(y))           # (N x N) pairwise distances
  
  # 2. initialization
  eigy   = base::eigen(cov(y))       
  X0     = y%*%eigy$vectors    # (N x ndim) rotated
  gamma0 = diag(X0)                 # variances ?
  sigg0  = compute_SSR(x, Delta)/m; # 
  beta0  = apply(X0,2,var)/2
  
  # 3. run the main part
  runcpp <- main_bmds(x, X0, sigg0, par.a, par.alpha, mc.iter, par.step, verbose, beta0)
  Xsol   <- runcpp$solX
  Xdist  <- as.matrix(stats::dist(Xsol))
  
  output = list()
  output$embed  = Xsol
  output$stress = compute_stress(x, Xdist)
  return(output)
  # return Rcpp::List::create(Rcpp::Named("solX")=Xsol,Rcpp::Named("solSSR")=SSRsol);
}